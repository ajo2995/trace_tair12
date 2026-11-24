#!/usr/local/bin/perl
# Usage: fetch_ensembl_canonical.pl REGISTRY_FILE DBNAME OUT_DIR
#
# Pulls from an Ensembl core database what is needed to see which transcript
# of each gene it marks as canonical, and how its gene models are built:
#
#   meta.tsv          every row of the meta table (provenance of the load)
#   transcripts.tsv   one row per transcript: gene, biotypes, coordinates,
#                     whether it is the gene's canonical transcript
#   exons.tsv         one row per exon of each transcript, in rank order
#   translations.tsv  where each translation starts and ends
#   attribs.tsv       how many genes and transcripts carry each attribute
#   genes_without_transcripts.tsv   genes that have no transcript row at all
#
# The connection settings (host, port, user, password) are read from an
# Ensembl registry file that defines $def_host, $def_port, $def_user and
# $def_pass, so no password is typed on a command line or copied elsewhere.
# The session is switched to read-only before any query, and only SELECT
# statements are issued.
use strict;
use warnings;
use DBI;

my ($registry, $dbname, $out_dir) = @ARGV;
defined $out_dir or die "usage: $0 registry_file dbname out_dir\n";

my %conn;
open(my $reg, '<', $registry) or die "cannot read $registry: $!\n";
while (<$reg>) {
  $conn{$1} = $2 if /^\s*my\s+\$def_(user|pass|host|port)\s*=\s*'?([^';]+?)'?\s*;/;
}
close $reg;
for my $key (qw(user pass host port)) {
  defined $conn{$key} or die "$registry does not define \$def_$key\n";
}

my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$conn{host};port=$conn{port}",
                       $conn{user}, $conn{pass}, { RaiseError => 1, PrintError => 0 });
eval { $dbh->do("SET SESSION TRANSACTION READ ONLY"); 1 }
  or warn "could not set the session read-only (old server?); continuing with SELECT statements only\n";
my ($version) = $dbh->selectrow_array("SELECT VERSION()");
print STDERR "connected to $dbname on $conn{host} (server $version)\n";

# dump(FILE, HEADER, SQL): write a query result as a tab-separated file.
sub dump_query {
  my ($file, $header, $sql) = @_;
  my $sth = $dbh->prepare($sql);
  $sth->execute;
  open(my $fh, '>', "$out_dir/$file") or die "cannot write $out_dir/$file: $!\n";
  print $fh join("\t", @$header), "\n";
  my $rows = 0;
  while (my @row = $sth->fetchrow_array) {
    print $fh join("\t", map { defined $_ ? $_ : 'NULL' } @row), "\n";
    $rows++;
  }
  close $fh;
  print STDERR sprintf("%-18s %d rows\n", $file, $rows);
}

dump_query('meta.tsv', [qw(species_id meta_key meta_value)],
  "SELECT species_id, meta_key, meta_value FROM meta ORDER BY meta_key, meta_id");

dump_query('transcripts.tsv',
  # is_canonical        the gene table names this transcript as the gene's canonical
  # is_canonical_attrib the transcript carries the separate "is_canonical" attribute
  # The two are stored independently, and a script that updates one does not
  # necessarily update the other.
  [qw(gene_id gene_biotype gene_source gene_analysis seq_region strand transcript_id transcript_biotype
      transcript_source transcript_start transcript_end is_canonical translation_id is_canonical_attrib)],
  "SELECT g.stable_id, g.biotype, g.source, a.logic_name, sr.name, g.seq_region_strand,
          t.stable_id, t.biotype, t.source, t.seq_region_start, t.seq_region_end,
          IF(g.canonical_transcript_id = t.transcript_id, 1, 0), tl.stable_id,
          IF(EXISTS (SELECT 1 FROM transcript_attrib ta
                       JOIN attrib_type at ON at.attrib_type_id = ta.attrib_type_id
                      WHERE ta.transcript_id = t.transcript_id AND at.code = 'is_canonical'), 1, 0)
     FROM gene g
     JOIN transcript t ON t.gene_id = g.gene_id
     JOIN seq_region sr ON sr.seq_region_id = g.seq_region_id
     LEFT JOIN analysis a ON a.analysis_id = g.analysis_id
     LEFT JOIN translation tl ON tl.transcript_id = t.transcript_id
    ORDER BY sr.name, g.seq_region_start, g.stable_id, t.stable_id");

dump_query('exons.tsv', [qw(transcript_id rank exon_dbid exon_start exon_end phase end_phase)],
  "SELECT t.stable_id, et.rank, e.exon_id, e.seq_region_start, e.seq_region_end, e.phase, e.end_phase
     FROM transcript t
     JOIN exon_transcript et ON et.transcript_id = t.transcript_id
     JOIN exon e ON e.exon_id = et.exon_id
    ORDER BY t.stable_id, et.rank");

dump_query('translations.tsv', [qw(transcript_id translation_id start_exon_dbid seq_start end_exon_dbid seq_end)],
  "SELECT t.stable_id, tl.stable_id, tl.start_exon_id, tl.seq_start, tl.end_exon_id, tl.seq_end
     FROM translation tl JOIN transcript t ON t.transcript_id = tl.transcript_id
    ORDER BY t.stable_id");

dump_query('attribs.tsv', [qw(feature code name count distinct_values example_value)],
  "SELECT 'transcript', at.code, at.name, COUNT(*), COUNT(DISTINCT ta.value), MIN(ta.value)
     FROM transcript_attrib ta JOIN attrib_type at ON at.attrib_type_id = ta.attrib_type_id
    GROUP BY at.attrib_type_id
   UNION ALL
   SELECT 'gene', at.code, at.name, COUNT(*), COUNT(DISTINCT ga.value), MIN(ga.value)
     FROM gene_attrib ga JOIN attrib_type at ON at.attrib_type_id = ga.attrib_type_id
    GROUP BY at.attrib_type_id");

dump_query('seq_regions.tsv', [qw(seq_region length coord_system)],
  "SELECT sr.name, sr.length, cs.name
     FROM seq_region sr JOIN coord_system cs ON cs.coord_system_id = sr.coord_system_id
    WHERE cs.name = 'chromosome' OR sr.name IN ('Mt', 'Pt')
    ORDER BY sr.name");

dump_query('genes_without_transcripts.tsv', [qw(gene_id gene_biotype gene_analysis seq_region gene_start gene_end)],
  "SELECT g.stable_id, g.biotype, a.logic_name, sr.name, g.seq_region_start, g.seq_region_end
     FROM gene g
     LEFT JOIN transcript t ON t.gene_id = g.gene_id
     JOIN seq_region sr ON sr.seq_region_id = g.seq_region_id
     LEFT JOIN analysis a ON a.analysis_id = g.analysis_id
    WHERE t.transcript_id IS NULL
    ORDER BY sr.name, g.seq_region_start");

my ($genes, $no_canonical) = $dbh->selectrow_array(
  "SELECT COUNT(*), SUM(canonical_transcript_id IS NULL OR canonical_transcript_id = 0) FROM gene");
print STDERR "genes: $genes; without a canonical transcript: " . ($no_canonical // 0) . "\n";
$dbh->disconnect;
