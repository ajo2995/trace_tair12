#!/usr/local/bin/perl
# Usage: inspect_ensembl_organelles.pl REGISTRY_FILE DBNAME OUT_DIR
#
# Read-only look at the organelle sequence regions (named Mt and Pt) of an
# Ensembl core database: what each region is called, how long it is, which
# attributes and synonyms hang off it, how it is assembled, how many features
# of each kind sit on it, and the sequence itself (written as FASTA).
#
# Everything in the core schema points at a region by its numeric
# seq_region_id, so this shows what would and would not follow if the two
# names were swapped.
#
# Connection settings are read from the registry file ($def_host, $def_port,
# $def_user, $def_pass). The session is set read-only and only SELECT and
# SHOW statements are issued.
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
defined $conn{$_} or die "$registry does not define \$def_$_\n" for qw(user pass host port);

my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$conn{host};port=$conn{port}",
                       $conn{user}, $conn{pass}, { RaiseError => 1, PrintError => 0 });
eval { $dbh->do("SET SESSION TRANSACTION READ ONLY"); 1 } or warn "could not set the session read-only\n";

sub show {
  my ($title, $sql, @bind) = @_;
  print "== $title ==\n";
  my $sth = $dbh->prepare($sql);
  $sth->execute(@bind);
  print join("\t", @{ $sth->{NAME} }), "\n";
  my $rows = 0;
  while (my @row = $sth->fetchrow_array) {
    print join("\t", map { defined $_ ? $_ : 'NULL' } @row), "\n";
    $rows++;
  }
  print "(no rows)\n" unless $rows;
  print "\n";
}

show("sequence regions named Mt or Pt",
  "SELECT sr.seq_region_id, sr.name, sr.length, cs.name AS coord_system, cs.version, cs.rank
     FROM seq_region sr JOIN coord_system cs ON cs.coord_system_id = sr.coord_system_id
    WHERE sr.name IN ('Mt', 'Pt') ORDER BY sr.name");

show("attributes of Mt, Pt and, for comparison, chromosome 1",
  "SELECT sr.name AS region, at.code, sra.value
     FROM seq_region_attrib sra
     JOIN seq_region sr ON sr.seq_region_id = sra.seq_region_id
     JOIN coord_system cs ON cs.coord_system_id = sr.coord_system_id
     JOIN attrib_type at ON at.attrib_type_id = sra.attrib_type_id
    WHERE sr.name IN ('Mt', 'Pt', '1') AND cs.name = 'chromosome'
    ORDER BY sr.name, at.code");

show("synonyms of Mt, Pt and chromosome 1",
  "SELECT sr.name AS region, s.synonym, ed.db_name
     FROM seq_region_synonym s
     JOIN seq_region sr ON sr.seq_region_id = s.seq_region_id
     LEFT JOIN external_db ed ON ed.external_db_id = s.external_db_id
    WHERE sr.name IN ('Mt', 'Pt', '1') ORDER BY sr.name, s.synonym");

show("how Mt and Pt are assembled from lower-level pieces",
  "SELECT asm.name AS region, cs.name AS piece_level, COUNT(*) AS pieces,
          MIN(cmp.name) AS first_piece_name, MAX(cmp.name) AS last_piece_name,
          MIN(a.asm_start) AS from_bp, MAX(a.asm_end) AS to_bp
     FROM assembly a
     JOIN seq_region asm ON asm.seq_region_id = a.asm_seq_region_id
     JOIN seq_region cmp ON cmp.seq_region_id = a.cmp_seq_region_id
     JOIN coord_system cs ON cs.coord_system_id = cmp.coord_system_id
    WHERE asm.name IN ('Mt', 'Pt')
    GROUP BY asm.name, cs.name ORDER BY asm.name");

show("every piece of Mt and Pt, by name",
  "SELECT asm.name AS region, cmp.name AS piece, a.asm_start, a.asm_end
     FROM assembly a
     JOIN seq_region asm ON asm.seq_region_id = a.asm_seq_region_id
     JOIN seq_region cmp ON cmp.seq_region_id = a.cmp_seq_region_id
    WHERE asm.name IN ('Mt', 'Pt') ORDER BY asm.name, a.asm_start");

show("all chromosome-level regions",
  "SELECT sr.seq_region_id, sr.name, sr.length
     FROM seq_region sr JOIN coord_system cs ON cs.coord_system_id = sr.coord_system_id
    WHERE cs.name = 'chromosome' ORDER BY sr.seq_region_id");

show("regions with a leftover temporary name",
  "SELECT sr.seq_region_id, sr.name, sr.length, cs.name AS coord_system
     FROM seq_region sr JOIN coord_system cs ON cs.coord_system_id = sr.coord_system_id
    WHERE sr.name LIKE '%tmp%' OR sr.name LIKE '%swap%'");

show("meta rows that mention Mt or Pt",
  "SELECT meta_key, meta_value FROM meta
    WHERE meta_value REGEXP '(^|[^[:alnum:]])(Mt|Pt|MT|PT|ChrM|ChrC)([^[:alnum:]]|\$)'");

# Feature counts per region, for every table that carries a seq_region_id.
print "== features on each region, by table ==\n";
my $tables = $dbh->selectcol_arrayref(
  "SELECT table_name FROM information_schema.columns
    WHERE table_schema = ? AND column_name = 'seq_region_id'
      AND table_name NOT IN ('seq_region', 'seq_region_attrib', 'seq_region_synonym', 'seq_region_mapping', 'dna')
    ORDER BY table_name", undef, $dbname);
print join("\t", qw(table Mt Pt)), "\n";
for my $table (@$tables) {
  my %n = map { $_->[0] => $_->[1] } @{ $dbh->selectall_arrayref(
    "SELECT sr.name, COUNT(*) FROM `$table` f JOIN seq_region sr ON sr.seq_region_id = f.seq_region_id
      WHERE sr.name IN ('Mt', 'Pt') GROUP BY sr.name") };
  next unless %n;
  print join("\t", $table, $n{Mt} // 0, $n{Pt} // 0), "\n";
}
print "\n";

show("other databases on this server for the same species",
  "SHOW DATABASES LIKE 'arabidopsis\\_thaliana%'");

# The sequence of each region, rebuilt from the dna table through the assembly table.
for my $name (qw(Mt Pt)) {
  my ($id, $length) = $dbh->selectrow_array(
    "SELECT sr.seq_region_id, sr.length FROM seq_region sr JOIN coord_system cs ON cs.coord_system_id = sr.coord_system_id
      WHERE sr.name = ? AND cs.name = 'chromosome'", undef, $name);
  my $seq = 'N' x $length;
  my ($own) = $dbh->selectrow_array("SELECT sequence FROM dna WHERE seq_region_id = ?", undef, $id);
  if (defined $own) {
    $seq = $own;
  } else {
    my $pieces = $dbh->selectall_arrayref(
      "SELECT a.asm_start, a.asm_end, a.cmp_start, a.cmp_end, a.ori, d.sequence
         FROM assembly a JOIN dna d ON d.seq_region_id = a.cmp_seq_region_id
        WHERE a.asm_seq_region_id = ? ORDER BY a.asm_start", undef, $id);
    for my $piece (@$pieces) {
      my ($asm_start, $asm_end, $cmp_start, $cmp_end, $ori, $sequence) = @$piece;
      my $part = substr($sequence, $cmp_start - 1, $cmp_end - $cmp_start + 1);
      if ($ori == -1) { $part = reverse $part; $part =~ tr/ACGTacgt/TGCAtgca/; }
      substr($seq, $asm_start - 1, $asm_end - $asm_start + 1) = $part;
    }
  }
  open(my $fa, '>', "$out_dir/region_$name.fa") or die "cannot write: $!\n";
  print $fa ">$name\n";
  print $fa substr($seq, $_ * 60, 60), "\n" for 0 .. int((length($seq) - 1) / 60);
  close $fa;
  my $n = ($seq =~ tr/Nn//);
  print "wrote region_$name.fa: " . length($seq) . " bp, $n N\n";
}
$dbh->disconnect;
