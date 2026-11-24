# Shared settings for the TAIR12 TRaCE run. Every script sources this file.

# The project directory is the directory this file is in.
export PROJ=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export LC_ALL=C

# Every tool puts its temporary files under the project. On the machine this
# was run on, the root filesystem (and so /tmp) had about 1 GB free.
export TMPDIR=$PROJ/tmp
export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$TMPDIR -XX:-UsePerfData"

export PATH=$PROJ/env/bin:$PATH

# TRaCE version used throughout.
export TRACE_COMMIT=836479f6bf55ca823ad7246134bdcb1aa070b267
export TRACE_PL=$PROJ/sw/trace/TRaCE.pl

# Election thresholds (TRaCE README defaults).
export MIN_TPM=0.5 MIN_EVI_COV=0.5 MIN_REF_COV=0.5 MAX_AED=0.5

export IPRSCAN_VERSION=5.78-109.0
export IPRSCAN_HOME=$PROJ/sw/interproscan-$IPRSCAN_VERSION

# The three settings below are specific to the machine. Each keeps the value
# used for this run unless it is already set in the environment.

# TRaCE.pl runs under the system perl. The conda env carries its own perl,
# and a different build can order hash keys differently.
export SYS_PERL=${SYS_PERL:-/usr/local/bin/perl}

# Directory holding the TAIR12 genome and annotation as fetched from TAIR on
# 2026-11-18, with their SHA256SUMS. It is only read.
export SRC=${SRC:-/scratch/olson/gramene-kg/sources/tair12/2026-11-18}

# mamba, used once to create the conda env in env/.
export MAMBA=${MAMBA:-/scratch/olson/sw/miniforge3/bin/mamba}
