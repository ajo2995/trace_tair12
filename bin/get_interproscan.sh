#!/usr/bin/env bash
# Download and unpack InterProScan.
#
# The EBI server gives this host about 1.4 MB/s per connection, so the 7 GB
# tarball is fetched as parallel byte ranges and stitched back together. The
# md5 published next to the tarball is checked before unpacking.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
cd "$PROJ/sw"

tarball=interproscan-$IPRSCAN_VERSION-64-bit.tar.gz
url=https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/$IPRSCAN_VERSION/$tarball
streams=8

if [ ! -e "$tarball.ok" ]; then
  total=$(curl -fsSI "$url" | awk 'tolower($1)=="content-length:"{print $2}' | tr -d '\r')
  chunk=$(( (total + streams - 1) / streams ))
  echo "fetching $tarball: $total bytes in $streams ranges"
  pids=()
  for i in $(seq 0 $((streams - 1))); do
    first=$(( i * chunk ))
    last=$(( first + chunk - 1 ))
    [ "$last" -ge "$total" ] && last=$(( total - 1 ))
    want=$(( last - first + 1 ))
    part=$tarball.part$i
    # A part of the right size was finished by an earlier attempt.
    if [ -e "$part" ] && [ "$(stat -c %s "$part")" -eq "$want" ]; then continue; fi
    curl -fsS --retry 10 --retry-delay 15 -r "$first-$last" -o "$part" "$url" &
    pids+=($!)
  done
  for p in ${pids[@]+"${pids[@]}"}; do wait "$p"; done
  for i in $(seq 0 $((streams - 1))); do cat "$tarball.part$i"; done > "$tarball"
  curl -fsS "$url.md5" | md5sum -c -
  rm -f "$tarball".part*
  touch "$tarball.ok"
fi

if [ ! -e "$IPRSCAN_HOME/.unpacked" ]; then
  echo "unpacking"
  tar -pxzf "$tarball"
  touch "$IPRSCAN_HOME/.unpacked"
fi
echo "InterProScan ready in $IPRSCAN_HOME"
du -sh "$IPRSCAN_HOME" | cut -f1
ls "$IPRSCAN_HOME/data/pfam"
