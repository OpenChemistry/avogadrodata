#!/bin/sh
# Regenerate the compressed-I/O test fixtures.
#
# These are deliberately produced by the reference command line tools rather
# than by Avogadro's own writer, so the tests prove we can read what the rest
# of the world produces. Run from this directory:
#
#   sh regenerate.sh
#
# Requires: gzip, bzip2, xz, zstd. Everything here is tiny; do not add large
# fixtures, the repository is cloned by every contributor.
set -e

data=../..                 # avogadrodata root, i.e. .../avogadrodata
xyz=$data/data/xyz
sdf=$data/data/sdf
cjson=$data/data/cjson
pdb=$data/data/pdb

rm -f ./*.gz ./*.bz2 ./*.xz ./*.zst ./misnamed.xyz ./magic-only.gz

# --- one file per codec, same source, so decoded bytes must be identical -----
gzip  -c -n "$xyz/methane.xyz" > methane.xyz.gz
bzip2 -c    "$xyz/methane.xyz" > methane.xyz.bz2
xz    -c    "$xyz/methane.xyz" > methane.xyz.xz
zstd  -q -c "$xyz/methane.xyz" > methane.xyz.zst

# --- gzip header variations --------------------------------------------------
# -N stores the original file name and mtime in the header.
gzip -c -N "$xyz/methane.xyz" > methane-named.xyz.gz

# --- concatenated members: cat of two independent streams --------------------
# Decoding must yield the concatenation of both payloads, not just the first.
gzip  -c -n "$xyz/methane.xyz" >  multi-members.xyz.gz
gzip  -c -n "$xyz/H2O.xyz"     >> multi-members.xyz.gz
bzip2 -c    "$xyz/methane.xyz" >  multi-members.xyz.bz2
bzip2 -c    "$xyz/H2O.xyz"     >> multi-members.xyz.bz2
xz    -c    "$xyz/methane.xyz" >  multi-members.xyz.xz
xz    -c    "$xyz/H2O.xyz"     >> multi-members.xyz.xz
zstd  -q -c "$xyz/methane.xyz" >  multi-members.xyz.zst
zstd  -q -c "$xyz/H2O.xyz"     >> multi-members.xyz.zst

# --- damaged inputs: these must fail cleanly, never hang or return junk ------
# Truncated: the trailer (and some of the body) is missing.
gzip -c -n "$xyz/nanotube.xyz" > truncated.xyz.gz.tmp
xz   -c    "$xyz/nanotube.xyz" > truncated.xyz.xz.tmp
zstd -q -c "$xyz/nanotube.xyz" > truncated.xyz.zst.tmp
for c in gz xz zst; do
  size=$(wc -c < truncated.xyz.$c.tmp)
  dd if=truncated.xyz.$c.tmp of=truncated.xyz.$c bs=1 count=$((size / 2)) 2>/dev/null
  rm -f truncated.xyz.$c.tmp
done

# Corrupt body: a valid gzip stream with one bit flipped well inside the
# deflate data. zlib must reject this on the CRC; libarchive would not, which
# is why gzip is decoded with zlib.
gzip -c -n "$xyz/nanotube.xyz" > corrupt.xyz.gz
python3 - <<'PY'
import os
p = 'corrupt.xyz.gz'
b = bytearray(open(p, 'rb').read())
b[len(b) // 2] ^= 0x01
open(p, 'wb').write(bytes(b))
PY

# Magic bytes followed by garbage: must not be mistaken for a real stream.
printf '\037\213\010garbage-not-a-deflate-stream' > magic-only.gz

# --- gzip content under a plain name: content sniffing must win over suffix --
gzip -c -n "$xyz/methane.xyz" > misnamed.xyz

# --- empty payload -----------------------------------------------------------
: > empty.xyz
gzip -c -n empty.xyz > empty.xyz.gz
rm -f empty.xyz

# --- decompression bomb: ~4 MiB of zeros in a few KB -------------------------
# Exercises the size cap and multi-chunk decoding without a large fixture.
dd if=/dev/zero bs=1024 count=4096 2>/dev/null | gzip -9 -c -n > zeros-4mib.gz

# --- other inner formats, including names that mimic real downloads ----------
zstd -q -c "$sdf/multi.sdf"        > multi.sdf.zst
gzip -c -n "$pdb/1CRN.pdb"         > 1crn.pdb.gz
xz    -c   "$cjson/ethane.cjson"   > ethane.cjson.xz

ls -l
