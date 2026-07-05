#!/bin/bash
set -euo pipefail

# 1. Compile afin with the conda toolchain and install it onto PATH.
#    The makefile honors $CXX/$CXXFLAGS/$LDFLAGS and the PREFIX/DESTDIR vars.
pushd afin
make CXX="${CXX}"
make install PREFIX="${PREFIX}"
popd

# 2. Install Fast-Plast's own data + helper scripts under share, driver in bin.
SHARE="${PREFIX}/share/fast-plast"
mkdir -p "${SHARE}/bin" "${SHARE}/Coverage_Analysis"
cp -R bin/.            "${SHARE}/bin/"
cp -R Coverage_Analysis/. "${SHARE}/Coverage_Analysis/"

# Decompress the bundled reference plastomes if shipped gzipped.
if [ -f "${SHARE}/bin/GenBank_Plastomes.gz" ]; then
    gunzip -f "${SHARE}/bin/GenBank_Plastomes.gz"
fi

# Drop any vendored third-party binaries that may have been committed to bin/
# (conda provides these as dependencies; we never ship our own copies).
rm -rf "${SHARE}/bin/"{Trimmomatic-*,fastp-*,pigz-*,bowtie2-*,bowtie-*,SPAdes-*,ncbi-blast-*,sspace_basic-*,jellyfish-*} 2>/dev/null || true

# 3. Install the driver and a convenience symlink.
install -m 0755 fast-plast.pl "${PREFIX}/bin/fast-plast.pl"
ln -sf fast-plast.pl "${PREFIX}/bin/fast-plast"

# 4. Belt-and-suspenders: the driver already finds its data via the ../share
#    layout, but also export FASTPLAST_SHARE on activation for robustness.
mkdir -p "${PREFIX}/etc/conda/activate.d"
cat > "${PREFIX}/etc/conda/activate.d/fast-plast.sh" <<'EOF'
export FASTPLAST_SHARE="${CONDA_PREFIX}/share/fast-plast"
EOF
