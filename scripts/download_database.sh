#!/usr/bin/env bash
set -euo pipefail

# Variables
DB_URL="https://zenodo.org/records/5571251/files/checkm2_database.tar.gz?download=1"
DB_DIR="module_db"
DB_TAR="checkm2_database.tar.gz"

# Create target directory
mkdir -p "${DB_DIR}"
cd "${DB_DIR}"

# Download database
echo "Downloading CheckM2 database..."
curl -L "${DB_URL}" -o "${DB_TAR}"

# Extract database
echo "Extracting database..."
tar -xvzf "${DB_TAR}"

# Optional: remove tarball to save space
rm -f "${DB_TAR}"

echo "CheckM2 database setup complete."
