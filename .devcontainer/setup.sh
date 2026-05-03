#!/usr/bin/env bash
set -e

echo "Installing system dependencies..."
sudo apt-get update
sudo apt-get install -y \
  default-jre \
  curl \
  git \
  graphviz

# Dashboard dependencies
pip install -r scripts/requirements.txt

echo "Installing Nextflow..."
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

echo "Setup complete."
