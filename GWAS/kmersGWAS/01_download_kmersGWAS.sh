#!/bin/bash
set -e

# Download kmersGWAS v0.2-beta
wget https://github.com/voichek/kmersGWAS/releases/download/v0.2-beta/v0_2_beta.zip
unzip v0_2_beta.zip

# Move to a standard directory name for easier scripting
if [ -d v0_2_beta ]; then
    mv v0_2_beta kmersGWAS
elif [ -d v0_2_beta/ ]; then
    mv v0_2_beta/ kmersGWAS
fi

# Make sure KMC and GEMMA are executable
chmod +x external_programs/kmc_v3
chmod +x external_programs/gemma_0_96

# Optionally add to PATH (for current session)
#export PATH="$(pwd)/kmersGWAS/external_programs/KMC/bin:$(pwd)/kmersGWAS/external_programs/GEMMA/bin:$PATH"

echo "kmersGWAS, KMC, and GEMMA are set up."
