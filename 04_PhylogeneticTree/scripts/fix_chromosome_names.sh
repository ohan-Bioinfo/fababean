#!/bin/bash

echo "Fixing chromosome names in BIM file..."

# Create a backup
cp data/faba_fingerprint.bim data/faba_fingerprint.bim.backup

# Fix chromosome names - convert to simple numeric codes
awk '{
    # Remove "chr" prefix if present
    chrom = $1;
    gsub(/^chr/, "", chrom);
    
    # Map chromosome names to numeric codes
    if (chrom == "1L") $1 = "1";
    else if (chrom == "1S") $1 = "2";
    else if (chrom == "2") $1 = "3"; 
    else if (chrom == "3") $1 = "4";
    else if (chrom == "4") $1 = "5";
    else if (chrom == "5") $1 = "6";
    else if (chrom == "6") $1 = "7";
    else if (chrom ~ /^sca/) $1 = "99";  # Put all scaffolds in chr99
    else $1 = chrom;
    
    print $0
}' data/faba_fingerprint.bim.backup > data/faba_fingerprint.bim

echo "Chromosome names fixed. Checking results:"
cut -f1 data/faba_fingerprint.bim | sort | uniq -c

echo "Creating new FAM file with proper sample names..."
# Ensure FAM file has proper sample names
cp data/faba_fingerprint.fam data/faba_fingerprint.fam.backup
awk '{print $1, $2, $3, $4, $5, $6}' data/faba_fingerprint.fam.backup > data/faba_fingerprint.fam

echo "Files updated successfully!"
