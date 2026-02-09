#!/usr/bin/env python3
"""
Convert VCF file to PHYLIP format for phylogenetic analysis
"""

import sys
import os

def vcf_to_phylip(vcf_file, phylip_file):
    """
    Convert VCF file to PHYLIP format
    """
    print(f"Converting {vcf_file} to {phylip_file}")
    
    samples = []
    sequences = {}
    snp_count = 0
    
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                continue  # Skip header lines
            elif line.startswith('#CHROM'):
                # This is the column header line
                parts = line.strip().split('\t')
                samples = parts[9:]  # Get sample names
                print(f"Found {len(samples)} samples: {samples}")
                # Initialize sequences for each sample
                for sample in samples:
                    sequences[sample] = []
                continue
            else:
                # This is a SNP line
                parts = line.strip().split('\t')
                if len(parts) < 10:
                    continue
                
                chrom, pos, snp_id, ref, alt = parts[0:5]
                genotypes = parts[9:]
                
                # Process genotypes for this SNP
                for i, gt in enumerate(genotypes):
                    if i >= len(samples):
                        continue
                    
                    # Extract genotype (first part before :)
                    gt_code = gt.split(':')[0]
                    
                    # Convert to single character code
                    if gt_code == '0/0':
                        sequences[samples[i]].append('A')  # Homozygous reference
                    elif gt_code == '1/1':
                        sequences[samples[i]].append('T')  # Homozygous alternate
                    elif gt_code == '0/1' or gt_code == '1/0':
                        sequences[samples[i]].append('G')  # Heterozygous
                    elif gt_code == './.' or gt_code == '.':
                        sequences[samples[i]].append('N')  # Missing data
                    else:
                        sequences[samples[i]].append('N')  # Unknown
                
                snp_count += 1
    
    print(f"Processed {snp_count} SNPs")
    
    # Write PHYLIP format
    with open(phylip_file, 'w') as f:
        # Header: number of samples and sequence length
        f.write(f"  {len(samples)} {snp_count}\n")
        
        # Write sequences
        for sample, seq in sequences.items():
            # PHYLIP format: 10-character sample name followed by sequence
            sample_name = sample[:10].ljust(10)
            sequence_str = ''.join(seq)
            f.write(f"{sample_name} {sequence_str}\n")
    
    print(f"PHYLIP file created: {phylip_file}")
    
    # Also create a simple FASTA file for alternative methods
    fasta_file = phylip_file.replace('.phy', '.fasta')
    with open(fasta_file, 'w') as f:
        for sample, seq in sequences.items():
            f.write(f">{sample}\n")
            sequence_str = ''.join(seq)
            # Write in blocks of 60 characters
            for i in range(0, len(sequence_str), 60):
                f.write(sequence_str[i:i+60] + '\n')
    
    print(f"FASTA file created: {fasta_file}")
    return samples, snp_count

def main():
    if len(sys.argv) != 3:
        print("Usage: python vcf_to_phylip.py <input.vcf> <output.phy>")
        sys.exit(1)
    
    vcf_file = sys.argv[1]
    phylip_file = sys.argv[2]
    
    if not os.path.exists(vcf_file):
        print(f"Error: VCF file {vcf_file} not found")
        sys.exit(1)
    
    samples, snp_count = vcf_to_phylip(vcf_file, phylip_file)
    print(f"Conversion completed successfully!")
    print(f"Samples: {len(samples)}, SNPs: {snp_count}")

if __name__ == "__main__":
    main()
