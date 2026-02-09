# scripts/calculate_pic.py
import pandas as pd
import numpy as np

def calculate_pic(p, q):
    """Calculate Polymorphic Information Content for biallelic SNP"""
    return 1 - (p**2 + q**2) - 2*(p**2)*(q**2)

# Read allele frequencies (we need to generate this first)
print("First, let's generate allele frequencies with PLINK...")
