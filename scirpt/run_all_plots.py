#!/usr/bin/env python3
"""
Main script to run all SNP visualization plots
"""

import subprocess
import sys
import os

def run_script(script_name):
    """Run a Python script and return success status"""
    try:
        print(f"\n{'='*50}")
        print(f"Running {script_name}...")
        print(f"{'='*50}")
        result = subprocess.run([sys.executable, script_name], 
                              capture_output=True, text=True)
        if result.returncode == 0:
            print(f"✓ {script_name} completed successfully")
            return True
        else:
            print(f"✗ {script_name} failed with error:")
            print(result.stderr)
            return False
    except Exception as e:
        print(f"✗ Error running {script_name}: {e}")
        return False

def main():
    """Run all visualization scripts"""
    scripts = [
        'plot_snp_distribution.py',
        'plot_filtering_pipeline.py', 
        'plot_missingness.py',
        'plot_heterozygosity.py',
        'plot_dashboard.py'
    ]
    
    # Check which scripts exist
    available_scripts = [s for s in scripts if os.path.exists(s)]
    missing_scripts = [s for s in scripts if not os.path.exists(s)]
    
    if missing_scripts:
        print("The following scripts are missing:")
        for script in missing_scripts:
            print(f"  - {script}")
        print("\nPlease create the missing scripts first.")
    
    # Run available scripts
    success_count = 0
    for script in available_scripts:
        if run_script(script):
            success_count += 1
    
    print(f"\n{'='*50}")
    print(f"Summary: {success_count}/{len(available_scripts)} scripts completed successfully")
    print(f"{'='*50}")

if __name__ == "__main__":
    main()
