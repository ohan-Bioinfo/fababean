import numpy as np

print("=== DEBUGGING KINSHIP FORMAT ===")

# Read sample IDs
with open('Diversity/Kinship/Faba_Kinship.king.id', 'r') as f:
    samples = [line.strip().split()[1] for line in f.readlines()]

print(f"Number of samples in .king.id: {len(samples)}")
print(f"Sample names: {samples}")

# Read kinship file line by line
print("\n=== KINSHIP FILE ANALYSIS ===")
kinship_lines = []
with open('Diversity/Kinship/Faba_Kinship.king', 'r') as f:
    for i, line in enumerate(f):
        values = [float(x) for x in line.strip().split()]
        kinship_lines.append(values)
        print(f"Line {i+1}: {len(values)} values - first 3: {values[:3]}")

print(f"\nTotal kinship lines: {len(kinship_lines)}")

# Check if this matches the expected triangular number
n_expected = len(samples)
triangular_number = n_expected * (n_expected - 1) // 2
print(f"Expected lines for {n_expected} samples (triangular): {triangular_number}")

# Let's try to understand the structure
print("\n=== MATRIX RECONSTRUCTION ATTEMPT ===")
# The kinship file might be just the lower triangular part without diagonal
# Let's try to build the matrix assuming it's for n-1 samples?

# Alternative: Maybe the first sample is missing from kinship?
if len(kinship_lines) == len(samples) - 1:
    print("Possible: kinship file has n-1 lines")
    # Try building matrix for n-1 samples
    n_kinship = len(kinship_lines)
    test_matrix = np.zeros((n_kinship, n_kinship))
    for i in range(n_kinship):
        for j in range(len(kinship_lines[i])):
            test_matrix[i, j] = kinship_lines[i][j]
            test_matrix[j, i] = kinship_lines[i][j]
    print(f"Test matrix shape: {test_matrix.shape}")
    
# Alternative: Maybe it's a different format?
print(f"\nLine length pattern:")
for i in range(min(5, len(kinship_lines))):
    print(f"  Line {i+1}: {len(kinship_lines[i])} elements")

# Check if it's a proper lower triangular matrix
is_triangular = all(len(kinship_lines[i]) == i+1 for i in range(len(kinship_lines)))
print(f"Is proper lower triangular: {is_triangular}")

# Count total elements
total_elements = sum(len(line) for line in kinship_lines)
print(f"Total kinship values: {total_elements}")
