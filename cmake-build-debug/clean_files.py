import numpy as np

xPath = 'xGrid.txt'
yPath = 'yGrid.txt'

with open(xPath, 'r') as file:
    print("xGrid.txt:")
    for i, line in enumerate(file):
        print(f"Row {i}: {repr(line)}")

with open(yPath, 'r') as file:
    print("\nyGrid.txt:")
    for i, line in enumerate(file):
        print(f"Row {i}: {repr(line)}")

def clean_file(filepath):
    with open(filepath, 'r') as file:
        cleaned_lines = [line.strip() for line in file if line.strip()]  # Remove empty lines and strip spaces
    with open(filepath, 'w') as file:
        file.write('\n'.join(cleaned_lines))

clean_file('xGrid.txt')
clean_file('yGrid.txt')

def load_clean_grid(filepath):
    with open(filepath, 'r') as file:
        return np.array([
            list(map(float, line.split()))  # Convert each space-separated value to float
            for line in file if line.strip()  # Skip empty lines
        ])

X = load_clean_grid(xPath)
Y = load_clean_grid(yPath)

print("Shape of X:", X.shape)
print("Shape of Y:", Y.shape)
print("X:\n", X)
print("Y:\n", Y)

