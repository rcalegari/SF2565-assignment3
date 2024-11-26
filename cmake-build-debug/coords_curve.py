import matplotlib.pyplot as plt

# File containing bottom curve coordinates
filename = "bottom_curve_coordinates.txt"

# Read the data
x_coords, y_coords = [], []
with open(filename, 'r') as file:
    for line in file:
        x, y = map(float, line.strip().split())
        x_coords.append(x)
        y_coords.append(y)

# Plot the curve
plt.figure(figsize=(8, 6))
plt.plot(x_coords, y_coords, label="Bottom Curve", color="blue")
plt.title("Bottom Curve Visualization")
plt.xlabel("X")
plt.ylabel("Y")
plt.grid(True)
plt.legend()
plt.show()
