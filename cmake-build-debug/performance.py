import matplotlib.pyplot as plt
import numpy as np

# Load the results
dataCache = np.loadtxt('performanceCache')
dataNoCache = np.loadtxt('performanceWithoutCache')
# Extract n and time
n = dataCache[:, 0]
timeNoCache = dataNoCache[:, 1]
timeCache = dataCache[:, 1]

# Plot
plt.figure(figsize=(10, 6))
plt.plot(n, timeNoCache, marker='o', label='without Cache')
plt.plot(n, timeCache, marker='o', label='with Cache')
plt.xlabel('Number of grid points (n)')
plt.ylabel('Elapsed time (sec)')
plt.title('Performance Test')
plt.legend()
plt.grid(True)
plt.savefig('performance_plot.png')  # Save the plot as a PNG file
plt.show()
