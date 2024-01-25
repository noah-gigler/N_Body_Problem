import matplotlib.pyplot as plt
import numpy as np

# Load the data from the file
data = np.loadtxt('output/direct_vs_approx.txt')

# Split the data into three vectors
x = data[:, 0]
s0 = data[:, 1]
s0_5 = data[:, 2]
s1 = data[:, 3]
s1_5 = data[:, 4]
n = data[:, 5]

plt.plot(x, s0, label='s = 0.0 * mean distance')
plt.plot(x, s0_5, label='s = 0.05 * mean distance')
plt.plot(x, s1, label='s = 0.1 * mean distance')
plt.plot(x, s1_5, label='s = 1.0 * mean distance')
plt.plot(x, n, label='newtonian approximation')
plt.xlabel('radius')
plt.ylabel('force')
ax = plt.gca()
plt.yscale('log')
plt.legend()
plt.title('Direct vs. Approximate Force Calculation')

plt.savefig('output/newton_bins500_r3.png', bbox_inches='tight', dpi = 300)
plt.show()