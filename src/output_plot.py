import numpy as np
import matplotlib.pyplot as plt

file = open("Pressure.txt")

x, y = file.readline().split(" ")
x = int(x)
y = int(y)

X = np.zeros((x, y))
Y = np.zeros((x, y))
Z = np.zeros((x, y))

for i in range(x):
    for j in range(y):
        px, py, pz = file.readline().split(" ")
        X[i, j] = float(px)
        Y[i, j] = float(py)
        Z[i, j] = float(pz)

fig1, ax2 = plt.subplots(constrained_layout=True)
CS = ax2.contourf(X, Y, Z, 20, origin='lower')

ax2.set_xlabel('x / L')
ax2.set_ylabel('y / L')

CS.cmap.set_under("magenta")
CS.cmap.set_over("yellow")
cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel('P/P_inf')

plt.show()