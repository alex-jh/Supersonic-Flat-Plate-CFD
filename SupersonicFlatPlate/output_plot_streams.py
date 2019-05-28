import numpy as np
import matplotlib.pyplot as plt

file = open("../SupersonicFlowGpu/Output2/VelocityX.txt")
#file = open("Pressure.txt")

x, y = file.readline().split(" ")
x = int(x)
y = int(y)

X = np.zeros((x, y))
Y = np.zeros((x, y))
Z = np.zeros((x, y))

for i in range(x):
    for j in range(y):
        px, py, pz = file.readline().split(" ")
        X[i, j] = float(py)
        Y[i, j] = float(px)
        Z[i, j] = float(pz)

u = Z

file = open("../SupersonicFlowGpu/Output2/VelocityY.txt")
#file = open("Pressure.txt")

x, y = file.readline().split(" ")
x = int(x)
y = int(y)

Z = np.zeros((x, y))

for i in range(x):
    for j in range(y):
        px, py, pz = file.readline().split(" ")
        Z[i, j] = float(pz)

v = Z

fig, ax = plt.subplots()

print (X, Y)

ax.streamplot(X,Y,u,v, density = 1)
#ax.axis([0.5,2.1,0,2])
#ax.xaxis.set_ticks([])
#ax.yaxis.set_ticks([])
ax.set_title('Stream Plot of Field Lines')


plt.show()