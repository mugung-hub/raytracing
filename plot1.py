import numpy as np 
import matplotlib.pyplot as plt

data = np.genfromtxt('with_gravity.csv',delimiter = ',')
data1 = np.genfromtxt('without_gravity.csv',delimiter = ',')


fig = plt.figure()
ax=fig.add_subplot(111)

ax.imshow(data,cmap='gray')

fig.savefig("with_gravity.png")

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.imshow(data1, cmap='gray')
fig1.savefig("without_gravity.png")
