import numpy as np 
import matplotlib.pyplot as plt

data = np.genfromtxt('test.csv',delimiter = ',')


fig = plt.figure()
ax=fig.add_subplot(111)

ax.imshow(data,cmap='gray')

fig.savefig("ring.png")