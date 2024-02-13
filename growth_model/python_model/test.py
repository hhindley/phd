# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import numpy as np 

x = np.linspace(0, 20, 100 )
fig = plt.plot(x, np.sin(x))
plt.show() 

# fig.savefig("test.png")

