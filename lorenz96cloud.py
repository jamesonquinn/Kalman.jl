from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

# these are our constants
N = 36  # number of variables
F = 8  # forcing

def Lorenz96(x,t):

  # compute state derivatives
  d = np.zeros(N)
  # first the 3 edge cases: i=1,2,N
  d[0] = (x[1] - x[N-2]) * x[N-1] - x[0]
  d[1] = (x[2] - x[N-1]) * x[0]- x[1]
  d[N-1] = (x[0] - x[N-3]) * x[N-2] - x[N-1]
  # then the general case
  for i in range(2, N-1):
      d[i] = (x[i+1] - x[i-2]) * x[i-1] - x[i]
  # add the forcing term
  d = d + F

  # return the state derivatives
  return d

x0 = F*np.ones(N) # initial state (equilibrium)
x0[19] += 0.01 # add small perturbation to 20th



cloudsd = .1
cloudn = 100
cloudSteps = 5
cloudStepSize = .125
results = np.zeros([cloudn*cloudSteps,N])
results.shape

for i in range(0,cloudn):
    x00 = np.copy(x0)
    for d in range(1,5):
        x00[d] += np.random.normal(0.,cloudsd)
    t = np.arange(0.0, cloudSteps*cloudStepSize, cloudStepSize)
    xs = odeint(Lorenz96, x00, t)
    for j in range(0,cloudSteps):
        results[j*cloudn+i,:] = xs[j,:]

t = np.arange(0.0, 30.0, 0.1)

x = odeint(Lorenz96, x0, t)

# plot first three variables
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.scatter(results[:,0],results[:,1],results[:,2],
        c=np.repeat(range(0,cloudSteps),cloudn))
ax.set_xlabel('$x_1$')
ax.set_ylabel('$x_2$')
ax.set_zlabel('$x_3$')
plt.show()
x.shape



# rotate the axes and update
for angle in range(0, 360):
    ax.view_init(30, angle)
    plt.draw()
    plt.pause(.001)
