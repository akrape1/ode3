from scipy.integrate import odeint
import numpy as np
from math import sqrt


print("*** This example used a depricated interface.  Please use projScPY2 instead ***")


# m, air_k specified as parameters
params=[1.0, 0.12]

def func(y,t):
    g = 9.81
    m = params[0]
    air_k = params[1]
    v = sqrt(y[1]*y[1] + y[3]*y[3])
    f0 = y[1]                                # f_ri
    f1 = -air_k * v * y[1] / m         # f_vi
    f2 = y[3]                                # f_rj
    f3 = -air_k * v * y[3] / m - g    # f_vj
    return [f0,f1,f2,f3]

#starting coordinates
y0=[0,10,0,10]   # x=0, vx=10 m/s, y=0, vy=10 m/s

t = np.linspace(0,1.6,200)
sol = odeint(func,y0,t)
#for y in us:
#    print(y)

from matplotlib import pyplot as plt

#y=np.array(y)
x=np.array(sol[ : , 0])
vx=np.array(sol[ : , 1])
y=np.array(sol[ : , 2])
vy=np.array(sol[ : , 3])
ke = 0.5 * params[0] * (vx*vx+vy*vy)

plt.subplot(1, 2, 1)
plt.plot(x,y,'-')
plt.xlabel('x values [m]')
plt.ylabel('y values [m]')
plt.title('Height vs distance')

plt.subplot(1, 2, 2)
plt.plot(t,ke,'-')
plt.xlabel('time [s]')
plt.ylabel('KE [J]')
plt.title('KE vs t')

plt.tight_layout() 
print("Close plot window to exit")
plt.show()


