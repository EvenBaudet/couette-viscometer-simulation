import matplotlib.pyplot as plt
import numpy as np

fig,ax = plt.subplots()
ax.clear
fig.set_size_inches(10,7)
ax1 = fig.add_axes([2/3-1/20, 5/8-1/16, 1/4+1/10, 1/5+1/20])
ax2 = fig.add_axes([2/3-1/20, 3/8-1/16, 1/4+1/10, 1/5-1/20])
ax3 = fig.add_axes([2/3-1/20, 1/8-1/16, 1/4+1/10, 1/5-1/20])

ax1.set_title("évolution de v(r)")
ax1.set_ylabel("v(r) (m/s)")
ax1.set_xlabel("r (cm)")
ax1.grid(True)

ax2.set_title("évolution de ω1(t)")
ax2.set_ylabel("ω1(t) (rad)")
ax2.set_xlabel("t (s)")
ax2.grid(True)

ax3.set_title("évolution de ω2(t)")
ax3.set_ylabel("ω2(t) (rad)")
ax3.set_xlabel("t (s)")
ax3.grid(True)


Omega = [0,0]
V = np.exp([-i/10 for i in range (100)])
t = 100
dt = 1
T = 100

def graph_omega(ax, Omega, t, dt, T, Couleur):
    ax.plot([t-dt,t], Omega, color = Couleur)
    ax.set_xlim(t-4*T,t+2*T)

def graph_v(ax,V,N,R1,R2,couleur):
    R = np.linspace(R1,R2,N)
    ax.plot(R,V,color = couleur)


graph_v(ax1,V,100,20,40,'red')

for i in range (850):
    Omega[0],Omega[1] = Omega[1],np.sin(2*np.pi/100*i)
    t = t + dt
    graph_omega(ax2, Omega, t, dt, T/2, 'blue')
    graph_omega(ax3, Omega, t, dt, T/2, 'green')

plt.show()
