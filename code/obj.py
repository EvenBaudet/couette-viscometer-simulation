import numpy as np
from eulerback import BackwardsEuler

class Viscosimetre(object):
    def __init__(self, f1=np.cos, A=0, Resolv=BackwardsEuler,N=24,nu=None,eta=5e-2,rho=1e3, R1=2e-2,R2=5e-2, dt=5e-4, r_diff=1):
        self.Resolv = Resolv # garde la classe pour rebuild

        self.nr = N
        self.A = A
        self.t = 0
        self.nu = nu if nu else eta/rho
        self.R1,self.R2 = R1,R2
        self.dt = dt
        self.r_diff = r_diff

        self.f1 = f1

        # self.r_diff = self.nu * dt / self.dr**2
        # self.resolv = Resolv(self)

        self.liste_t = []
        self.liste_v = []
        self.l_omega_1 = []
        self.l_omega_2 = []

        self.rebuild()
        self.append_listes()

    def reload(self):
        self.resolv = self.Resolv(self)

    def rebuild(self):
        self.dr = (self.R2-self.R1)/self.nr
        # self.nu = self.eta/self.rho

        self.r_grid_1d = np.linspace(self.R1, self.R2, self.nr)
        v0 = 0*self.r_grid_1d # vitesse initiale
        self.v = v0 # pour compatibilite
        self.reload()

    def append_listes(self):
        self.liste_t.append(self.t)
        # self.liste_v.append(self.v.copy())
        self.l_omega_1.append(self.current_omega1())
        self.l_omega_2.append(self.current_omega2())

    def current_omega1(self):
        # evite de tout calculer
        return self.v[0]/self.R1

    def current_omega2(self):
        # evite de tout calculer
        return self.v[-1]/self.R2

    def omega(self):
        return self.v/self.r_grid_1d

    def currentv1(self):
        return self.A*self.R1*self.f1(self.t)

    def update(self):
        self.v = self.resolv.calc_next()
        self.t += self.dt

if __name__ == "__main__":
    visc = Viscosimetre()
    print(visc.v)
    for k in range(3):
        visc.update()
        print(visc.v)