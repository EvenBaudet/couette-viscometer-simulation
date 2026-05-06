import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import identity, csc_matrix
from scipy.sparse.linalg import spsolve
from matplotlib.animation import FuncAnimation
from matplotlib.colors import Normalize

from abstractresolver import Resolver

class BackwardsEuler(Resolver):
    # params d'init
    # nu = 1    # viscosité cinématique
    # R1 = 0.3    # Rayon Intérieur (r1)
    # R2 = 1.0    # Rayon Extérieur (r2)
    # nr = 100    # Nombre de points radiaux
    # dr = (R2 - R1) / (nr - 1)    # Pas spatial radial
    # nt = 400    # Nombre de pas de temps
    # dt = 0.0005 # Pas de temps

    fps = 30         # Frames par seconde pour l'animation
    skip_frames = 5  # Pas de temps à sauter pour chaque frame d'animation

    # Grille radiale 1D (commence à R1)
    def __init__(self, visc):

        self.visc = visc
        #self.r_grid_1d = np.linspace(visc.R1, visc.R2, nr)

        # Nombre de diffusion 'r_diff'
        #self.r_diff = visc.nu * dt / visc.dr**2

        #Condition Initiale 
        #self.u_initial = 0*self.r_grid_1d # vitesse nulle partout à t = 0 
        #self.u = self.u_initial.copy()

        self.reload()
        # Construction de la Matrice du Système M pour Euler Rétrograde

        # M est initialisé avec l'identité (le terme u^n)

        # Remplir la matrice M avec les coefficients du schéma implicite
        # i va de 1 à nr-2 (points internes, non CL)

    def build_matrice(self):
        M = identity(self.visc.nr, format='csc') 

        r_diff = self.visc.r_diff
        dr = self.visc.dr
        nr  = self.visc.nr
        for i in range(1, nr - 1):
            r_i = self.visc.r_grid_1d[i]
            
            # 1. Sous-diagonale (i, i-1)
            M[i, i-1] = -r_diff * (1 - dr / (2 * r_i))
            
            # 2. Diagonale principale (i, i)
            M[i, i] = 1 + 2 * r_diff
            
            # 3. Sur-diagonale (i, i+1)
            M[i, i+1] = -r_diff * (1 + dr / (2 * r_i))

        # Conditions aux Limites
        # 1. Bord Intérieur (r = R1, i = 0)
        M[0, :] = 0
        M[0, 0] = 1

        # 2. Bord Extérieur (r = R2, i = nr-1)
        # M[nr - 1, :] = 0
        # M[nr - 1, nr - 1] = 1

        self.M = M

    def reload(self):
        self.build_matrice()

    def calc_next(self):
        b = self.visc.v.copy()
        
        # Imposer les conditions aux limites sur le vecteur RHS 'b' (T=0 aux bords)
        b[0] = self.visc.currentv1()        # Bord Intérieur r=R1
        b[self.visc.nr - 1] = b[self.visc.nr - 2]   # Bord Extérieur r=R2

        return spsolve(self.M, b)
        
        # if t % skip_frames == 0 and t > 0:
        #     vitesses.append(u.copy())