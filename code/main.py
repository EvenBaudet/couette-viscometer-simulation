import matplotlib.pyplot as plt
from matplotlib.quiver import Quiver
import numpy as np
import matplotlib.animation as animation
from matplotlib.widgets import Button, Slider, TextBox, CheckButtons
import itertools

from obj import Viscosimetre

from graphing import graph_omega, graph_v

cos, sin, tanh, pi = np.cos, np.sin, np.tanh, np.pi

fig, ax = plt.subplots()
x_lim, y_lim = 6e-2, 6e-2
fig.set_size_inches(11,7)
fig.subplots_adjust(top = 0.92, bottom = 0.22, left = 0.05, right = 0.5)
ax.set_xlim(-x_lim, x_lim)
ax.set_ylim(-y_lim, y_lim)

fig.text(0.02, 0.95, "Simulation d'un viscosimètre de Couette |", size = 13, fontweight='bold', color = "#A22020")
fig.text(0.035, 0.2, "vitesse du cylindre intérieur", size = 11, fontweight='bold', color = "#A22020")
fig.text(0.295, 0.2, "paramètres de la simulation", size = 11, fontweight='bold', color = "#A22020")
fig.text(0.68, 0.95, "paramètres temporels fixés", size = 10, fontweight='bold', color = "#A22020")
fig.text(0.295, 0.055, "discrétisation radiale", size = 10, fontweight='bold', color = "#A22020")

#recentrage des axes
ax.spines['left'].set_position('center')
ax.spines['bottom'].set_position('center')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')

ax.set_xticklabels([])
ax.set_yticklabels([])
ax.tick_params(axis='x', which='both', length=6, direction='inout')  
ax.tick_params(axis='y', which='both', length=6, direction='inout')

# paramètres du problème
t = 0
nu = 1.5 #e-5 m2/s
R1, R2 = 2e-2, 5e-2
rho = 1000 # masse volumique (en kg/m3)
eta = 0.05 # viscosite (en Pa*s)
l_p = 1e-3 #longueur d'un pixel dans l'interface
r_diff = 1
nr = 30
dr = (R2 - R1)/nr
dt = (r_diff*dr**2)/nu
updating = False #variable décrite plus loin (fonction signal_click et visual_click)
lancement = False #variable décrivant l'état (demarre/pause) de l'animation

# fonction imposée au premier cylindre
def signal_cosinus(t):
    return np.cos(2*np.pi*omega_frequence*t)

def signal_carre(t):
    t2 = (t*omega_frequence)%1
    return np.sign(0.5-t2) # 1 0 si t2=0.5, 1 si t2 < 0.5, -1 sinon

omega_amplitude = 10.28
omega_frequence = 0

# positions des cylindres
theta_1, theta_2 = 0,0

visc = Viscosimetre(f1=signal_cosinus,A=0,N=nr, eta=eta,rho=rho, R1=R1, R2=R2, dt=dt, r_diff=r_diff)

# représentation des deux cylindres
cercle_int = plt.Circle((0,0), R1, edgecolor="#104077", fill=False, linewidth = 2)
cercle_ext = plt.Circle((0,0), R2, edgecolor="#104077", fill=False, linewidth = 2)
ax.add_patch(cercle_int)
ax.add_patch(cercle_ext)

# représentation de segments collés aux cylindres qui traduiront la rotation du cylindre
repere_1, = plt.plot([-R1, -R1 - (R2-R1)/8], [0, 0], color = "#104077", linewidth = 2)
repere_2, = plt.plot([R1, R1 + (R2-R1)/8], [0, 0], color = "#104077", linewidth = 2)
repere_3, = plt.plot([-R2, -R2 - (R2-R1)/8], [0, 0], color = "#104077", linewidth = 2)
repere_4, = plt.plot([R2, R2 + (R2-R1)/8], [0, 0], color = "#104077", linewidth = 2)


#ajout des sliders pour choisir les rayons, l'amplitude, la frequence et le nombre de rayons modélisé entre R1 et R2
def cree_slider(axes, title, min, max, init, valfmt):
    ax_slider = plt.axes(axes)
    return Slider(ax_slider, title, min, max, valinit=init, valfmt=valfmt)

slider_R1 = cree_slider([0.32, 0.16, 0.15, 0.02], 'R1 ', 1, 2.5, 2, valfmt='%1.2f cm')
slider_R2 = cree_slider([0.32, 0.13, 0.15, 0.02], 'R2 ', 3, 5.5, 5, valfmt='%1.2f cm')
#l'amplitude peut être choisie négative car si la frequence est nulle cela nous permet d'avoir un cylindre tournant en sens inverse
slider_A = cree_slider([0.1, 1/40, 0.08, 0.02], 'amplitude ', -10, 10, 0, valfmt='%1.2f rad.s-1')
slider_f = cree_slider([0.1, 3/40, 0.1, 0.03], 'frequence ', 0, 1, 0, valfmt='%1.2f Hz')

ax_slider = plt.axes([0.32, 0.02, 0.08, 0.02])
slider_nr = Slider(ax_slider, 'nr ', 15, 50, valinit = visc.nr, valfmt='%1.2f points', valstep = 1)

# fonctions pour la visualisation eulerienne
def reload_points(visc):
    global theta, r, points
    R1,R2,nr,dr = visc.R1,visc.R2,visc.nr,visc.dr
    theta, r = np.array([(k*np.pi/2)%2*np.pi for k in range(nr - 4)]), np.linspace(R1 + 2*dr,R2 - 2*dr,nr - 4)
    points = [ax.plot(r[i]*cos(theta[i]), r[i]*sin(theta[i]), 'o', color = "#473636", markersize = 5) for i in range(nr - 4)]

def reload_lignes(visc):
    global lignes_champ
    lignes_champ = [[] for _ in range(visc.nr - 4)]

reload_points(visc)
lignes_champ = [[]for _ in range(nr - 4)]
vecteurs = []

# mise à jour des rayons/discretisation spatiale
def changement_R(val):
    global lignes_champ, points, r, theta
    visc.R1, visc.R2, visc.nr = slider_R1.val*1e-2, slider_R2.val*1e-2, int(slider_nr.val)
    for i in range(len(points)):
        for ligne in lignes_champ[i]:
            ligne[0].remove()
        points[i][0].remove()
    ax1.set_xlim(visc.R1,visc.R2)
    ax1.set_xticks(np.linspace(visc.R1, visc.R2, 5))

    # reconstruction lourde (plus les memes)
    visc.rebuild()
    if check_visuel.get_status()[0]:
        reload_points(visc)
        reload_lignes(visc)
        
slider_R1.on_changed(changement_R)
slider_R2.on_changed(changement_R)
slider_nr.on_changed(changement_R)

#paramètres temporels

text_dt = fig.text(0.7, 0.915, 'pas de temps: '+str(round(visc.dt*1e3,3))+' ms', size = 10)
text_t = fig.text(0.7, 0.885, 'temps: '+str(round(visc.t,3))+' s', size = 10)

#Ajout du bouton play et stop 

ax_icone = plt.axes([0.91, 0.88, 0.05, 0.05])
icone = Button(ax_icone, "▶", color = "#88B4E5")

def demarre(event):
    global lancement
    if lancement:
        ani.event_source.stop()
        icone.label.set_text("▶")
    else:
        ani.event_source.start()
        icone.label.set_text("| |")
    lancement = not lancement

icone.on_clicked(demarre)

#variable du fluide
def changement_nu(val):
    try:
        visc.nu = float(val)
    except:
        pass
    else:
        visc.reload()


ax_boite = plt.axes([0.32, 0.085, 0.1, 0.03])
boite_nu = TextBox(ax_boite, 'ν ', initial = str(nu))
boite_nu.on_text_change(changement_nu)
fig.text(0.43, 0.095, "e-5 m2.s-1")

# Ajout du choix des fonctions de vitesse à imposer au cylindre intérieur
# la variable updating permet d'éviter les appels récursifs lié aux méthodes set_active des objets checkButtons
def signal_click(label):
    global updating

    if updating:
        return
    
    updating = True
    if label == "sinusoïdal" and check_signal.get_status()[0]:
        if check_signal.get_status()[1]:
            check_signal.set_active(1)
            visc.f1 = signal_cosinus

    elif label == "carré" and check_signal.get_status()[0]:
        if check_signal.get_status()[0]:
            check_signal.set_active(0)
            visc.f1 = signal_carre
    updating = False
    
def changement_signal(val):
    global omega_frequence
    visc.A = slider_A.val
    omega_frequence = slider_f.val
    # changement du signal d'entree, pas besoin de reload le resolveur

slider_A.on_changed(changement_signal)
slider_f.on_changed(changement_signal)


rax = plt.axes([0.1, 0.105, 0.1, 0.1])
rax.set_frame_on(False)
check_signal = CheckButtons(ax = rax, labels = ["sinusoïdal", "carré"])
check_signal.set_active(0)
check_signal.on_clicked(signal_click)

# Ajout des boutons qui permettent de changer le visuel (lagrangien/eulerien)
# la variable updating permet d'éviter les appels récursifs lié aux méthodes set_active des objets checkButtons
def visuel_click(label):
    global lignes_champ, points, r, theta, updating

    if updating:
        return

    updating = True
    if label == "Lagrangienne" and check_visuel.get_status()[1] or label == "Eulerienne" and not check_visuel.get_status()[0]:
        for i in range(len(points)):
            for ligne in lignes_champ[i]:
                ligne[0].remove()
            points[i][0].remove()
        points = []
            
        if check_visuel.get_status()[0]:
            check_visuel.set_active(0)

    elif label == "Eulerienne" and check_visuel.get_status()[0]:
        for objet in ax.collections[:]: 
            if isinstance(objet, Quiver):
                objet.remove()
        reload_lignes(visc)
        reload_points(visc)
        if check_visuel.get_status()[1]:
            check_visuel.set_active(1)

    updating = False

fig.text(0.415, 0.95, "description du mouvement: ", size = 10, fontweight='bold', color = "#2052A2")
rax = plt.axes([0.45, 0.86, 0.1, 0.1])
rax.set_frame_on(False)
check_visuel = CheckButtons(ax = rax, labels = ["Eulerienne", "Lagrangienne"], label_props = {'color': ["#2052A2"]})
check_visuel.set_active(0)
check_visuel.on_clicked(visuel_click)

# initialisation graphiques a droite de l'interface

ax1 = fig.add_axes([0.62, 0.6, 0.35, 0.22])
ax2 = fig.add_axes([0.62, 0.34, 0.35, 0.15])
ax3 = fig.add_axes([0.62, 0.08, 0.35, 0.15])

ax1.set_title("évolution de v(r)")
ax1.set_ylabel("v(r) (m.s-1)")
ax1.set_xlabel("r (cm)")
ax1.set_xlim(R1,R2)
ax1.set_xticks(np.linspace(R1, R2, 5))
ax1.grid(True)

ax2.set_title("évolution de ω1(t)")
ax2.set_ylabel("ω₁(t) (rad.s-1)")
ax2.set_xlabel("t (s)")
ax2.grid(True)

ax3.set_title("évolution de ω2(t)")
ax3.set_ylabel("ω₂(t) (rad.s-1)")
ax3.set_xlabel("t (s)")
ax3.grid(True)

#fonction couleur qui renvoie du rouge si la vitesse est très négative, grise si nulle et bleu si très positive 
#(utilisation de tanh pratique car entre 0 et 1)
def couleur(v):
    r = '#%02x' % int(min(220*(1 - tanh(v/(omega_amplitude*R1*0.35))), 220))
    g = '%02x' % int(max(0,220*(1 - tanh(abs(v)/(omega_amplitude*R1*0.35))))) 
    b = '%02x' % int(min(220*(1 + tanh(v/(omega_amplitude*R1*0.35))), 220))
    return r+g+b

def omega_points():
        return visc.v * visc.r_grid_1d

def affiche_temps(t,dt):
    text_t.set_text('temps: '+str(round(t,3))+' s')
    text_dt.set_text('pas de temps: '+str(round(dt*1e3,3))+' ms')

# si le point de vue choisi est eulerien, on ajoute une ligne de champ initiale, puis on rajoute jusqu'à 12 lignes qui
# suivent chaque point
def eulerien(points, lignes_champ, theta):
    r = visc.r_grid_1d[1:]
    theta += visc.omega()[:-4]*visc.dt
    for i in range(len(points)):
            points[i][0].set_data([r[i]*cos(theta[i])], [r[i]*sin(theta[i])])
            points[i][0].set_color(couleur(visc.v[i]))

            if len(lignes_champ[i]) == 0:
                lignes_champ[i].append(ax.plot([r[i]*cos(theta[i]), r[i]*cos(theta[i])], [r[i]*sin(theta[i]), r[i]*sin(theta[i])], color = couleur(visc.v[i])))
            else:
                x, y = lignes_champ[i][-1][0].get_data()
                x, y = x[-1], y[-1]
                if (x - r[i]*cos(theta[i]))**2 + (y - r[i]*sin(theta[i]))**2 > l_p**2:
                    lignes_champ[i].append(ax.plot([x, r[i]*cos(theta[i])], [y, r[i]*sin(theta[i])], color = couleur(visc.v[i])))
                if len(lignes_champ[i]) > 12:
                    lignes_champ[i][0][0].remove() 
                    lignes_champ[i].pop(0)

# si le point de vue choisi est lagrangien, on enlève les anciens vecteurs et on rajoute ceux au temps t+dt
def lagrangien(vecteurs):
    for objet in ax.collections[:]: 
        if isinstance(objet, Quiver):
            objet.remove()
    R = np.linspace(visc.R1, visc.R2, 10)
    V = np.interp(R, np.linspace(visc.R1, visc.R2, visc.nr), visc.v)
    for i in range(len(R)):
        vecteurs.append(ax.quiver(0, R[i], -visc.R2*V[i]/(2*10*visc.R1), 0, angles='xy', scale_units='xy', scale = 1, color = couleur(V[i])))
        vecteurs.append(ax.quiver(0, -R[i], visc.R2*V[i]/(2*10*visc.R1), 0, angles='xy', scale_units='xy', scale = 1, color = couleur(V[i])))

# fonction appelée a chaque pas de temps dt par l'animation
def f(frames):
    if lancement:
        global theta, dr, theta_1, theta_2, lignes_champ, points, nr, nu
        
        # mise à jour des paramètres à chaque dt
        visc.update()
        affiche_temps(visc.t,visc.dt)
        if boite_nu.text != '':
            nu = float(boite_nu.text)*1e-5
        
        if nu != 0:
            visc.dt = min(visc.r_diff*visc.dr**2/nu, 2*pi/300)
        
        R1, R2 = visc.R1, visc.R2

        # pour éviter les périodes trop longues (voir infinies..)
        if omega_frequence < 1/2:
            T = 2
        else:
            T = 1/omega_frequence

        # calcul des vitesses angulaires et angles des cylindres
        omega_1, omega_2 = visc.current_omega1(),visc.current_omega2()

        theta_1 += omega_1*visc.dt
        theta_2 += omega_2*visc.dt
        
        cercle_int.set_radius(visc.R1)
        cercle_ext.set_radius(visc.R2)
        
        # mise à jour des positions des repères 
        repere_1.set_data([R1*cos(theta_1), (R1 + (R2-R1)/8)*cos(theta_1)],[R1*sin(theta_1), (R1 + (R2-R1)/8)*sin(theta_1)])
        repere_2.set_data([-R1*cos(-theta_1), -(R1 + (R2-R1)/8)*cos(-theta_1)],[R1*sin(-theta_1), (R1 + (R2-R1)/8)*sin(-theta_1)])
        repere_3.set_data([R2*cos(theta_2), (R2 + (R2-R1)/8)*cos(theta_2)],[R2*sin(theta_2), (R2 + (R2-R1)/8)*sin(theta_2)])
        repere_4.set_data([-R2*cos(-theta_2), -(R2 + (R2-R1)/8)*cos(-theta_2)],[R2*sin(-theta_2), (R2 + (R2-R1)/8)*sin(-theta_2)])
        
        # mise a jour du visuel eulerien/lagrangien
        if check_visuel.get_status()[0]:
            eulerien(points, lignes_champ, theta)
        elif check_visuel.get_status()[1]:
            lagrangien(vecteurs)
            
        # mise à jour des graphes, recentrage de l'axe x après un certain temps (4*T)
        graph_omega(ax2, omega_1, visc.t, visc.dt, T, "#35b041")
        graph_omega(ax3, omega_2, visc.t, visc.dt, T, "#c03d3d")
        graph_v(ax1, visc.v, visc.r_grid_1d, R1, R2, "#2b8eff")
    
        if visc.t%(4*T) < visc.dt:
            ax2.set_xlim(visc.t,visc.t + 4*T)
            ax3.set_xlim(visc.t,visc.t + 4*T)


ani = animation.FuncAnimation(fig, f, frames = itertools.count(), interval = 1, blit = False)
plt.show()
