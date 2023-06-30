import matplotlib.pyplot as plt
import math
import numpy as np
#Variables générales
T = 16  # Nombre de générations

# Variables liées à l'espèce
N0 = 1  # Effectif initial
R = 1.4  # Taux de croissance maximal
MU = 0.1  # Taux de mortalité

# Variables liées au site
K = 200  # Capacité d'accueil

A = np.array([[1, 2], [2, 1]])

def generation(n) :
    """ Fonction de croissance de deux populations en compétition sur une génération

    Cette fonction calcule les effectifs à la génération suivante de deux populations suivant une croissance logistique, 
    en compétition de Lotka-Volterra

    Paramètres
    ----------
    n : list
        Effectifs des populations au temps t
    
    Return
    ------
    list
        Effectifs des populations au temps t+1
    
    """
    n1 = n[0]  # Effectif de la population 1 au temps t
    n2 = n[1]  # Effectif de la population 2 au temps t
    n1_plus = int((n1 + n1*R*(1-((A[0][0]*n1+A[0][1]*n2)/(K*0.8)))))  # Effectif de la population 1 au temps t+1
    n2_plus = int((n2 + n2*R*(1-((A[1][1]*n2+A[1][0]*n1)/K))))  # Effectif de la population 2 au temps t+1

    return [n1_plus,n2_plus]

N01 = 2
N02 = 1

list_n1 = [N01]  # Liste des effectifs de la population 1 à chaque pas de temps
list_n2 = [N02]  # Liste des effectifs de la population 2 à chaque pas de temps

for gen in range(T):
    n_plus = generation([list_n1[gen],list_n2[gen]])
    list_n1 = list_n1 + [n_plus[0]]
    list_n2 = list_n2 + [n_plus[1]]

N01_2 = 1
N02_2 = 2

list_n1_2 = [N01_2]  # Liste des effectifs de la population 1 à chaque pas de temps
list_n2_2 = [N02_2]  # Liste des effectifs de la population 2 à chaque pas de temps

for gen in range(T):
    n_plus = generation([list_n1_2[gen],list_n2_2[gen]])
    list_n1_2 = list_n1_2 + [n_plus[0]]
    list_n2_2 = list_n2_2 + [n_plus[1]]

fig, ax = plt.subplots(1, 2, sharey = 'all')
ax[0].set_xlabel("Génération")
ax[1].set_xlabel("Génération")
ax[0].set_ylabel("Effectifs")
ax[1].set_ylabel("Effectifs")
ax[0].set_title("A.")
ax[1].set_title("B.")

ax[0].plot(list_n1, 'b')
ax[0].plot(list_n2, 'r')
ax[0].plot([6,6],[0,200], 'grey', linestyle = "dotted")
ax[0].plot([13,13],[0,200], 'grey', linestyle = "dotted")
ax[0].annotate(str('Phase 1'), xy=(0.5,190), xycoords='data')
ax[0].annotate(str('Phase 2'), xy=(7.5,190), xycoords='data')
ax[0].annotate(str('Phase 3'), xy=(13.5,190), xycoords='data')

ax[1].plot(list_n1_2, 'b')
ax[1].plot(list_n2_2, 'r')
ax[1].plot([5,5],[0,200], 'grey', linestyle = "dotted")
ax[1].plot([8,8],[0,200], 'grey', linestyle = "dotted")
ax[1].annotate(str('Phase 1'), xy=(0.5,190), xycoords='data')
ax[1].annotate(str('Phase 2'), xy=(5.5,190), xycoords='data')
ax[1].annotate(str('Phase 3'), xy=(10.5,190), xycoords='data')
plt.show()