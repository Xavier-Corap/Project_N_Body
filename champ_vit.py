import numpy as np
import matplotlib.pyplot as plt

# Chargement des données de positions et de vitesses depuis les fichiers txt
positions = np.loadtxt('images/pos_01200.dat')
vitesses = np.loadtxt('vitesses/vit_01200.dat')
print(positions.shape)
# Définir les paramètres de la grille 2D
Nx, Ny = 32, 32 # Nombre de cellules en x et y
xmin, xmax = -0*3.08e16, 30000*3.08e16  # Limites en x
ymin, ymax = -0*3.08e16, 30000*3.08e16  # Limites en y
dx = (xmax - xmin) / Nx  # Taille de la cellule en x
dy = (ymax - ymin) / Ny  # Taille de la cellule en y

# Initialisation des tableaux pour stocker les vitesses moyennes
vx_moyen = np.zeros((Nx, Ny))
vy_moyen = np.zeros((Nx, Ny))
max_speed = np.zeros((Nx, Ny))  # Tableau pour stocker la vitesse maximale dans chaque cellule
nb_part = np.ones((Nx, Ny))
# Calcul des vitesses moyennes dans chaque cellule de la grille
for i in range(len(positions)):
    # Trouver les indices de la cellule correspondante
    ix = int((positions[i, 0] - xmin) % (xmax - xmin) / dx)
    iy = int((positions[i, 1] - ymin) % (ymax - ymin) / dy)
    
    # Ajouter la vitesse de la particule à la cellule correspondante
    vx_moyen[ix, iy] += vitesses[i, 0]
    vy_moyen[ix, iy] += vitesses[i, 1]
    
    # Mettre à jour la vitesse maximale dans la cellule
    max_speed[ix, iy] = max(max_speed[ix, iy], np.sqrt(vitesses[i, 0]**2 + vitesses[i, 1]**2))
    if i == 0 :
      print((ix, iy))
      nb_part[ix, iy] = 100000
    nb_part[ix, iy] +=1

print(nb_part)
# Attribuer un poids aux vitesses moyennes
for i in range(Nx):
    for j in range(Ny):
        if max_speed[i, j] != 0:
            vx_moyen[i, j] *= (nb_part[i,j]/1e10)
            vy_moyen[i, j] *= (nb_part[i,j]/1e10)

# Créer une carte de chaleur des vitesses moyennes
plt.style.use('dark_background')
plt.figure(figsize=(8, 6))

density = np.loadtxt("images_density/pos_01200.dat")
density2 = np.loadtxt("images_density/pos_01200.dat")

plt.imshow(density, norm = 'log', extent =(xmin, xmax, ymin, ymax), vmin =1e-23, vmax = 1e-17, origin = "lower")
cax = [1e-16,1e-18,1e-20,1e-22,1e-24, 1e-50]
plt.colorbar(label='Densité intégrée', ticks = cax).minorticks_on()

plt.xlabel('Axe Y')
plt.ylabel('Axe X')



# Ajouter les vecteurs de champ de vitesse
xv, yv = np.meshgrid(np.linspace(xmin, xmax, Nx), np.linspace(ymin, ymax, Ny), indexing='xy')
plt.quiver(xv, yv, vx_moyen, vy_moyen, color='white', scale=50)
plt.title('Carte des vitesses moyennes avec vecteurs (normalisés par la vitesse maximale)')
#plt.gca().invert_yaxis()
#plt.gca().invert_yaxis()
plt.show()

