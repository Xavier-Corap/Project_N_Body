import os
import numpy as np

# Fonction pour lire les positions à partir du fichier
def read(file_path):
    data = np.loadtxt(file_path)
    return data

# Répertoire où sont stockés les fichiers
# data_directory = "./"

# Ouvrir un fichier pour écrire les distances
fichier = open("separation.dat", "w")

# Nombre total d'itérations
total_iterations = 4000

# Pas de temps
dt = 3 * 1e-4

# Liste pour stocker les distances et les temps
dist = []
time = []

# Afficher chaque particule à chaque pas de temps
for pas in range(0, total_iterations + 1, 25):
    if pas == 0:
        pas = 1
    
    file_path = "vitesses/" + f"vit_0000{pas}.dat"
    if pas > 9:
        file_path = "vitesses/" + f"vit_000{pas}.dat"
    if pas > 99:
        file_path = "vitesses/" + f"vit_00{pas}.dat"
    if pas > 999:
        file_path = "vitesses/" + f"vit_0{pas}.dat"
        
    file_path2 = "images/" + f"pos_0000{pas}.dat"
    if pas > 9:
        file_path2 = "images/" + f"pos_000{pas}.dat"
    if pas > 99:
        file_path2 = "images/" + f"pos_00{pas}.dat"
    if pas > 999:
        file_path2 = "images/" + f"pos_0{pas}.dat"

    try:
        pos = read(file_path2)
    
        # Calcul de la distance entre deux particules
        galx1 = np.mean(pos[1:,0])
        galy1 = np.mean(pos[1:,1])
        galz1 = np.mean(pos[1:,2])
        
        galx2 = pos[0, 0]
        galy2 = pos[0, 1]
        galz2 = pos[0, 2]
        d = np.sqrt((galx1-galx2)**2 + (galy1-galy2)**2)
        dist.append(round(d/1e20, 2))
        
        time.append(dt*pas)

        # Écriture de la distance dans le fichier
        fichier.write(str(round(d/1e20, 2)) + " ")

    except FileNotFoundError:
        print(f"Le fichier {file_path} n'a pas été trouvé. Arrêt de l'affichage.")
        break

# Fermeture du fichier
fichier.close()

