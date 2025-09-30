import os
import numpy as np
import matplotlib.pyplot as plt
import imageio

# Fonction pour lire les positions à partir du fichier
def read_positions(file_path):
    data = np.loadtxt(file_path)
    return data

# Répertoire où sont stockés les fichiers
# data_directory = "./"

# Nombre total d'itérations
total_iterations = 4000
axis_bounds = (0, 30000 * 3.08e16)
G = 6.67e-11

# Liste pour stocker les noms de fichiers temporaires
temp_files = []

# Créer un dossier temporaire pour les images
temp_dir = 'temp_images_density'
os.makedirs(temp_dir, exist_ok=True)

# Afficher chaque particule à chaque pas de temps
for pas in range(0, total_iterations + 1, 25):
    if pas == 0:
       pas = 1

    file_path = "images_density/" + f"pos_0000{pas}.dat"
    if pas > 9:
        file_path = "images_density/" + f"pos_000{pas}.dat"
    if pas > 99:
        file_path = "images_density/" + f"pos_00{pas}.dat"
    if pas > 999:
        file_path = "images_density/" + f"pos_0{pas}.dat"
        
    file_path2 = "images/" + f"pos_0000{pas}.dat"
    if pas > 9:
        file_path2 = "images/" + f"pos_000{pas}.dat"
    if pas > 99:
        file_path2 = "images/" + f"pos_00{pas}.dat"
    if pas > 999:
        file_path2 = "images/" + f"pos_0{pas}.dat"
    try:
        plt.figure()
        pos = read_positions(file_path2)
        rho = np.zeros(300)
        r = np.linspace(0, 299, 300)
        # Créer une carte 2D de densité
        plt.figure()
        
        galx1 = np.mean(pos[1:,0])
        galy1 = np.mean(pos[1:,1])
        galz1 = np.mean(pos[1:,2])
        
        for i in range(1, len(pos[:-1,0])):
            j = int((np.sqrt((galx1-pos[i,0])**2 + (galy1-pos[i,1])**2 + (galz1-pos[i,2])**2))/1e20)
            if j > 299:
               break
            a = j
            
            if j == 0:
               a = 1
            rho[j] += 2e30*1e5/(((4/3)*np.pi*G*(a*1000*3.08e16)**3)-((4/3)*np.pi*G*((a-1)*1000*3.08e16)**3))


        plt.loglog(r, rho, label="numérique")
        #plt.loglog(r[1:]+mean, rho_theo(r[1:], M), label="théorique")


        plt.xlabel('r (kpc)')
        plt.ylabel(r'$\rho$ (kg.m$^{-3}$)')

        plt.title(r'$\rho$ en fonction de r au pas ' + str(pas))

        
        # Enregistrer le plot comme une image temporaire
        temp_file = os.path.join(temp_dir, f"temp_images_rho{pas:04d}.png")
        temp_files.append(temp_file)
        plt.savefig(temp_file)
        plt.close()

    except FileNotFoundError:
        print(f"Le fichier {file_path} n'a pas été trouvé. Arrêt de l'affichage.")
        break

# Créer un GIF à partir des images temporelles
gif_path = 'output_animation_rho.gif'
imageio.mimsave(gif_path, [imageio.imread(img) for img in temp_files], duration=0.0001)

# Supprimer les images temporaires
for temp_file in temp_files:
    os.remove(temp_file)

# Supprimer le dossier temporaire
os.rmdir(temp_dir)

# Afficher le chemin du GIF généré
print(f"GIF généré avec succès : {gif_path}")

