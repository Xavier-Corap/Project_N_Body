import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import imageio
from scipy.stats import kde

# Fonction pour lire les positions à partir du fichier
def read_positions(file_path):
    data = np.loadtxt(file_path)
    return data

# Répertoire où sont stockés les fichiers
#data_directory = "./"

# Nombre total d'itérations
total_iterations = 4000
axis_bounds = (0, 30000*3.08e16)

# Liste pour stocker les noms de fichiers temporaires
temp_files = []

# Créer un dossier temporaire pour les images
temp_dir = 'temp_images'
os.makedirs(temp_dir, exist_ok=True)

# Afficher chaque particule à chaque pas de temps
for pas in range(0, total_iterations + 2, 25):
    if pas == 0:
       pas = 1
       

    plt.style.use('dark_background')
        

    file_path = "images/"+f"pos_0000{pas}.dat"
    if pas > 9 : 
      file_path = "images/"+f"pos_000{pas}.dat"
    if pas > 99 : 
      file_path = "images/"+f"pos_00{pas}.dat"
    if pas > 999 : 
      file_path = "images/"+f"pos_0{pas}.dat"
      
    try:
        positions = read_positions(file_path)

        # Créer un plot 3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        ax.xaxis.set_pane_color((1,1,1,0))
        ax.yaxis.set_pane_color((1,1,1,0))
        ax.zaxis.set_pane_color((1,1,1,0))
           

        # Afficher les positions des particules
        ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], marker='o', label=f'Pas de temps {pas}', color = 'white', s = 0.05)
        ax.scatter(positions[0, 0], positions[0, 1], positions[0, 2], marker='o', color = 'red', s = 100)
        # Ajouter des labels et un titre
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(f'Positions des particules au pas de temps {pas}')

        # Afficher la légende
        ax.legend()
        

        # Définir les bornes des axes
        ax.set_xlim(axis_bounds)
        ax.set_ylim(axis_bounds)
        ax.set_zlim(axis_bounds)

        # Enregistrer le plot comme une image temporaire
        temp_file = os.path.join(temp_dir, f"temp_image_{pas:04d}.png")
        temp_files.append(temp_file)
        plt.savefig(temp_file)
        plt.close()
        

    except FileNotFoundError:
        print(f"Le fichier {file_path} n'a pas été trouvé. Arrêt de l'affichage.")
        break

# Créer un GIF à partir des images temporelles
gif_path = 'output_animation.gif'
imageio.mimsave(gif_path, [imageio.imread(img) for img in temp_files], duration=0.5)

# Supprimer les images temporaires
for temp_file in temp_files:
    os.remove(temp_file)

# Supprimer le dossier temporaire
os.rmdir(temp_dir)



# Afficher le chemin du GIF généré
print(f"GIF généré avec succès : {gif_path}")

