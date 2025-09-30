import os
import numpy as np
import matplotlib.pyplot as plt
import imageio

# Fonction pour lire les positions à partir du fichier
def read_positions(file_path):
    data = np.loadtxt(file_path)
    return data

# Nombre total d'itérations
total_iterations = 4000

# Liste pour stocker les noms de fichiers temporaires
temp_files = []

# Créer un dossier temporaire pour les images
temp_dir = 'temp_images_density'
os.makedirs(temp_dir, exist_ok=True)

# Afficher chaque particule à chaque pas de temps
for pas in range(0, total_iterations + 1, 25):
    if pas == 0:
       pas = 1

    file_path = f"images_density/pos_{pas:05d}.dat"

    try:
        plt.figure()
        positions = read_positions(file_path)

        # Créer une carte 2D de densité
        plt.imshow(positions, norm='log', extent=(0, 30000 * 3.08e16, 0, 30000 * 3.08e16), vmin=1e-23, vmax=1e-17)
        cax = [1e-16, 1e-18, 1e-20, 1e-22, 1e-24, 1e-50]
        plt.colorbar(label='Densité intégrée', ticks=cax).minorticks_on()
        plt.gca().invert_yaxis()
        plt.xlabel('Axe Y')
        plt.ylabel('Axe X')
        plt.title(f'Carte de densité intégrée selon l\'axe Z au pas {pas}')

        # Enregistrer le plot comme une image temporaire
        temp_file = os.path.join(temp_dir, f"temp_images_density{pas:05d}.png")
        temp_files.append(temp_file)
        plt.savefig(temp_file)
        plt.close()

    except FileNotFoundError:
        print(f"Le fichier {file_path} n'a pas été trouvé. Arrêt de l'affichage.")
        break

# Créer un GIF à partir des images temporelles
gif_path = 'output_animation_density.gif'
imageio.mimsave(gif_path, [imageio.imread(img) for img in temp_files], duration=0.0001)

# Supprimer les images temporaires
for temp_file in temp_files:
    os.remove(temp_file)

# Supprimer le dossier temporaire
os.rmdir(temp_dir)

# Afficher le chemin du GIF généré
print(f"GIF généré avec succès : {gif_path}")

