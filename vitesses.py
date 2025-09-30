import os
import numpy as np
import matplotlib.pyplot as plt
import imageio

# Définition des constantes
a = 0.05 * 15000 * 3.08e16
G = 6.67e-11
Mp = 2e30 * 1e5
size = 30000 * 3.08e16
xini = 5000 * 3.08e16 * 0
yini = 5000 * 3.08e16 * 0

# Fonction pour lire les positions à partir du fichier
def read(file_path):
    data = np.loadtxt(file_path)
    return data

# Liste pour stocker les noms de fichiers temporaires
temp_files = []

# Nombre total d'itérations
total_iterations = 4000

# Créer un dossier temporaire pour les images
temp_dir = 'temp_vitesses'
os.makedirs(temp_dir, exist_ok=True)

# Afficher chaque particule à chaque pas de temps
for pas in range(0, total_iterations + 1, 25):
    if pas == 0:
       pas = 1

    # Construction des chemins vers les fichiers de positions et de vitesses
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
        # Lecture des données de vitesses et de positions
        v_f = read(file_path)
        pos_f = read(file_path2)

        # Création d'un nouveau graphique
        plt.figure()

        # Calcul des distances et de la masse
        ray_f = np.sqrt((pos_f[:,0] - 0.5 * size + xini)**2 + (pos_f[:,1] - 0.5 * size - yini)**2 + (pos_f[:,2] - 0.5 * size)**2)
        M = len(ray_f) * Mp

        # Plot de la loi de probabilité g(q)
        q = np.linspace(0, 1, 1000)
        plt.plot(q, (1 - q**2)**(7/2) * q**2 / 0.04295, label="loi de probabilité g(q)")

        # Plot de l'histogramme des vitesses
        q = np.sqrt(v_f[:,0]**2 + v_f[:,1]**2 + v_f[:,2]**2) / ((ray_f**2 + a**2)**(-1./4.) * np.sqrt(2 * G * M))
        plt.hist(q, bins=200, density=True, range=(0, 1))
        plt.xlim(0, 1)
        plt.legend()

        # Titre du graphique
        plt.title('Distribution des vitesses au pas ' + str(pas))

        # Enregistrement du graphique comme une image temporaire
        temp_file = os.path.join(temp_dir, f"temp_vitesses{pas:04d}.png")
        temp_files.append(temp_file)
        plt.savefig(temp_file)
        plt.close()

    except FileNotFoundError:
        print(f"Le fichier {file_path} n'a pas été trouvé. Arrêt de l'affichage.")
        break

# Création d'un GIF à partir des images temporaires
gif_path = 'output_vitesses.gif'
imageio.mimsave(gif_path, [imageio.imread(img) for img in temp_files], duration=0.0001)

# Suppression des images temporaires
for temp_file in temp_files:
    os.remove(temp_file)

# Suppression du dossier temporaire
os.rmdir(temp_dir)

# Affichage du chemin du GIF généré
print(f"GIF généré avec succès : {gif_path}")

