import os
import subprocess
import shutil
import json
import numpy as np

# Définition des paramètres de simulation
n_gal = 2
N_grille = 128

parametres = {
    "galaxies": n_gal,
    "particules_galaxie_1": 100000,
    "particules_galaxie_2": 870000,
    "distance entre les 2" : 2*5000*3.08e16,
    "N_grille": N_grille,
    "n_iter_max": [2500, 40],
    "a": 0.1 * 15000 * 3.08e16,
    "vitesse galaxie2" : "120km/s",
    "vitesse galaxie1" : "1200km/s"
}

# Fonction pour exécuter un fichier Python
def execute_python_file(file_path):
   try:
      completed_process = subprocess.run(['python3',file_path], capture_output=True, text=True)
      if completed_process.returncode == 0:
         print("Execution successful.")
         print("Output:")
         print(completed_process.stdout)
      else:
         print(f"Error: Failed to execute '{file_path}'.")
         print("Error output:")
         print(completed_process.stderr)
   except FileNotFoundError:
      print(f"Error: The file '{file_path}' does not exist.")

# Exécution des fichiers Python pour générer les données de simulation
file_path = 'vitesses.py'
execute_python_file(file_path)

file_path = 'images.py'
execute_python_file(file_path)

file_path = 'density.py'
execute_python_file(file_path)

file_path = 'rhobis.py'
execute_python_file(file_path)

file_path = 'separation.py'
execute_python_file(file_path)

# Création du dossier de destination pour les résultats
dst_folder = 'résultats/'+str(n_gal)+'galaxie_'+str(N_grille)+'test_31/'
isExist = os.path.exists(dst_folder)
if not isExist:
    os.makedirs(dst_folder, mode = 0o777, exist_ok = True)
        
# Copie des données générées dans le dossier de destination
data = open("separation.dat", "r").read()
src = 'separation.dat'
dst = os.path.join(dst_folder, src) 
fichier = open(dst, "w")
fichier.write(data)

# Déplacement des fichiers d'animation dans le dossier de destination
src = 'output_vitesses.gif'
dst = os.path.join(dst_folder, 'output_vitesses.gif') 
if os.path.exists(dst):
    os.remove(dst)
shutil.move(src, dst)

src = 'output_animation.gif'
dst = os.path.join(dst_folder, 'output_animation.gif') 
if os.path.exists(dst):
    os.remove(dst)
shutil.move(src, dst)

src = 'output_animation_density.gif'
dst = os.path.join(dst_folder, 'output_animation_density.gif')
if os.path.exists(dst):
    os.remove(dst)
shutil.move(src, dst)

src = 'output_animation_rho.gif'
dst = os.path.join(dst_folder, 'output_animation_rho.gif')
if os.path.exists(dst):
    os.remove(dst)
shutil.move(src, dst)

# Enregistrement des paramètres de simulation dans un fichier texte
with open('paramètres.txt', 'w') as param: 
    param.write(json.dumps(parametres))

src = 'paramètres.txt'
dst = os.path.join(dst_folder, 'paramètres.txt')
if os.path.exists(dst):
    os.remove(dst)
shutil.move(src, dst)

param.close()

