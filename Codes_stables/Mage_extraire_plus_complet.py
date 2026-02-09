import os
import pandas as pd
import numpy as np
import sys

# Chemin vers ton fichier texte
path = "/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Case_studies/Synthetic_case/Simplified_ks_rectangular_MC/PamHyr/_PAMHYR_/"

mes_case = "Synt_rec_2"

mes_variables = ["Z"]
mes_derive = "X"
reach_number = 3


# Fonction pour convertir "DDD:HH:MM:SS" → secondes
def time_to_seconds(time_str):
    days, hours, minutes, seconds = map(int, time_str.split(":"))
    return (((days * 24 + hours) * 60 + minutes) * 60) + seconds


# Define the variable to extract:
if str(mes_derive) == "X":
    # if spatial variable (dX), a temporal grid is necessary
    file_path = str(path) + str(mes_case) + "/default-mage/" + str(mes_case) + ".PAR"

    params = {}

    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("init_time"):
                params["init_time"] = line.split()[1]
            elif line.startswith("final_time"):
                params["final_time"] = line.split()[1]
            elif line.startswith("timestep_bin"):
                params["timestep_bin"] = int(line.split()[1])
    # Conversion en secondes
    start = time_to_seconds(params["init_time"])
    end = time_to_seconds(params["final_time"])
    step = params["timestep_bin"]

    # Temporal grid
    grid = np.arange(start, end + 1, step)

elif str(mes_derive) == "T":
    # if temporal variable (dT), a spatial grid is necessary

    file_path = str(path) + "smooth_bed_14_PK.txt"

    # Lire la ligne comme une liste de valeurs séparées par des virgules
    pk_values = pd.read_csv(file_path, header=None, sep=";")
    # On récupère la première (et seule) ligne
    ligne = pk_values.iloc[0]

    # On la convertit en liste Python
    pk_list = ligne.tolist()

    grid = [float(x) for x in pk_list]

# Définition du répertoire de travail
workdir = str(path) + str(mes_case) + "/default-mage/"

os.chdir(workdir)

for variable in mes_variables:
    for id in grid:
        f = open("setting.txt", "w")
        f.writelines(str(mes_case) + "\n")
        f.writelines(variable + "d" + str(mes_derive) + "\n")
        f.writelines(str(reach_number) + "\n")
        if (variable + "d" + str(mes_derive)) == "ZdX":
            # Mage question : Date de T0 au format AAAA:MM:JJTHH:MM:SS ou en secondes
            f.writelines("secondes\n")
            # Zdt : Mage question : Date de la ligne d'eau à extraire en secondes
            f.writelines(str(id) + "\n")
        elif (variable + "d" + str(mes_derive)) == "QdT":
            f.writelines(str(id) + "\n")  # QdT
            f.writelines("0 1\n")  # QdT
        else:
            sys.exit(
                "Variable extraction format is new, be sure to put the get information before to use the script, could be necessary to adapt it"
            )
        f.close()
        os.system(
            "/home/famendezrios/Documents/Softwares/pamhyr2/mage8/mage_Extraire < setting.txt"
        )
        os.rename(
            str(mes_case) + ".res",
            str(mes_case) + variable + str(id) + "_reach_" + str(reach_number) + ".res",
        )
