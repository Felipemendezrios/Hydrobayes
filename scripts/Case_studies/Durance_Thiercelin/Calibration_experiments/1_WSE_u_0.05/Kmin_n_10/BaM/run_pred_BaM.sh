#!/bin/bash
exec > /dev/null 2>&1

nohup '/home/famendezrios/Documents/Git/BaM/makefile/BaM' -cf '/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/scripts/Case_studies/Durance_Thiercelin/Calibration_experiments/1_WSE_u_0.05/Kmin_n_10/BaM/Config_BaM_ParamU.txt' > /dev/null 2>&1 &
nohup '/home/famendezrios/Documents/Git/BaM/makefile/BaM' -cf '/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/scripts/Case_studies/Durance_Thiercelin/Calibration_experiments/1_WSE_u_0.05/Kmin_n_10/BaM/Config_BaM_Maxpost.txt' > /dev/null 2>&1 &
nohup '/home/famendezrios/Documents/Git/BaM/makefile/BaM' -cf '/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/scripts/Case_studies/Durance_Thiercelin/Calibration_experiments/1_WSE_u_0.05/Kmin_n_10/BaM/Config_BaM_TotalU.txt' > /dev/null 2>&1 &
