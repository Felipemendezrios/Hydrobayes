#!/bin/bash
/home/famendezrios/Documents/Git/BaM_dev/makefile/BaM -cf /home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Development/New_features_Ben/Config_BaM.txt --seed 1  > run1.log 2>&1 &
/home/famendezrios/Documents/Git/BaM_dev/makefile/BaM -cf /home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/Development/New_features_Ben_spatial/Config_BaM.txt --seed 2 > run2.log 2>&1 &
wait

