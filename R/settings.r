dir_exe_BaM <- "/home/famendezrios/Documents/Git/BaM/makefile/"
# MAGE_executable <- "/home/famendezrios/Documents/Softwares/pamhyr2/mage8/mage"
MAGE_executable <- "/home/famendezrios/Documents/Git/mage/exe/mage"


# Setting jump standard deviation for MCMC sampling
jump_MCMC_theta_param_user <- 8
jump_MCMC_error_model_user <- 0.001
threshold_jump_MCMC_error_model <- 0.5

# Calibration
command_line_MAGE <- ""
nCycles <- 100
nAdapt <- 100
burn <- 0.5
nSlim <- 10
ID_model_BaM <- "MAGE_ZQV"
nX_BaM <- 4
nY_BaM <- 5
