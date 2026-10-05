
Codes:

Python: 1_Open_transects: only work in windows because of QRevInt
This code allows estimating the uncertainty of each gauging under the assumption that they are independent.
No mean is calculated because river dynamics cause the discharge to vary over time, and the length of the cross-section is too large to assume that the discharge is constant in time.
That is why the estimation is performed using transects.

Oursin uncertainty is calculated following the BT signal. More information in Oursin_uncertainty_analysis.pdf

This code generate csv files stored here: /home/famendezrios/Documents/papiers/Qmec/Data/Collect_information_Raw_Data/Lower_Seine/Gaugings/Piney_2015/Intermedial_results/

Results of each discharge measurement either performed using M9 or Rio Grande with the associated uncertainty.

R: 2_Oursin_Uncertainty_estimation_more_assign_Gaugings.r

This script let to add in the uncertainty estimation following Oursin_uncertainty_analysis.pdf.
Besides, it let to filter incoherent discharge measurements and correct Time zone of devices

This code generate RData files stored in processed data folder: /home/famendezrios/Documents/papiers/Qmec/Data/Processed_data/Lower_Seine/Gaugings/Piney_2018/

