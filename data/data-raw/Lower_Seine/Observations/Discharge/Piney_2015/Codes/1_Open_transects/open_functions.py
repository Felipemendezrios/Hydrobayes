# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 14:28:25 2022

@author: blais
"""
# ========================================
# External imports
# ========================================
import os
import glob
import tkinter as tk
import tkinter.filedialog
import numpy as np
import scipy.io as sio
import warnings

# ========================================
# Internal imports
# ========================================
from Classes.Measurement import Measurement


# =============================================================================
# Functions
# =============================================================================
def select_file():
    # Open a window to select a measurement
    root = tk.Tk()
    root.withdraw()
    root.attributes("-topmost", True)
    path_meas = tk.filedialog.askopenfilenames(parent=root, title='Select ADCP file',
                                               filetypes=[('ADCP data', '.mat .mmt')])
    if path_meas:
        keys = [elem.split('.')[-1].split('.')[0] for elem in path_meas]
        name_folders = [path_meas[0].split('/')[-2].split('.')[0]][0]

        if 'mmt' in keys:
            type_meas = 'TRDI'
            file_path = path_meas[0]
        elif 'mat' in keys:
            qrev_data = False
            for path in path_meas:
                mat_data = sio.loadmat(path, struct_as_record=False, squeeze_me=True)
                if 'version' in mat_data:
                    type_meas = 'QRev'
                    file_path = mat_data
                    qrev_data = True
                    break
            if not qrev_data:
                type_meas = 'SonTek'
                ind_meas = [i for i, s in enumerate(keys) if "mat" in s]
                file_path = [(path_meas[x]) for x in ind_meas]

        return file_path, type_meas, name_folders
    else:
        warnings.warn('No file selected - end')
        return None, None, None


def select_directory():
    # Open a window to select a folder which contains measurements
    root = tk.Tk()
    root.withdraw()
    root.attributes("-topmost", True)
    path_folder = tk.filedialog.askdirectory(parent=root, title='Select folder')
    if path_folder:
        # ADCP folders path
        path_folder = '\\'.join(path_folder.split('/'))
        # path_folders = np.array(glob.glob(path_folder + "/*"))
        path_folders = [f.path for f in os.scandir(path_folder) if f.is_dir()]

        # Load their name
        name_folders = np.array([os.path.basename((x)) for x in path_folders])
        # Exclude files
        # excluded_folders = [s.find('.') == -1 for s in name_folders]
        # path_folders = path_folders[excluded_folders]
        # name_folders = name_folders[excluded_folders]

        # Open measurement
        type_meas = list()
        path_meas = list()
        name_meas = list()
        no_adcp = list()
        for id_meas in range(len(path_folders)):
            list_files = os.listdir(path_folders[id_meas])
            exte_files = list([i.split('.', 1)[-1] for i in list_files])
            if 'mmt' in exte_files or 'mat' in exte_files:
                if 'mat' in exte_files:
                    loca_meas = [i for i, s in enumerate(exte_files) if "mat" in s]  # transects index
                    fileNameRaw = [(list_files[x]) for x in loca_meas]
                    qrev_data = False
                    path_qrev_meas = []
                    for name in fileNameRaw:
                        path = os.path.join(path_folders[id_meas], name)
                        mat_data = sio.loadmat(path, struct_as_record=False, squeeze_me=True)
                        if 'version' in mat_data:
                            path_qrev_meas.append(os.path.join(path_folders[id_meas], name))
                            qrev_data = True
                    if qrev_data:
                        recent_time = None
                        recent_id = 0
                        for i in range(len(path_qrev_meas)):
                            qrev_time = os.path.getmtime(path_qrev_meas[i])
                            if recent_time is None or qrev_time < recent_time:
                                recent_time = qrev_time
                                recent_id = i
                        path = path_qrev_meas[recent_id]
                        mat_data = sio.loadmat(path, struct_as_record=False, squeeze_me=True)
                        print('QRev file')
                        type_meas.append('QRev')
                        name_meas.append(name_folders[id_meas])
                        path_meas.append(mat_data)

                    else:
                        print('SonTek file')
                        type_meas.append('SonTek')
                        fileName = [s for i, s in enumerate(fileNameRaw)
                                    if "QRev.mat" not in s]
                        name_meas.append(name_folders[id_meas])
                        path_meas.append(
                            [os.path.join(path_folders[id_meas], fileName[x]) for x in
                             range(len(fileName))])

                else:
                    print('TRDI file')
                    type_meas.append('TRDI')
                    loca_meas = exte_files.index("mmt")  # index of the measurement file
                    fileName = list_files[loca_meas]
                    path_meas.append(os.path.join(path_folders[id_meas], fileName))
                    name_meas.append(name_folders[id_meas])
            else:
                no_adcp.append(name_folders[id_meas])
                print(f"No ADCP : {name_folders[id_meas]}")
        return path_folder, path_meas, type_meas, name_meas, no_adcp
    else:
        warnings.warn('No folder selected - end')
        return None, None, None, None, None