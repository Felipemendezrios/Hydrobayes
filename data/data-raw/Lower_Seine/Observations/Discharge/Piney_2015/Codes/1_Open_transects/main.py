import os

import pandas as pd

from open_functions import select_file
import numpy as np
from MiscLibs.common_functions import cart2pol, pol2cart
from datetime import datetime
from Classes.Measurement import Measurement
import pandas as pd

navigation_reference_user = 'BT'

# Load path of the measurement
path_meas, type_meas, name_meas = select_file()

if path_meas is not None:
    # Open ADCP measurement
    meas = Measurement(in_file=path_meas, source=type_meas, proc_type='QRev', run_oursin=True, use_weighted=True)

    settings = meas.current_settings()
    settings_change = False

    # Change navigation reference
    if navigation_reference_user == 'BT' and meas.current_settings()['NavRef'] != 'bt_vel':
        settings['NavRef'] = 'BT'
    elif navigation_reference_user == 'GGA':
        settings['NavRef'] = 'GGA'
        settings['ggaDiffQualFilter'] = 1
    meas.apply_settings(settings)

    date_start = []
    date_end = []
    q_list = []
    oursin_u = []
    checked_list = meas.checked_transect_idx
    for id_transect in checked_list:
        print(id_transect)

        # Change selected transects
        checked_transect_idx = [id_transect]
        meas.checked_transect_idx = []
        for n in range(len(meas.transects)):
            if n in checked_transect_idx:
                meas.transects[n].checked = True
                meas.checked_transect_idx.append(n)
            else:
                meas.transects[n].checked = False
        meas.selected_transects_changed(checked_transect_idx)

        # Date data
        date_start.append(datetime.utcfromtimestamp(meas.transects[id_transect].date_time.start_serial_time).strftime('%Y-%m-%d %H:%M:%S'))
        date_end.append(datetime.utcfromtimestamp(meas.transects[id_transect].date_time.end_serial_time).strftime('%Y-%m-%d %H:%M:%S'))

        # Discharge
        q_list.append(meas.mean_discharges(meas)['total_mean'])

        # uncertainty
        oursin_u.append(np.median(meas.oursin.u['total_95']))

    df = pd.DataFrame({'start': date_start, 'end':date_end, 'discharge': q_list, 'oursin': oursin_u})
    if type(path_meas) is list:
        name = path_meas[0]
    else:
        name = path_meas
    df.to_csv(os.getcwd()+os.sep+name.split('/')[-2]+'.csv', sep=';')
    print(os.getcwd()+os.sep+name.split('/')[-2]+'.csv')


    # # Get transect data
    # id_transect = 1
    # transect = meas.transects[id_transect]
    # transect.boat_vel.bt_vel.u_processed_mps
    # depth = transect.depths.bt_depths.depth_processed_m
    # ship_data = transect.boat_vel.compute_boat_track(transect, ref=meas.current_settings()['NavRef'])
    #
    #
    # if transect.orig_start_edge == 'Right':
    #     # Reverse transects in ordred to start at 0 on left edge
    #     # Valid data
    #     valid = transect.depths.bt_depths.valid_data[::-1]
    #     # Track
    #     dmg_ind = np.where(abs(ship_data['dmg_m']) == max(abs(ship_data['dmg_m'])))[0][0]
    #     x_track = ship_data['track_x_m'] - ship_data['track_x_m'][dmg_ind]
    #     y_track = ship_data['track_y_m'] - ship_data['track_y_m'][dmg_ind]
    #     x_transect = x_track[::-1]
    #     y_transect = y_track[::-1]
    #     dmg_transect = ship_data['dmg_m'][::-1]
    #     distance = abs(ship_data['distance_m'][::-1] - max(ship_data['distance_m']))
    #     # Time
    #     timestamp = (np.nancumsum(transect.date_time.ens_duration_sec) + transect.date_time.start_serial_time)[::-1]
    #     # Depth
    #     depth_transect = transect.depths.bt_depths.depth_processed_m[::-1]
    #     cells_depth = transect.depths.bt_depths.depth_cell_depth_m[:, ::-1]
    #     # Velocity data
    #     vel_x = np.copy(transect.w_vel.u_processed_mps[:, ::-1])
    #     vel_y = np.copy(transect.w_vel.v_processed_mps[:, ::-1])
    #     vel_z = np.copy(transect.w_vel.w_mps[:, ::-1])
    #     x_velocity = vel_x[:, valid]
    #     y_velocity = vel_y[:, valid]
    #     z_velocity = vel_z[:, valid]
    # else:
    #     # Valid data
    #     valid = transect.depths.bt_depths.valid_data
    #     # Track
    #     x_transect = ship_data['track_x_m']
    #     y_transect = ship_data['track_y_m']
    #     dmg_transect = ship_data['dmg_m']
    #     distance = ship_data['distance_m'][::-1]
    #     # Time
    #     timestamp = np.nancumsum(transect.date_time.ens_duration_sec) + transect.date_time.start_serial_time
    #     # Depth
    #     depth_transect = transect.depths.bt_depths.depth_processed_m
    #     cells_depth = transect.depths.bt_depths.depth_cell_depth_m
    #     # Velocity data
    #     vel_x = np.copy(transect.w_vel.u_processed_mps)
    #     vel_y = np.copy(transect.w_vel.v_processed_mps)
    #     vel_z = np.copy(transect.w_vel.w_mps)
    #     x_velocity = vel_x[:, valid]
    #     y_velocity = vel_y[:, valid]
    #     z_velocity = vel_z[:, valid]
    #
    # # Convert time tu utc
    # t = []
    # for stamp in timestamp:
    #     t.append(datetime.utcfromtimestamp(stamp))
    #
    # # Edge shape
    # left_dist = meas.transects[id_transect].edges.left.distance_m
    # left_coef = meas.discharge[id_transect].edge_coef('left', meas.transects[id_transect])
    # right_dist = meas.transects[id_transect].edges.right.distance_m
    # right_coef = meas.discharge[id_transect].edge_coef('right', meas.transects[id_transect])
    #
    # # GPS
    # # lat = transect.gps.gga_lat_ens_deg
    # # lon = transect.gps.gga_lon_ens_deg
    # # utm_position = transect.gps.utm_ens_m
    #
    # # Compute mean velocity components in each ensemble
    # w_vel_mean_1 = np.nanmean(x_velocity, 0)
    # w_vel_mean_2 = np.nanmean(y_velocity, 0)
    #
    # # Compute a unit vector
    # direction, _ = cart2pol(w_vel_mean_1, w_vel_mean_2)
    # unit_vec_1, unit_vec_2 = pol2cart(direction, 1)
    # unit_vec = np.vstack([unit_vec_1, unit_vec_2])
    #
    # # Compute the velocity magnitude in the direction of the mean velocity of each
    # # ensemble using the dot product and unit vector
    # w_vel_prim = np.tile([np.nan], x_velocity.shape)
    # w_vel_sec = np.tile([np.nan], x_velocity.shape)
    # for i in range(x_velocity.shape[0]):
    #     w_vel_prim[i, :] = np.sum(np.vstack([x_velocity[i, :], y_velocity[i, :]]) * unit_vec, 0)
    #     w_vel_sec[i, :] = unit_vec_2 * x_velocity[i, :] - unit_vec_1 * y_velocity[i, :]
    #




