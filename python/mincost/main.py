import numpy as np
import pandas as pd
import scipy.sparse
import scipy.sparse.linalg
import multiprocessing
import glob
from tqdm.notebook import tqdm

import io_f
import signal_f
import design_filter
import transform
import qpsolver
import map_matching
import area_prediction

INPUT_PATH    = '../input/google-smartphone-decimeter-challenge'
VELOCITY_PATH = '../input/gsdc-vehicle-speed-estimation-result/_doppler_velocity'
DT_X = 1.0

VELOCITY_DF = pd.concat([pd.read_csv(f'{VELOCITY_PATH}/doppler_velocity_train.csv'),
                         pd.read_csv(f'{VELOCITY_PATH}/doppler_velocity_test.csv'),
                         ], axis=0)

def get_optimization_constants(base_df, velocity_df, sensor_df, params, use_sensor):
    const = dict()
    dt    = DT_X
    TIME_y = base_df['Time'].values
    TIME_d = velocity_df['Time'].values
    N_y = TIME_y.shape[0]
    N_d = TIME_d.shape[0]
    N_x = int(np.ceil(np.max(TIME_y) / dt) + 1)
    const['N_y'] = N_y
    const['N_d'] = N_d
    const['N_x'] = N_x

    a = np.array([[1, dt, (1/2)*dt**2],
                  [0,  1,  dt],
                  [0,  0,  1]])
    e3 = scipy.sparse.eye(3)
    A = np.empty(shape=(2*(N_x-1), 2*N_x), dtype=np.object)
    for i_x in range(N_x-1):
        A[2*i_x  , 2*i_x  ] = a
        A[2*i_x+1, 2*i_x+1] = a
        A[2*i_x  , 2*i_x+2] = -e3
        A[2*i_x+1, 2*i_x+3] = -e3
    const['A'] = scipy.sparse.bmat(A, format='csr')
    
    b = np.array([[(1/6)*dt**3,
                   (1/2)*dt**2,
                   dt]]).T
    const['B'] = scipy.sparse.block_diag([b for _ in range(2*(N_x-1))], format='csr')

    diag_R  = np.full(2*N_x - 2, params['sigma_u']**(-2) * dt)
    const['R'] = scipy.sparse.spdiags(diag_R, [0], 2*N_x - 2, 2*N_x - 2, format='csc')
    
    x_index  = np.floor(TIME_y / dt).astype(int)
    alpha    = (TIME_y / dt) - x_index
    coeff_y0 = 1 - 3*alpha**2 + 2*alpha**3
    coeff_y1 =     3*alpha**2 - 2*alpha**3
    coeff_v0 = dt * alpha * (alpha - 1)**2
    coeff_v1 = dt * alpha**2 * (alpha - 1)
    C = np.empty(shape=(2*N_y, 2*N_x), dtype=np.object)
    for i_x in range(N_x):
        C[0, 2*i_x  ] = scipy.sparse.coo_matrix((1, 3))
        C[0, 2*i_x+1] = scipy.sparse.coo_matrix((1, 3))
    for i_y in range(N_y):
        i_x = x_index[i_y]
        c_i = np.array([[coeff_y0[i_y], coeff_v0[i_y], 0]])
        C[2*i_y,   2*i_x]   = c_i
        C[2*i_y+1, 2*i_x+1] = c_i
        if i_x < N_x - 1:
            c_iplus = np.array([[coeff_y1[i_y], coeff_v1[i_y], 0]])
            C[2*i_y,   2*i_x+2] = c_iplus
            C[2*i_y+1, 2*i_x+3] = c_iplus
    const['Cp_orig']  = scipy.sparse.bmat(C, format='csr')

    diag_Lp = np.full(2*N_y, params['sigma_p']**(-2))
    const['Lp_orig'] = scipy.sparse.spdiags(diag_Lp, [0], 2*N_y, 2*N_y, format='csr')
    const['Yp_orig'] = base_df[['latDeg', 'lngDeg']].values.flatten()

    BLH = transform.BLH(
        lat=np.deg2rad(base_df['latDeg'].values),
        lng=np.deg2rad(base_df['lngDeg'].values),
        hgt=np.zeros(N_y),
    )
    DEG2RAD = np.pi / 180.0
    J = transform.jacobian_BL_to_EN(BLH) * DEG2RAD
    J = np.mean(J, axis=0)
    J[0, 0] = 0
    J[1, 1] = 0
    JJ = scipy.sparse.block_diag([J, J], format='csr')
    const['J'] = J

    # ドップラ速度に関するパラメータ
    x_index  = np.floor(TIME_d / dt).astype(int)
    alpha    = (TIME_d / dt) - x_index
    coeff_y0 = 1 - 3*alpha**2 + 2*alpha**3
    coeff_y1 =     3*alpha**2 - 2*alpha**3
    coeff_v0 = dt * alpha * (alpha - 1)**2
    coeff_v1 = dt * alpha**2 * (alpha - 1)
    C = np.empty(shape=(N_d, N_x), dtype=np.object)
    for i_x in range(N_x):
        C[0, i_x] = scipy.sparse.coo_matrix((2, 6))
    for i_d in range(N_d):
        i_x = x_index[i_d]
        c = np.array([[0, coeff_y0[i_d], coeff_v0[i_d], 0, 0, 0],
                      [0, 0, 0, 0, coeff_y0[i_d], coeff_v0[i_d]]])
        C[i_d, i_x] = J @ c
        if i_x < N_x - 1:
            c = np.array([[0, coeff_y1[i_d], coeff_v1[i_d], 0, 0, 0],
                          [0, 0, 0, 0, coeff_y1[i_d], coeff_v1[i_d]]])
            C[i_d, i_x+1] = J @ c
    const['Cd_orig']  = scipy.sparse.bmat(C, format='csr')

    diag_Ld = np.full(2*N_d, params['sigma_d']**(-2))
    const['Ld_orig'] = scipy.sparse.spdiags(diag_Ld, [0], 2*N_d, 2*N_d, format='csr')
    const['Yd_orig'] = velocity_df[['v_east', 'v_north']].values.flatten()

    if sensor_df is None:
        const['use_sensor'] = False
        const['use_inquality'] = False
        return const

    TIME_s = sensor_df['Time'].values
    N_s = TIME_s.shape[0]
    const['N_s'] = N_s
    const['use_sensor'] = use_sensor
    const['use_inquality'] = (use_sensor and params['use_not_go_back_constraint'])
    x_index = np.round(TIME_s / dt).astype(int)
    const['x_index_sensor'] = x_index
    if not use_sensor:
        return const

    # 速度制約・速度コストに関するパラメータ
    COS_TH = sensor_df['cos_th'].values
    SIN_TH = sensor_df['sin_th'].values
    CV = np.empty(shape=(N_s, N_x), dtype=np.object)
    GV = np.empty(shape=(N_s, N_x), dtype=np.object)
    cv = np.array([[0, 1, 0, 0, 0, 0],
                   [0, 0, 0, 0, 1, 0]], dtype=np.float64)
    for i_x in range(N_x):
        CV[0, i_x] = scipy.sparse.coo_matrix((1, 6))
        GV[0, i_x] = scipy.sparse.coo_matrix((1, 6))
    for i_s in range(N_s):
        i_x = x_index[i_s]
        k = np.array([[SIN_TH[i_s], -COS_TH[i_s]]])
        CV[i_s, i_x] = k @ J @ cv
        k = np.array([[-COS_TH[i_s], -SIN_TH[i_s]]])
        GV[i_s, i_x] = k @ J @ cv
    const['Cv'] = scipy.sparse.bmat(CV, format='csr')
    const['Gv'] = scipy.sparse.bmat(GV, format='csr')
    const['hv'] = np.full((N_s, ), -params['vmin'])

    diag_Lv = np.full(N_s, params['sigma_v']**(-2))
    const['Lv'] = scipy.sparse.spdiags(diag_Lv, [0], N_s, N_s, format='csr')

    # 加速度コストに関するパラメータ
    DOT_V_COS_TH = sensor_df['dotV'] * sensor_df['cos_th'].values
    DOT_V_SIN_TH = sensor_df['dotV'] * sensor_df['sin_th'].values
    OMEGA = sensor_df['omega'].values
    CA = np.empty(shape=(N_s, N_x), dtype=np.object)
    ca = np.array([[0, 1, 0, 0, 0, 0],
                   [0, 0, 0, 0, 1, 0],
                   [0, 0, 1, 0, 0, 0],
                   [0, 0, 0, 0, 0, 1]], dtype=np.float64)
    for i_x in range(N_x):
        CA[0, i_x] = scipy.sparse.coo_matrix((2, 6))
    for i_s in range(N_s):
        i_x = x_index[i_s]
        k = np.array([[0,  OMEGA[i_s], 1, 0],
                      [-OMEGA[i_s], 0, 0, 1]])
        CA[i_s, i_x] = k @ JJ @ ca
    const['Ca'] = scipy.sparse.bmat(CA, format='csr')
    const['Ya'] = np.stack([DOT_V_COS_TH, DOT_V_SIN_TH], axis=1).flatten()

    diag_La = np.full(2*N_s, params['sigma_a']**(-2))
    const['La'] = scipy.sparse.spdiags(diag_La, [0], 2*N_s, 2*N_s, format='csr')

    return const

def solve_QP(const, p_valid, d_valid):
    A = const['A']
    B = const['B']
    R = const['R']
    Cp_orig = const['Cp_orig']
    Lp_orig = const['Lp_orig']
    Yp_orig = const['Yp_orig']
    Cd_orig = const['Cd_orig']
    Ld_orig = const['Ld_orig']
    Yd_orig = const['Yd_orig']
    
    p_valid2 = np.stack([p_valid, p_valid], axis=1).flatten()
    Cp = Cp_orig[p_valid2, :]
    Lp = Lp_orig[np.ix_(p_valid2, p_valid2)]
    Yp = Yp_orig[p_valid2]

    d_valid2 = np.stack([d_valid, d_valid], axis=1).flatten()
    Cd = Cd_orig[d_valid2, :]
    Ld = Ld_orig[np.ix_(d_valid2, d_valid2)]
    Yd = Yd_orig[d_valid2]

    CLC_p = Cp.T @ (Lp @ Cp)
    CLC_d = Cd.T @ (Ld @ Cd)
    
    CLY_p = Cp.T @ (Lp @ Yp)
    CLY_d = Cd.T @ (Ld @ Yd)
    
    if const['use_sensor']:
        Cv = const['Cv']
        Lv = const['Lv']
        Ca = const['Ca']
        La = const['La']
        Ya = const['Ya']
        CLC_v = Cv.T @ (Lv @ Cv)
        CLC_a = Ca.T @ (La @ Ca)
        Q     = CLC_p + CLC_d + CLC_v + CLC_a

        CLY_a = Ca.T @ (La @ Ya)
        q     = CLY_p + CLY_d + CLY_a
    else:
        Q = CLC_p + CLC_d
        q = CLY_p + CLY_d

    if const['use_inquality']:
        G = const['Gv']
        h = const['hv']
        X_star = qpsolver.solve_qp_with_inequality(R=R, Q=Q, q=q, A=A, B=B, G=G, h=h)
    else:
        X_star = qpsolver.solve_qp(R=R, Q=Q, q=q, A=A, B=B)
    return X_star

def get_baseline(collection_name):
    df = BASELINE_DF[BASELINE_DF['collectionName'] == collection_name].copy()
    df.reset_index(drop=True, inplace=True)
    return df

def get_velocity(collection_name):
    df = VELOCITY_DF[VELOCITY_DF['collectionName'] == collection_name].copy()
    df.reset_index(drop=True, inplace=True)
    return df

def apply_costmin(base_df, velocity_df, sensor_df, params, N_LOOP):
    const = get_optimization_constants(base_df, velocity_df, sensor_df, params, use_sensor=True)

    if params['use_map']:
        distance = map_matching.distance_to_nearest_neighbor(base_df)
        default_p_valid = (distance < params['threshold_distance_to_nearest_neighbor'])
        p_valid = default_p_valid
    else:
        default_p_valid = np.full(const['N_y'], True)
        p_valid = default_p_valid

    V = np.sqrt(np.sum(velocity_df[['v_east', 'v_north']].values**2, axis=1))
    default_d_valid = (V < params['vmax'])
    d_valid = default_d_valid

    for loop in range(N_LOOP):
        X_star = solve_QP(const, p_valid, d_valid)
        Y_star = const['Cp_orig'] @ X_star
        Y_star = np.reshape(Y_star, (-1, 2))
        pp_df  = base_df.copy()
        pp_df['latDeg'] = Y_star[:, 0]
        pp_df['lngDeg'] = Y_star[:, 1]

        # 计算位置误差并更新有效点
        distance = transform.pd_haversine_distance(pp_df, base_df)
        p_valid = default_p_valid & (distance < params['reject_p'])

        # 计算速度误差并更新有效点
        dXYdt = const['Cd_orig'] @ X_star
        dXYdt = np.reshape(dXYdt, (-1, 2))
        v_err = dXYdt - velocity_df[['v_east', 'v_north']].values
        v_err = np.sqrt(np.sum(v_err**2, axis=1))
        d_valid = default_d_valid & (v_err < params['reject_d'])

    return pp_df

def recalibrate_sensor_by_vehicle_motion(base_df, velocity_df, sensor_df, params):
    const = get_optimization_constants(base_df, velocity_df, sensor_df, params, use_sensor=False)

    if params['use_map']:
        distance = map_matching.distance_to_nearest_neighbor(base_df)
        p_valid  = (distance < params['threshold_distance_to_nearest_neighbor'])
    else:
        p_valid = np.full(const['N_y'], True)

    V = np.sqrt(np.sum(velocity_df[['v_east', 'v_north']].values**2, axis=1))
    d_valid = (V < params['vmax'])

    X_star = solve_QP(const, p_valid, d_valid)
    X_mat  = np.reshape(X_star, (-1, 6))
    dotB   = X_mat[const['x_index_sensor'], 1]
    dotL   = X_mat[const['x_index_sensor'], 4]
    dotXY  = const['J'] @ np.stack([dotB, dotL], axis=0) # shape = (2, N)
    dotX   = dotXY[0, :]
    dotY   = dotXY[1, :]
    V = np.sqrt(np.sum(dotXY**2, axis=0))
    cond = (V > (20 / 3.6))
    trig_moving_direction = signal_f.Trig.from_data(dotX[cond], dotY[cond])
    trig_theta   = signal_f.Trig.from_rad(sensor_df['theta'].values[cond])
    trig_offset  = trig_theta - trig_moving_direction
    angle_offset = np.arctan2(np.mean(trig_offset.sin), np.mean(trig_offset.cos))
    sensor_df['theta']  = sensor_df['theta'] - angle_offset
    sensor_df['cos_th'] = np.cos(sensor_df['theta'].values)
    sensor_df['sin_th'] = np.sin(sensor_df['theta'].values)

    return sensor_df

def do_postprocess(args):
    train_or_test, collection, params = args

    base_df = get_baseline(collection)
    t_ref   = base_df['millisSinceGpsEpoch'].min()
    base_df['Time'] = 1e-3 * (base_df['millisSinceGpsEpoch'] - t_ref).values

    velocity_df = get_velocity(collection)
    velocity_df['Time'] = (1e-3 * (velocity_df['millisSinceGpsEpoch'] - t_ref).values
                           -  params['Mi8_velocity_timeshift'] * (velocity_df['phoneName'] == 'Mi8').astype(float)
                           )
    velocity_df = velocity_df[(  velocity_df['Time'] >= base_df['Time'].min())
                              & (velocity_df['Time'] <= base_df['Time'].max())]
    velocity_df.reset_index(drop=True, inplace=True)
    
    phone_list = [path.split('/')[-1] for path in sorted(glob.glob(f'{INPUT_PATH}/{train_or_test}/{collection}/*'))]
    sensor_df_list   = []
    dt_up   = 2.5 * 1e-3
    dt_down = DT_X
    FLT = design_filter.make_sinc_filter(F_cutoff=2.0, dt=dt_up)
    for phone in phone_list:
        gnss_log_filename = f'{INPUT_PATH}/{train_or_test}/{collection}/{phone}/{phone}_GnssLog.txt'
        sensor_df_orig = io_f.read_GnssLog_sensors(gnss_log_filename)
        if signal_f.check_sensor_availability(sensor_df_orig):
            sensor_df = signal_f.preprocess_sensor_data(sensor_df_orig, t_ref, dt_up, dt_down, FLT)
            sensor_df = signal_f.remove_different_posture(sensor_df)
            sensor_df = sensor_df[(  sensor_df['Time'] >= base_df['Time'].min())
                                  & (sensor_df['Time'] <= base_df['Time'].max())].copy()
            sensor_df.reset_index(drop=True, inplace=True)
            sensor_df_list.append(sensor_df)
    if len(sensor_df_list) > 0:
        time_list = [df['Time'].max() - df['Time'].min() for df in sensor_df_list]
        idx = np.argmax(time_list)
        sensor_df = sensor_df_list[idx]
        sensor_df = signal_f.add_calibrated_signals(sensor_df, dt_down)
        sensor_df = recalibrate_sensor_by_vehicle_motion(base_df, velocity_df, sensor_df, params)
    else:
        sensor_df = None

    pp_df = base_df
    pp_df = apply_costmin(pp_df, velocity_df, sensor_df, params, N_LOOP=3)
    if params['use_map']:
        params_stage2 = dict(params)
        params_stage2['sigma_p'] = params['sigma_p_stage2']
        for _ in range(params['num_stage2_iterations']):
            pp_df = map_matching.snap_to_nearest_neighbor(pp_df)
            pp_df = apply_costmin(pp_df, velocity_df, sensor_df, params_stage2, N_LOOP=1)
    return pp_df

def make_postprocessing_df(train_or_test, config):
    args_list = []
    for collection_list, params in config:
        for collection_name in collection_list:
            args_list.append((train_or_test, collection_name, params))
    processes = multiprocessing.cpu_count()
    with multiprocessing.Pool(processes=processes) as pool:
        df_list = pool.imap_unordered(do_postprocess, args_list)
        df_list = tqdm(df_list, total=len(args_list))
        df_list = list(df_list)
    output_df = pd.concat(df_list, axis=0).sort_values(['phone', 'millisSinceGpsEpoch'])
    return output_df

def print_score(output_df):
    score_list = []
    for gid, phone_df in output_df.groupby('phone'):
        drive, phone = gid.split('_')
        gt_df = pd.read_csv(f'{INPUT_PATH}/train/{drive}/{phone}/ground_truth.csv')
        d = transform.pd_haversine_distance(phone_df, gt_df)
        score = np.mean([np.quantile(d, 0.50), np.quantile(d, 0.95)])
        score_list.append(score)
    score = np.mean(score_list)
    print(f'train score: {score:.3f}')
    return

def main():
    params_highway = { 'sigma_u'  : 1.0,
                       'sigma_p'  : 3.0,
                       'sigma_a'  : 2.0 * 1e+5,
                       'sigma_v'  : 4.0 * 1e+5,
                       'sigma_d'  : 0.16 * 1e+5,
                       'reject_p' : 7.0,   # [m]
                       'reject_d' : 1.0,   # [m/s]
                       'vmin'     : -0.05, # [m/s]
                       'vmax'     : 50.0,  # [m/s]
                       'Mi8_velocity_timeshift' : 0.46,
                       'use_not_go_back_constraint' : False,
                       'use_map'  : False,
                      }
    params_treeway = { 'sigma_u'  : 1.0,
                       'sigma_p'  : 6.0,
                       'sigma_a'  : 0.8 * 1e+5,
                       'sigma_v'  : 3.0 * 1e+5,
                       'sigma_d'  : 0.12 * 1e+5,
                       'reject_p' : 12.0,  # [m]
                       'reject_d' : 1.0,   # [m/s]
                       'vmin'     : -0.05, # [m/s]
                       'vmax'     : 50.0,  # [m/s]
                       'Mi8_velocity_timeshift' : 0.30,
                       'use_not_go_back_constraint' : False,
                       'use_map'  : False,
                      }
    params_downtown = { 'sigma_u'  : 1.0,
                        'sigma_p'  : 20.0,
                        'sigma_a'  : 0.4 * 1e+5,
                        'sigma_v'  : 1.0 * 1e+5,
                        'sigma_d'  : 1.3 * 1e+5,
                        'reject_p' : 20.0,  # [m]
                        'reject_d' : 3.0,   # [m/s]
                        'vmin'     : -0.05, # [m/s]
                        'vmax'     : 50.0,  # [m/s]
                        'Mi8_velocity_timeshift' : 0.0,
                        'use_not_go_back_constraint' : False,
                        'use_map'  : True,
                        'threshold_distance_to_nearest_neighbor' : 8.0,
                        'sigma_p_stage2' : 3.0,
                        'num_stage2_iterations' : 1,
                       }
    collection_list_all = np.array(sorted(path.split('/')[-1] for path in glob.glob(f'{INPUT_PATH}/train/*')))
    collection_list_highway  = collection_list_all[np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]) - 1]
    collection_list_treeway  = collection_list_all[np.array([22,23,25,26,28]) - 1]
    collection_list_downtown = collection_list_all[np.array([24,27,29]) - 1]
    config = [
        (collection_list_highway,  params_highway),
        (collection_list_treeway,  params_treeway),
        (collection_list_downtown, params_downtown),
    ]
    train_pp_df = make_postprocessing_df('train', config)
    print_score(train_pp_df)

    test_base = pd.read_csv(f'{INPUT_PATH}/baseline_locations_test.csv')
    collection_list_highway, collection_list_treeway, collection_list_downtown = area_prediction.predict_area(test_base)
    config = [
        (collection_list_highway,  params_highway),
        (collection_list_treeway,  params_treeway),
        (collection_list_downtown, params_downtown),
    ]
    test_pp_df = make_postprocessing_df('test', config)

    train_pp_df.to_csv('smoothing_2nd_train.csv', index=False)
    test_pp_df.to_csv('smoothing_2nd_test.csv', index=False)

    columns = ['phone', 'millisSinceGpsEpoch', 'latDeg', 'lngDeg']
    sub_df = test_pp_df[columns]
    sub_df.to_csv('submission.csv', index=False)
    return