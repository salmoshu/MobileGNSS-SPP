import transform
import qpsolver

import numpy as np
import pandas as pd
import scipy.sparse

DT_X = 1.0

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
    A = np.empty(shape=(2*(N_x-1), 2*N_x), dtype=object)
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
    C = np.empty(shape=(2*N_y, 2*N_x), dtype=object)
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

    # 多普勒速度相关参数
    x_index  = np.floor(TIME_d / dt).astype(int)
    alpha    = (TIME_d / dt) - x_index
    coeff_y0 = 1 - 3*alpha**2 + 2*alpha**3
    coeff_y1 =     3*alpha**2 - 2*alpha**3
    coeff_v0 = dt * alpha * (alpha - 1)**2
    coeff_v1 = dt * alpha**2 * (alpha - 1)
    C = np.empty(shape=(N_d, N_x), dtype=object)
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

if __name__ == "__main__":
    # 生成测试数据
    np.random.seed(42)
    N = 100  # 数据点数量
    
    # 基础位置数据
    base_df = pd.DataFrame({
        'millisSinceGpsEpoch': np.arange(N) * 1000,
        'Time': np.arange(N) * 1.0,
        'latDeg': 35.68 + np.cumsum(np.random.normal(0, 0.0001, N)),
        'lngDeg': 139.76 + np.cumsum(np.random.normal(0, 0.0001, N))
    })
    
    # 速度数据
    velocity_df = pd.DataFrame({
        'millisSinceGpsEpoch': np.arange(N) * 1000,
        'Time': np.arange(N) * 1.0,
        'v_east': np.random.normal(0, 0.5, N),
        'v_north': np.random.normal(0, 0.5, N)
    })
    
    # 传感器数据
    sensor_df = None
    
    # 优化参数
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
                        'use_map'  : False,
                        'threshold_distance_to_nearest_neighbor' : 8.0,
                        'sigma_p_stage2' : 3.0,
                        'num_stage2_iterations' : 1,
                       }
    ppdf = apply_costmin(base_df, velocity_df, sensor_df, params_downtown, N_LOOP=3)
    print(ppdf.head())