from dataclasses import dataclass
import numpy as np
import pandas as pd
from scipy.interpolate import InterpolatedUnivariateSpline

import constants as C
import design_filter

@dataclass
class Trig:
    cos : np.array
    sin : np.array

    def __add__(self, other):
        cos = (self.cos * other.cos) - (self.sin * other.sin)
        sin = (self.sin * other.cos) + (self.cos * other.sin)
        return Trig(cos=cos, sin=sin)

    def __sub__(self, other):
        cos = (self.cos * other.cos) + (self.sin * other.sin)
        sin = (self.sin * other.cos) - (self.cos * other.sin)
        return Trig(cos=cos, sin=sin)

    def __neg__(self):
        return Trig(cos=self.cos, sin=-self.sin)

    def to_rad(self):
        return np.arctan2(self.sin, self.cos)

    @staticmethod
    def from_data(x, y):
        r = np.sqrt(x**2 + y**2)
        cos = x / r
        sin = y / r
        return Trig(cos=cos, sin=sin)

    @staticmethod
    def from_rad(th):
        cos = np.cos(th)
        sin = np.sin(th)
        return Trig(cos=cos, sin=sin)

def apply_filter_padding_edge(X, flt):
    assert(len(flt) % 2 == 1)
    n_pad = (len(flt) - 1) // 2
    X_pad = np.concatenate([np.full(n_pad, X[0]), X, np.full(n_pad, X[-1])], axis=0)
    return np.convolve(X_pad, flt, mode='valid')

def preprocess_sensor_data(sensor_dfs, t_ref, dt_up, dt_down, flt):
    assert(len(flt) % 2 == 1)
    acce_df = sensor_dfs['acce']
    gyro_df = sensor_dfs['gyro']
    magn_df = sensor_dfs['magn']
    acce_time = 1e-3 * (acce_df['millisSinceGpsEpoch'] - t_ref).values
    gyro_time = 1e-3 * (gyro_df['millisSinceGpsEpoch'] - t_ref).values
    magn_time = 1e-3 * (magn_df['millisSinceGpsEpoch'] - t_ref).values
    t_up_min  = np.ceil( np.max([acce_time[ 0], gyro_time[ 0], magn_time[ 0]]) / dt_up) * dt_up
    t_up_max  = np.floor(np.min([acce_time[-1], gyro_time[-1], magn_time[-1]]) / dt_up) * dt_up
    time_up   = np.arange(t_up_min, t_up_max + dt_up / 10, dt_up)
    
    acce_x_up = InterpolatedUnivariateSpline(acce_time, acce_df['UncalAccelXMps2'].values, k=1)(time_up)
    acce_y_up = InterpolatedUnivariateSpline(acce_time, acce_df['UncalAccelYMps2'].values, k=1)(time_up)
    acce_z_up = InterpolatedUnivariateSpline(acce_time, acce_df['UncalAccelZMps2'].values, k=1)(time_up)
    gyro_x_up = InterpolatedUnivariateSpline(gyro_time, gyro_df['UncalGyroXRadPerSec'].values, k=1)(time_up)
    gyro_y_up = InterpolatedUnivariateSpline(gyro_time, gyro_df['UncalGyroYRadPerSec'].values, k=1)(time_up)
    gyro_z_up = InterpolatedUnivariateSpline(gyro_time, gyro_df['UncalGyroZRadPerSec'].values, k=1)(time_up)
    magn_x_up = InterpolatedUnivariateSpline(magn_time, magn_df['UncalMagXMicroT'].values, k=1)(time_up)
    magn_y_up = InterpolatedUnivariateSpline(magn_time, magn_df['UncalMagYMicroT'].values, k=1)(time_up)
    magn_z_up = InterpolatedUnivariateSpline(magn_time, magn_df['UncalMagZMicroT'].values, k=1)(time_up)
    sensor_flt = np.stack([
        np.convolve(acce_x_up, flt, mode='valid'),
        np.convolve(acce_y_up, flt, mode='valid'),
        np.convolve(acce_z_up, flt, mode='valid'),
        np.convolve(gyro_x_up, flt, mode='valid'),
        np.convolve(gyro_y_up, flt, mode='valid'),
        np.convolve(gyro_z_up, flt, mode='valid'),
        np.convolve(magn_x_up, flt, mode='valid'),
        np.convolve(magn_y_up, flt, mode='valid'),
        np.convolve(magn_z_up, flt, mode='valid'),
    ], axis=1)
    n_flt = (len(flt) - 1) // 2
    time_flt   = time_up[n_flt:-n_flt]

    roundint = lambda x : int(np.round(x))
    
    dt_ratio    = roundint(dt_down / dt_up)
    t_down_min  = np.ceil( time_flt[ 0] / dt_down) * dt_down
    t_down_max  = np.floor(time_flt[-1] / dt_down) * dt_down
    N_down      = roundint((t_down_max - t_down_min) / dt_down) + 1
    idx_offset  = roundint((t_down_min - time_flt[ 0]) / dt_up)
    idx_down    = dt_ratio * np.arange(N_down) + idx_offset
    time_down   = time_flt[idx_down]
    sensor_down = sensor_flt[idx_down, :]

    columns = ['UncalAccelXMps2', 'UncalAccelYMps2', 'UncalAccelZMps2',
               'UncalGyroXRadPerSec', 'UncalGyroYRadPerSec', 'UncalGyroZRadPerSec',
               'UncalMagXMicroT', 'UncalMagYMicroT', 'UncalMagZMicroT']
    df = pd.DataFrame(sensor_down, columns=columns)
    df['Time'] = time_down
    return df

def calibrate_magn_offset(sensor_df):
    N = sensor_df.shape[0]
    x = sensor_df['UncalMagXMicroT']
    y = sensor_df['UncalMagYMicroT']
    z = sensor_df['UncalMagZMicroT']
    x2 = x**2
    y2 = y**2
    z2 = z**2
    r2 = x2 + y2 + z2

    a11 = 2 * np.sum(x2)
    a12 = 2 * np.sum(x * y)
    a13 = 2 * np.sum(x * z)
    a14 = 2 * np.sum(x)
    a21 = a12
    a22 = 2 * np.sum(y2)
    a23 = 2 * np.sum(y * z)
    a24 = 2 * np.sum(y)
    a31 = a13
    a32 = a23
    a33 = 2 * np.sum(z2)
    a34 = 2 * np.sum(z)
    a41 = a14
    a42 = a24
    a43 = a34
    a44 = 2 * N
    
    b1 = np.sum(x * r2)
    b2 = np.sum(y * r2)
    b3 = np.sum(z * r2)
    b4 = np.sum(r2)

    A = np.array([[a11, a12, a13, a14],
                  [a21, a22, a23, a24],
                  [a31, a32, a33, a34],
                  [a41, a42, a43, a44]])
    b = np.array([b1, b2, b3, b4])
    mx0, my0, mz0, d = np.linalg.solve(A, b)
    mr = np.sqrt(mx0**2 + my0**2 + mz0**2 + 2*d)
    return mx0, my0, mz0, mr

def remove_gyro_drift(omega, threshold):
    mask  = (np.abs(omega) < threshold).astype(float)
    drift = np.sum(omega * mask) / np.sum(mask)
    return omega - drift

def trapezoidal_integration(V, dt):
    dX = (0.5 * dt) * (V[0:-1] + V[1:])
    return np.concatenate([[0], np.cumsum(dX)], axis=0)

def central_difference(X, dt):
    V = (1 / (2*dt)) * (X[2:] - X[0:-2])
    return np.concatenate([[0], V, [0]], axis=0)

def add_calibrated_signals(sensor_df, dt_down):
    mx0, my0, mz0, mr = calibrate_magn_offset(sensor_df)
    MX = sensor_df['UncalMagXMicroT'].values - mx0
    MZ = sensor_df['UncalMagZMicroT'].values - mz0
    trig_th_hat = Trig.from_data(-MX, -MZ)

    omega = sensor_df['UncalGyroYRadPerSec'].values
    omega = remove_gyro_drift(omega, 0.02)
    omega = remove_gyro_drift(omega, 0.01)
    integ_omega = trapezoidal_integration(omega, dt_down)
    trig_integ_omega = Trig.from_rad(integ_omega)

    trig_th_drift      = trig_th_hat - trig_integ_omega
    trig_th_drift_mean = Trig.from_data(np.mean(trig_th_drift.cos), np.mean(trig_th_drift.sin))
    trig_th_residual   = trig_th_drift - trig_th_drift_mean
    th_drift_mean      = trig_th_drift_mean.to_rad()
    FLT = design_filter.make_sinc_filter(F_cutoff=1/15.0, dt=dt_down)
    th_residual = apply_filter_padding_edge(trig_th_residual.to_rad(), FLT)
    th = integ_omega + th_residual + (th_drift_mean - C.MAGNETIC_DECLINATION)
    sensor_df['omega']  = omega
    sensor_df['theta']  = th
    sensor_df['cos_th'] = np.cos(th)
    sensor_df['sin_th'] = np.sin(th)
    sensor_df['dotV']   = - (sensor_df[f'UncalAccelZMps2'] - sensor_df[f'UncalAccelZMps2'].mean())
    return sensor_df

def check_sensor_availability(sensor_dfs):
    if sensor_dfs['acce'] is None:
        return False
    if sensor_dfs['gyro'] is None:
        return False
    if sensor_dfs['magn'] is None:
        return False
    acce_df = sensor_dfs['acce']
    gyro_df = sensor_dfs['gyro']
    magn_df = sensor_dfs['magn']
    if acce_df.shape[0] == 0:
        return False
    if gyro_df.shape[0] == 0:
        return False
    if magn_df.shape[0] == 0:
        return False
    XYZ = ['X', 'Y', 'Z']
    acce = acce_df[[f'UncalAccel{axis}Mps2'     for axis in XYZ]].values
    gyro = gyro_df[[f'UncalGyro{axis}RadPerSec' for axis in XYZ]].values
    magn = magn_df[[f'UncalMag{axis}MicroT'     for axis in XYZ]].values
    if np.sum(np.abs(np.diff(acce, axis=0))) == 0:
        return False
    if np.sum(np.abs(np.diff(acce, axis=0))) == 0:
        return False
    if np.sum(np.abs(np.diff(acce, axis=0))) == 0:
        return False
    acce_dt = 1e-3 * np.diff(acce_df['millisSinceGpsEpoch'].values)
    gyro_dt = 1e-3 * np.diff(acce_df['millisSinceGpsEpoch'].values)
    magn_dt = 1e-3 * np.diff(acce_df['millisSinceGpsEpoch'].values)
    if np.std(acce_dt) > 1e-3:
        return False
    if np.std(gyro_dt) > 1e-3:
        return False
    if np.std(magn_dt) > 1e-3:
        return False
    acce_mean     = np.mean(acce, axis=0)
    expected_acce = np.array([0, 9.8, 0])
    cos_angle_num = np.dot(acce_mean, expected_acce)
    cos_angle_den = np.linalg.norm(acce_mean) * np.linalg.norm(expected_acce)
    angle = np.arccos(cos_angle_num / cos_angle_den)
    if angle > np.deg2rad(20):
        return False
    return True

def remove_different_posture(sensor_df):
    ay = sensor_df['UncalAccelYMps2'] # gravitational acceleration
    valid_orig = np.abs(ay - np.mean(ay)) < 5
    valid = valid_orig
    # Delete the previous and next 5 seconds also.
    for i in range(1, 6):
        valid_shift_plus  = np.concatenate([np.full(i, True), valid_orig[0:-i]], axis=0)
        valid_shift_minus = np.concatenate([valid_orig[i:], np.full(i, True)], axis=0)
        valid = valid & valid_shift_plus & valid_shift_minus
    sensor_df = sensor_df[valid].copy()
    return sensor_df