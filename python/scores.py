import numpy as np
import pandas as pd
import pymap3d as pm
import pymap3d.vincenty as pmv

import os
import sys
import shutil
from datetime import datetime

current_dir = os.path.dirname(os.path.abspath(__file__))
rtkpkg_dir = os.path.join(current_dir, './rtklibpy')
sys.path.append(rtkpkg_dir)
cfgfile = os.path.join(rtkpkg_dir, 'config_spp.py')
shutil.copyfile(cfgfile, '__ppk_config.py')

import __ppk_config as cfg
import rtkcmn as gn

GPS_EPOCH = datetime(1970, 1, 1)  # GPS 时代的起始时间 (1970-01-01 00:00:00 UTC)

def dms_to_decimal(degrees_minutes, direction):
    """
    将 NMEA 格式的度分 (DDMM.MMMM) 转换为十进制度 (DD.DDDD)。
    
    Args:
        degrees_minutes (str): 度分格式，例如 "3959.89023018"
        direction (str): 方向，"N"/"S" 或 "E"/"W"
    
    Returns:
        float: 十进制度
    """
    if not degrees_minutes:
        return 0.0
    # 分离度和分
    degrees = int(float(degrees_minutes) // 100)  # 取整得到度
    minutes = float(degrees_minutes) % 100  # 取余得到分
    decimal = degrees + minutes / 60.0  # 转换为十进制度
    # 根据方向调整符号
    if direction in ["S", "W"]:
        decimal = -decimal
    return decimal

def get_res_from_nmea(nmea_file):
    """
    从 NMEA 文件中提取 GNSS 数据
    
    Args:
        nmea_file (str): 输入 NMEA 文件路径
    Returns:
        list: 包含时间、纬度、经度、速度、航向的 GT 数据列表
    """
    # 初始化数据列表
    date_obj = None
    gt_data = []
    
    # 临时存储 GGA 和 RMC 数据
    current_gga = None
    current_rmc = None
    
    with open(nmea_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('!'):
                continue
            
            # 解析 NMEA 语句
            parts = line.split(',')
            if len(parts) < 2:
                continue
            
            sentence_type = parts[0]
            # for GnssLogger
            if parts[0] == 'NMEA':
                parts = parts[1:]
                sentence_type = parts[0]

            # 解析 $GNRMC 或 $GPRMC
            if sentence_type in ['$GNRMC', '$GPRMC']:
                if len(parts) < 12:
                    continue
                # 检查有效性 (A 或 D)
                status = parts[2]
                if status == 'V':  # 无效定位
                    continue
                
                # 提取时间、纬度、经度、速度、航向
                time_str = parts[1]  # UTC 时间，格式 HHMMSS.ss
                date_str = parts[9]  # 日期，格式 DDMMYY
                lat = parts[3]  # 纬度，格式 DDMM.MMMM
                lat_dir = parts[4]  # 纬度方向，N/S
                lon = parts[5]  # 经度，格式 DDDMM.MMMM
                lon_dir = parts[6]  # 经度方向，E/W
                speed_knots = parts[7]  # 速度（节）
                bearing = parts[8]  # 航向（度）

                # 转换时间
                if not (time_str and date_str):
                    continue
                try:
                    date_obj = datetime.strptime(f"{date_str} {time_str}", "%d%m%y %H%M%S.%f")
                    unix_time_millis = int((date_obj - GPS_EPOCH).total_seconds() * 1000)
                except ValueError:
                    continue
                
                # 转换纬度和经度
                latitude = dms_to_decimal(lat, lat_dir)
                longitude = dms_to_decimal(lon, lon_dir)
                
                # 转换速度（从节到米/秒，1 节 = 0.514444 m/s）
                speed_mps = float(speed_knots) * 0.514444 if speed_knots else 0.0
                
                # 航向
                bearing_deg = float(bearing) if bearing else 0.0
                
                current_rmc = {
                    'unix_time_millis': unix_time_millis,
                    'latitude': latitude,
                    'longitude': longitude,
                    'speed_mps': speed_mps,
                    'bearing_deg': bearing_deg
                }
            
            # 解析 $GNGGA 或 $GPGGA
            elif sentence_type in ['$GNGGA', '$GPGGA']:
                if len(parts) < 15:
                    continue
                # 检查定位质量
                fix_quality = int(parts[6])
                if fix_quality == 0:  # 无效定位
                    continue

                # 提取时间、高度、卫星数量
                time_str = parts[1]  # UTC 时间，格式 HHMMSS.ss
                lat = parts[2]  # 纬度，格式 DDMM.MMMM
                lat_dir = parts[3]  # 纬度方向，N/S
                lon = parts[4]  # 经度，格式 DDDMM.MMMM
                lon_dir = parts[5]  # 经度方向，E/W
                altitude = parts[9]  # 高度（米）
                geoid_height = parts[11]  # 地球椭球高度（米）
                
                # 转换时间
                if not time_str:
                    continue
                # GGA 语句没有日期，使用 RMC 的日期
                if date_obj:
                    date_str = date_obj.strftime("%d%m%y")
                    try:
                        time_obj = datetime.strptime(f"{date_str} {time_str}", "%d%m%y %H%M%S.%f")
                        unix_time_millis = int((time_obj - GPS_EPOCH).total_seconds() * 1000)
                    except ValueError:
                        continue
                else:
                    continue
                
                # 转换高度（高度 = 海拔高度 + 地球椭球高度）
                altitude_m = float(altitude) if altitude else 0.0
                geoid_height_m = float(geoid_height) if geoid_height else 0.0
                altitude_m += geoid_height_m
                
                current_gga = {
                    'unix_time_millis': unix_time_millis,
                    'altitude_m': altitude_m
                }
            
            # 如果同时有 RMC 和 GGA 数据，生成 GT 记录
            if current_rmc and current_gga and current_rmc['unix_time_millis'] == current_gga['unix_time_millis']:
                gt_record = {
                    'MessageType': 'Fix',
                    'Provider': 'GT',
                    'LatitudeDegrees': round(current_rmc['latitude'], 9),  # 9 位小数
                    'LongitudeDegrees': round(current_rmc['longitude'], 9),  # 9 位小数
                    'AltitudeMeters': round(current_gga['altitude_m'], 4),  # 4 位小数
                    'SpeedMps': round(current_rmc['speed_mps'], 6),  # 6 位小数
                    'AccuracyMeters': 1.2,  # 固定为 1.2 米，1 位小数
                    'BearingDegrees': round(current_rmc['bearing_deg'], 1),  # 1 位小数
                    'UnixTimeMillis': current_rmc['unix_time_millis']  # 整数，无小数
                }
                gt_data.append(gt_record)
                # 重置临时数据
                current_rmc = None
                current_gga = None
    
    # 转换为 DataFrame 并保存
    if gt_data:
        df = pd.DataFrame(gt_data)
        df = df[['MessageType', 'Provider', 'LatitudeDegrees', 'LongitudeDegrees',
                 'AltitudeMeters', 'SpeedMps', 'AccuracyMeters', 'BearingDegrees',
                 'UnixTimeMillis']]
        # 设置显示格式
        pd.options.display.float_format = '{:.9f}'.format  # 统一浮点数显示格式
        print(f"Format {len(gt_data)} records from {nmea_file}")
        return df
    else:
        print("No valid records generated.")
        return None
    
def align_gt_with_utc(gt_df, utc):
    """
    根据 utc 的整秒时间戳，对齐 gt_df，确保两者行数相同。
    utc 是一维 NumPy 数组，表示时间戳（可能非整秒），将其四舍五入到整秒后进行匹配。
    
    Args:
        gt_df (pd.DataFrame): 真值数据，包含 UnixTimeMillis 列（整秒时间戳）
        utc (np.ndarray): 一维时间戳数组，形状为 (n_utc,)，可能非整秒
    
    Returns:
        pd.DataFrame: 裁剪后的 gt_df，行数与 utc 相同
    """
    # 1. 提取时间戳
    gt_timestamps = gt_df['UnixTimeMillis'].astype(np.int64).values  # gt_df 的时间戳
    utc_timestamps = utc.astype(np.int64)  # utc 是一维数组，直接使用

    # 2. 为 utc 新增整秒时间戳：四舍五入到最近的整秒
    utc_rounded_timestamps = (np.round(utc_timestamps / 1000) * 1000).astype(np.int64)

    # 3. 为每个 utc 时间戳找到 gt_df 中最近的时间戳
    matched_indices = []
    for utc_ts in utc_rounded_timestamps:
        # 找到 gt_timestamps 中最近的时间戳
        closest_idx = np.argmin(np.abs(gt_timestamps - utc_ts))
        matched_indices.append(closest_idx)

    # 4. 转换为 NumPy 数组
    matched_indices = np.array(matched_indices)

    # 5. 根据匹配的索引裁剪 gt_df
    aligned_gt_df = gt_df.iloc[matched_indices].copy()

    # 6. 确保行数一致
    if len(aligned_gt_df) != len(utc):
        raise ValueError(f"Row counts do not match: aligned_gt_df has {len(aligned_gt_df)} rows, utc has {len(utc)} rows")

    return aligned_gt_df

'''
pos data example:
%  GPST                      x-ecef(m)      y-ecef(m)      z-ecef(m)   Q  ns   sdx(m)   sdy(m)   sdz(m)  sdxy(m)  sdyz(m)  sdzx(m) age(s)  ratio    vx(m/s)    vy(m/s)    vz(m/s)      sdvx     sdvy     sdvz    sdvxy    sdvyz    sdvzx
2025/03/12 06:41:08.000  -2181703.6705   4379531.0601   4077857.2671   5  21   0.0073   0.0155   0.0103  -0.0098   0.0094  -0.0064   0.00    0.0   -0.40179    0.02403   -0.36281   0.14020  0.19425  0.23790 -0.12020  0.16007 -0.11952
2025/03/12 06:41:09.000  -2181704.1562   4379531.1373   4077849.4817   5  21   0.0145   0.0126   0.0345  -0.0114   0.0201  -0.0187   0.00    0.0   -0.16121    0.03221   -0.09058   0.13883  0.19411  0.23297 -0.12003  0.16043 -0.12339
2025/03/12 06:41:10.000  -2181703.4290   4379530.7210   4077848.6123   5  21   0.0261   0.0169   0.0187  -0.0146   0.0162  -0.0189   0.00    0.0   -0.13202    0.16239    0.11867   0.14127  0.19510  0.23363 -0.12215  0.16151 -0.12525

Construct (all results are numpy matrices):
utc:       Unix timestamp (note that the original data in the pos file is GPST)
position:  Position
velocity:  Velocity
pos_cov:   Position covariance
vel_cov:   Velocity covariance
'''
def get_res_from_pos(file_path):
    utc = []
    position = []
    velocity = []
    pos_cov = []
    vel_cov = []

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('%'):
                continue
            parts = line.split()

            # Convert GPST to UTC
            gpst_str = parts[0] + ' ' + parts[1]
            gpst = datetime.strptime(gpst_str, '%Y/%m/%d %H:%M:%S.%f')
            gpst = gn.epoch2time((gpst.year, gpst.month, gpst.day, gpst.hour, gpst.minute, gpst.second + gpst.microsecond / 1e6))
            utc_time = gn.gpst2utc(gpst)
            utc_timestamp = round((utc_time.time + utc_time.sec)*1000)
            utc.append(utc_timestamp)

            # Position
            position.append([float(parts[2]), float(parts[3]), float(parts[4])])

            # Velocity
            velocity.append([float(parts[15]), float(parts[16]), float(parts[17])])

            # Position covariance
            sdx, sdy, sdz, sdxy, sdyz, sdzx = map(float, parts[7:13])
            # print(gpst_str, sdx, sdy, sdz, sdxy, sdyz, sdzx)
            # raise ''
            # Convert standard deviations to covariances
            cov = np.array([
                [sdx**2, sdxy**2*np.sign(sdxy), sdzx**2*np.sign(sdzx)],
                [sdxy**2*np.sign(sdxy), sdy**2, sdyz**2*np.sign(sdyz)],
                [sdzx**2*np.sign(sdzx), sdyz**2*np.sign(sdyz), sdz**2]
            ])
            pos_cov.append(cov)

            # Velocity covariance
            sdvx, sdvy, sdvz, sdvxy, sdvyz, sdvzx = map(float, parts[18:])
            cov = np.array([
                [sdvx**2, sdvxy**2*np.sign(sdvxy), sdvzx**2*np.sign(sdvzx)],
                [sdvxy**2*np.sign(sdvxy), sdvy**2, sdvyz**2*np.sign(sdvyz)],
                [sdvzx**2*np.sign(sdvzx), sdvyz**2*np.sign(sdvyz), sdvz**2]
            ])
            vel_cov.append(cov)

    utc = np.array(utc)
    position = np.array(position)
    velocity = np.array(velocity)
    pos_cov = np.array(pos_cov)
    vel_cov = np.array(vel_cov)

    return utc, position, pos_cov, velocity, vel_cov

######### Score Computation #########
# Compute distance by Vincenty's formulae
def vincenty_distance(llh1, llh2):
    """
    Args:
        llh1 : [latitude,longitude] (deg)
        llh2 : [latitude,longitude] (deg)
    Returns:
        d : distance between llh1 and llh2 (m)
    """
    d, az = np.array(pmv.vdist(llh1[:, 0], llh1[:, 1], llh2[:, 0], llh2[:, 1]))

    return d


# Compute score
def calc_score(llh, llh_gt):
    """
    Args:
        llh : [latitude,longitude] (deg)
        llh_gt : [latitude,longitude] (deg)
    Returns:
        score : (m)
    """
    d = vincenty_distance(llh, llh_gt)
    score_50 = np.quantile(d, 0.50)
    score_95 = np.quantile(d, 0.95)
    score_100 = np.quantile(d, 1.00)
    score = np.mean([score_50, score_95])

    return score, [score_50, score_95, score_100]

if __name__ == '__main__':
    path = r'../data/01-opensky/data01'
    # path = r'../data/02-street/data01'
    # path = r'../data/03-downtown/data01'
    # path = r'../data/04-elevated/data01'

    obs_path = path + r'\rover.obs'
    nav_path = path + r'\rover.nav'
    scene, group = path.split('/')[-2:]

    utc, x_wls, cov_x, v_wls, cov_v = get_res_from_pos(path + r'\rover.pos') # opitimized result
    bl_df = get_res_from_nmea(path + r'\baseline.nmea') # baseline result
    gt_df = get_res_from_nmea(path + r'\rtk.nmea')      # ground truth result

    # Ground truth and baseline
    aligned_gt_df = align_gt_with_utc(gt_df, utc) # 对齐时间戳
    aligned_bl_df = align_gt_with_utc(bl_df, utc) # 对齐时间戳
    llh_ekf = np.array(pm.ecef2geodetic(x_wls[:, 0], x_wls[:, 1], x_wls[:, 2])).T
    llh_gt = aligned_gt_df[['LatitudeDegrees', 'LongitudeDegrees']].to_numpy()
    llh_bl = aligned_bl_df[['LatitudeDegrees', 'LongitudeDegrees']].to_numpy()

    # Distance from ground truth
    vd_bl = vincenty_distance(llh_bl, llh_gt)
    vd_ekf = vincenty_distance(llh_ekf, llh_gt)

    # Score
    score_bl, scores_bl = calc_score(llh_bl, llh_gt)
    score_ekf, scores_ekf = calc_score(llh_ekf, llh_gt)

    print("\033[1;32m" + f"{scene}, {group} [nepoch = {len(utc)}]: (CEP50+CEP95)/2, CEP50, CEP95, CEP100" + "\033[0m")
    print(f'Score Baseline   {score_bl:.2f} ({scores_bl[0]:.2f} {scores_bl[1]:.2f} {scores_bl[2]:.2f}) [m]')
    print(f'Score EKF       {score_ekf:.2f} ({scores_ekf[0]:.2f} {scores_ekf[1]:.2f} {scores_ekf[2]:.2f}) [m]')
