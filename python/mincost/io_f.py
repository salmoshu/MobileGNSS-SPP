import io
import datetime
from dataclasses import dataclass, asdict
import numpy as np
import pandas as pd
from scipy.interpolate import InterpolatedUnivariateSpline, RectBivariateSpline

UTC_TO_GPS_OFFSET_MS = ((datetime.date(1980, 1, 6) - datetime.date(1970, 1, 1)).days * 24 * 3600 - 18) * 1000

def read_GnssLog_sensors(filename):
    acce_lines = []
    gyro_lines = []
    magn_lines = []
    orie_lines = []
    with open(filename, 'r') as f:
        for line in f:
            if 'UncalAccel' in line:
                line = line.rstrip().lstrip('#')
                acce_lines.append(line)
                continue
            if 'UncalGyro' in line:
                line = line.rstrip().lstrip('#')
                gyro_lines.append(line)
                continue
            if 'UncalMag' in line:
                line = line.rstrip().lstrip('#')
                magn_lines.append(line)
                continue
            if 'OrientationDeg' in line:
                line = line.rstrip().lstrip('#')
                orie_lines.append(line)
                continue
    acce_df = pd.read_csv(io.StringIO('\n'.join(acce_lines))) if len(acce_lines) != 0 else None
    gyro_df = pd.read_csv(io.StringIO('\n'.join(gyro_lines))) if len(gyro_lines) != 0 else None
    magn_df = pd.read_csv(io.StringIO('\n'.join(magn_lines))) if len(magn_lines) != 0 else None
    orie_df = pd.read_csv(io.StringIO('\n'.join(orie_lines))) if len(orie_lines) != 0 else None
    def modify_df(df):
        if df is None:
            return None
        df.dropna(axis=0, inplace=True)
        df.drop_duplicates(subset='utcTimeMillis', inplace=True)
        if df.shape[0] < 2:
            return None
        dt_valid = np.concatenate([[True], np.diff(df['utcTimeMillis'].values) > 0])
        df = df[dt_valid].copy()
        df.reset_index(drop=True, inplace=True)
        df['millisSinceGpsEpoch'] = df['utcTimeMillis'] - UTC_TO_GPS_OFFSET_MS
        return df
    dfs = dict(
        acce = modify_df(acce_df),
        gyro = modify_df(gyro_df),
        magn = modify_df(magn_df),
        orie = modify_df(orie_df),
    )
    return dfs

@dataclass
class IONEX:
    iono_height : float
    base_radius : float
    lat_1       : float
    lat_2       : float
    lat_delta   : float
    lng_1       : float
    lng_2       : float
    lng_delta   : float
    time_1      : np.datetime64
    time_2      : np.datetime64
    time_delta  : np.timedelta64
    iono_map    : np.array
    lat_range   : np.array
    lng_range   : np.array

def concat_sp3(sp3_df_list):
    sp3_df  = pd.concat(sp3_df_list, axis=0)
    sat_set_list = [frozenset(sp3_df['SatName']) for sp3_df in sp3_df_list]
    sat_sum  = sat_set_list[0]
    sat_prod = sat_set_list[0]
    for sat_set in sat_set_list[1:]:
        sat_sum  = sat_sum  | sat_set
        sat_prod = sat_prod & sat_set
    sat_partial = sat_sum - sat_prod
    for sat in sat_partial:
        # print(f'concat_sp3: drop {sat}')
        sp3_df = sp3_df[sp3_df['SatName'] != sat]
    sp3_df = sp3_df.reset_index(drop=True)
    return sp3_df

def concat_ionex(ionex_list):
    assert(len(np.unique([ionex.iono_height for ionex in ionex_list])) == 1)
    assert(len(np.unique([ionex.base_radius for ionex in ionex_list])) == 1)
    assert(len(np.unique([ionex.lat_1       for ionex in ionex_list])) == 1)
    assert(len(np.unique([ionex.lat_2       for ionex in ionex_list])) == 1)
    assert(len(np.unique([ionex.lat_delta   for ionex in ionex_list])) == 1)
    assert(len(np.unique([ionex.lng_1       for ionex in ionex_list])) == 1)
    assert(len(np.unique([ionex.lng_2       for ionex in ionex_list])) == 1)
    assert(len(np.unique([ionex.lng_delta   for ionex in ionex_list])) == 1)
    assert(len(np.unique([ionex.time_delta  for ionex in ionex_list])) == 1)
    N = len(ionex_list)
    iono_map = []
    for i in range(N-1):
        assert(ionex_list[i].time_2 == ionex_list[i+1].time_1)
        iono_map.append(ionex_list[i].iono_map[0:-1, :, :])
    iono_map.append(ionex_list[-1].iono_map)
    kw = asdict(ionex_list[0])
    kw['time_2']   = ionex_list[-1].time_2
    kw['iono_map'] = np.concatenate(iono_map, axis=0)
    return IONEX(**kw)

def read_GnssLog_Raw(filename):
    lines = []
    with open(filename, 'r') as f:
        for line in f:
            if 'Raw' in line:
                line = line.rstrip().lstrip('#')
                lines.append(line)
    sio = io.StringIO('\n'.join(lines))
    return pd.read_csv(sio)

def read_clock_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    for index, line in enumerate(lines):
        if 'TIME SYSTEM ID' in line:
            assert(line.strip().split()[0] == 'GPS')
            continue
        if 'END OF HEADER' in line:
            start_index = index + 1
            break
    lines = lines[start_index:]
    SAT, EPOCH, DELTA_TSV = [], [], []
    for line in lines:
        if not line.startswith('AS '):
            continue
        tokens = line.rstrip().split()
        sat = tokens[1]
        epoch = datetime.datetime(year   = int(tokens[2]),
                                  month  = int(tokens[3]),
                                  day    = int(tokens[4]),
                                  hour   = int(tokens[5]),
                                  minute = int(tokens[6]),
                                  second = int(float(tokens[7])),
                                  )
        delta_tsv = float(tokens[9])
        SAT.append(sat)
        EPOCH.append(epoch)
        DELTA_TSV.append(delta_tsv)
    df = pd.DataFrame({
        'Epoch'    : EPOCH,
        'SatName'  : SAT,
        'DeltaTSV' : DELTA_TSV,
    })
    df = df[df['Epoch'] < (df['Epoch'].values[0] + pd.Timedelta(1, unit='day'))]
    df = df.reset_index(drop=True)
    return df

def read_sp3_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()        
    for index, line in enumerate(lines):
        if line.startswith('%c '):
            time_system = line.split()[3]
            assert((time_system == 'GPS') or (time_system == 'ccc'))
            continue
        if line.startswith('* '):
            start_index = index
            break
    lines = lines[start_index:]

    data = []
    for line in lines:
        if line.startswith('* '):
            tokens = line.rstrip().split()
            epoch = datetime.datetime(
                year   = int(tokens[1]),
                month  = int(tokens[2]),
                day    = int(tokens[3]),
                hour   = int(tokens[4]),
                minute = int(tokens[5]),
                second = int(float(tokens[6])),
            )
        elif line.startswith('P'):
            tokens = line.rstrip().split()
            sat = tokens[0][1:]
            x, y, z, delta_t = [float(s) for s in tokens[1:5]]
            x = x * 1e+3
            y = y * 1e+3
            z = z * 1e+3
            delta_t = delta_t * 1e-6
            data.append([epoch, sat, x, y, z, delta_t])
    columns = ['Epoch', 'SatName', 'X', 'Y', 'Z', 'DeltaTSV_SP3']
    df = pd.DataFrame(data, columns=columns)
    df = df[df['Epoch'] < (df['Epoch'].values[0] + pd.Timedelta(1, unit='day'))]
    err_df = df[(df['X'] == 0) & (df['Y'] == 0) & (df['Z'] == 0)]
    for sat in np.unique(err_df['SatName'].values):
        # print(f'read_sp3_file: drop {sat}')
        df = df[df['SatName'] != sat]
    df = df.reset_index(drop=True)
    return df

def read_SINEX_TRO_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    for index, line in enumerate(lines):
        if '+TROP/SOLUTION' in line:
            start_index = index + 2
            break
    lines = lines[start_index:]
    data = []
    for line in lines:
        if '-TROP/SOLUTION' in line:
            break
        tokens  = line.strip().split()
        y, d, s = [int(x) for x in tokens[1].split(':')]
        epoch = datetime.datetime(y+2000, 1, 1) + datetime.timedelta(days=d-1) + datetime.timedelta(seconds=s)
        data.append([epoch] + [1e-3 * float(x) for x in tokens[2:]])
    columns = ['Epoch',
               'TROTOT', 'TROTOT_STD',
               'TGNTOT', 'TGNTOT_STD',
               'TGETOT', 'TGETOT_STD']
    df = pd.DataFrame(data, columns=columns)
    df = df[df['Epoch'] < (df['Epoch'].values[0] + pd.Timedelta(1, unit='day'))]
    df = df.reset_index(drop=True)
    return df

def read_IONEX_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    kw = dict()
    #==============================
    # read header
    #==============================
    for index, line in enumerate(lines):
        tokens = line.strip().split()
        if 'EPOCH OF FIRST MAP' in line:
            kw['time_1'] = np.datetime64(datetime.datetime(
                year   = int(tokens[0]),
                month  = int(tokens[1]),
                day    = int(tokens[2]),
                hour   = int(tokens[3]),
                minute = int(tokens[4]),
                second = int(tokens[5]),
            ))
            continue
        if 'EPOCH OF LAST MAP' in line:
            kw['time_2'] = np.datetime64(datetime.datetime(
                year   = int(tokens[0]),
                month  = int(tokens[1]),
                day    = int(tokens[2]),
                hour   = int(tokens[3]),
                minute = int(tokens[4]),
                second = int(tokens[5]),                
            ))
            continue
        if 'INTERVAL' in line:
            kw['time_delta'] = np.timedelta64(datetime.timedelta(
                seconds=int(tokens[0]),
            ))
            continue
        if 'HGT1 / HGT2 / DHGT' in line:
            h1, h2, dh = [float(x) for x in tokens[0:3]]
            assert(h1 == h2)
            assert(dh == 0.0)
            kw['iono_height'] = h1 * 1000
            continue
        if 'LAT1 / LAT2 / DLAT' in line:
            lat_1, lat_2, lat_delta = [float(x) for x in tokens[0:3]]
            assert((lat_2 - lat_1) * lat_delta > 0)
            if (lat_1 > lat_2):
                flip_lat = True
                lat_1, lat_2 = lat_2, lat_1
                lat_delta = - lat_delta
            else:
                flip_lat = False
            kw['lat_1']     = np.deg2rad(lat_1)
            kw['lat_2']     = np.deg2rad(lat_2)
            kw['lat_delta'] = np.deg2rad(lat_delta)
            continue
        if 'LON1 / LON2 / DLON' in line:
            lng_1, lng_2, lng_delta = [float(x) for x in tokens[0:3]]
            assert((lng_2 - lng_1) * lng_delta > 0)
            if (lng_1 > lng_2):
                flip_lng = True
                lng_1, lng_2 = lng_2, lng_1
                lng_delta = - lng_delta
            else:
                flip_lng = False
            kw['lng_1']     = np.deg2rad(lng_1)
            kw['lng_2']     = np.deg2rad(lng_2)
            kw['lng_delta'] = np.deg2rad(lng_delta)
            continue
        if 'MAPPING FUNCTION' in line:
            assert(tokens[0] == 'COSZ')
            continue
        if 'BASE RADIUS' in line:
            kw['base_radius'] = 1000 * float(tokens[0])
            continue
        if 'EXPONENT' in line:
            TEC_coeff = 10**float(tokens[0])
            continue
        if 'MAP DIMENSION' in line:
            assert(int(tokens[0]) == 2)
            continue
        if 'END OF HEADER' in line:
            line_count = index + 1
            break
    #==============================
    # read data
    #==============================
    N_lat  = 1 + int((kw['lat_2'] - kw['lat_1']) / kw['lat_delta'])
    N_lng  = 1 + int((kw['lng_2'] - kw['lng_1']) / kw['lng_delta'])
    N_time = 1 + (kw['time_2'] - kw['time_1']) // kw['time_delta']
    iono_map = np.zeros((N_time, N_lat, N_lng), dtype=np.float64)

    data_per_line = 16
    lines_per_data = (N_lng + data_per_line - 1) // data_per_line
    
    for time_count in range(N_time):
        assert('START OF TEC MAP' in lines[line_count])
        assert(int(lines[line_count].strip().split()[0]) == time_count + 1)
        line_count += 1
        assert('EPOCH OF CURRENT MAP' in lines[line_count])
        line_count += 1
        for lat_count in range(N_lat):
            assert('LAT/LON1/LON2/DLON/H' in lines[line_count])
            line_count += 1
            values = []
            for i in range(lines_per_data):
                values.extend([int(x) for x in lines[line_count+i].strip().split()])
            if 9999 in values:
                print('Warning: There is non-available TEC values.')
            iono_map[time_count, lat_count, :] = np.array(values).astype(float)
            line_count += lines_per_data
        assert('END OF TEC MAP' in lines[line_count])
        assert(int(lines[line_count].strip().split()[0]) == time_count + 1)
        line_count += 1

    if flip_lat:
        iono_map = np.flip(iono_map, axis=1)
    if flip_lng:
        iono_map = np.flip(iono_map, axis=2)
    iono_map = iono_map * TEC_coeff
    kw['iono_map']  = iono_map
    kw['lat_range'] = np.linspace(kw['lat_1'], kw['lat_2'], N_lat)
    kw['lng_range'] = np.linspace(kw['lng_1'], kw['lng_2'], N_lng)
    return IONEX(**kw)