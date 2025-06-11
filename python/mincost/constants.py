import json
import datetime
from collections import defaultdict
import numpy as np

GPS_ORIGIN_DAY       = datetime.date(1980, 1, 6)
GPS_ORIGIN_DATETIME  = datetime.datetime(1980, 1, 6)
GLONASS_LEAP_SECONDS = 18
BEIDOU_LEAP_SECONDS  = 14
TZ_MSK = datetime.timezone(datetime.timedelta(hours=+3), 'MSK')

WGS84_SEMI_MAJOR_AXIS = 6378137.0
WGS84_SEMI_MINOR_AXIS = 6356752.314245
WGS84_SQUARED_FIRST_ECCENTRICITY  = 6.69437999013e-3
WGS84_SQUARED_SECOND_ECCENTRICITY = 6.73949674226e-3
WGS84_FIRST_ECCENTRICITY  = np.sqrt(WGS84_SQUARED_FIRST_ECCENTRICITY)
WGS84_SECOND_ECCENTRICITY = np.sqrt(WGS84_SQUARED_SECOND_ECCENTRICITY)

LIGHT_SPEED = 299792458.0

OMEGA_EARTH = 7.2921151467e-5
MU_EARTH    = 3.986005e+14

FREQ_GPS_L1  = 1.575420e+09
FREQ_GPS_L5  = 1.176450e+09
FREQ_GAL_E1  = FREQ_GPS_L1
FREQ_GAL_E5A = FREQ_GPS_L5
FREQ_QZS_J1  = FREQ_GPS_L1
FREQ_QZS_J5  = FREQ_GPS_L5
FREQ_BDS_B1I = 1.561098e+09
FREQ_GLO_G1_NOMINAL = 1602.00 * 1e+6
FREQ_GLO_G1_DELTA   = 562.5 * 1e+3

CONSTELLATION_TYPE_MAP = {
    'GPS'     : 1,
    'GLONASS' : 3,
    'QZSS'    : 4,
    'BEIDOU'  : 5,
    'GALILEO' : 6,
}

RAW_STATE_BIT_MAP = {
     0: "Code Lock",
     1: "Bit Sync",
     2: "Subframe Sync",
     3: "Time Of Week Decoded State",
     4: "Millisecond Ambiguity",
     5: "Symbol Sync",
     6: "GLONASS String Sync",
     7: "GLONASS Time Of Day Decoded",
     8: "BEIDOU D2 Bit Sync",
     9: "BEIDOU D2 Subframe Sync",
    10: "Galileo E1BC Code Lock",
    11: "Galileo E1C 2^nd^ Code Lock",
    12: "Galileo E1B Page Sync",
    13: "SBAS Sync",
    14: "Time Of Week Known",
    15: "GLONASS Time Of Day Known",
}
RAW_STATE_BIT_INV_MAP = { value : key for key, value in RAW_STATE_BIT_MAP.items() }

SYSTEM_NAME_MAP = {
    'GPS'     : 'G',
    'GLONASS' : 'R',
    'GALILEO' : 'E',
    'BEIDOU'  : 'C',
    'QZSS'    : 'J',
}

GLONASS_FREQ_CHANNEL_MAP = {
    1 : 1,
    2 : -4,
    3 : 5,
    4 : 6,
    5 : 1,
    6 : -4,
    7 : 5,
    8 : 6,
    9 : -2,
    10 : -7,
    11 : 0,
    12 : -1,
    13 : -2,
    14 : -7,
    15 : 0,
    16 : -1,
    17 : 4,
    18 : -3,
    19 : 3,
    20 : 2,
    21 : 4,
    22 : -3,
    23 : 3,
    24 : 2,
}

QZSS_PRN_SVID_MAP = {
    193 : 1,
    194 : 2,
    199 : 3,
    195 : 4,
}

INIT_B = np.deg2rad(  37.5)
INIT_L = np.deg2rad(-122.2)
INIT_H = 0.0

FREQ_TOL = 100.0
Cn0DbHz_THRESHOLD = 20.0
ReceivedSvTimeUncertaintyNanos_THRESHOLD = 100
RAW_PSEUDO_RANGE_THRESHOLD = 50_000 * 1e+3

CLOCK_TIME_MARGIN = datetime.timedelta(seconds=90)
ORBIT_TIME_MARGIN = datetime.timedelta(hours=3)
IONO_TIME_MARGIN  = datetime.timedelta(hours=2)

EPSILON_M = 0.01
ELEVATION_CUTOFF = np.deg2rad(7.0)
DEFAULT_TROPO_DELAY_M = 2.48

HAVERSINE_RADIUS = 6_371_000

MAGNETIC_DECLINATION = np.deg2rad(10.0)