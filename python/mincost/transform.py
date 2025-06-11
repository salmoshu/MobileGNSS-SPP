import numpy as np
from dataclasses import dataclass

import constants as C

@dataclass
class ECEF:
    x: np.array
    y: np.array
    z: np.array

    def to_numpy(self):
        return np.stack([self.x, self.y, self.z], axis=0)

    @staticmethod
    def from_numpy(pos):
        x, y, z = [np.squeeze(w) for w in np.split(pos, 3, axis=-1)]
        return ECEF(x=x, y=y, z=z)

@dataclass
class BLH:
    lat : np.array
    lng : np.array
    hgt : np.array

@dataclass
class ENU:
    east  : np.array
    north : np.array
    up    : np.array

@dataclass
class AZEL:
    elevation : np.array
    azimuth   : np.array
    zenith    : np.array

def BLH_to_ECEF(blh):
    a  = C.WGS84_SEMI_MAJOR_AXIS
    e2 = C.WGS84_SQUARED_FIRST_ECCENTRICITY
    sin_B = np.sin(blh.lat)
    cos_B = np.cos(blh.lat)
    sin_L = np.sin(blh.lng)
    cos_L = np.cos(blh.lng)
    n = a / np.sqrt(1 - e2*sin_B**2)
    x = (n + blh.hgt) * cos_B * cos_L
    y = (n + blh.hgt) * cos_B * sin_L
    z = ((1 - e2) * n + blh.hgt) * sin_B
    return ECEF(x=x, y=y, z=z)

def ECEF_to_BLH_approximate(ecef):
    a = C.WGS84_SEMI_MAJOR_AXIS
    b = C.WGS84_SEMI_MINOR_AXIS
    e2  = C.WGS84_SQUARED_FIRST_ECCENTRICITY
    e2_ = C.WGS84_SQUARED_SECOND_ECCENTRICITY
    x = ecef.x
    y = ecef.y
    z = ecef.z
    r = np.sqrt(x**2 + y**2)
    t = np.arctan2(z * (a/b), r)
    B = np.arctan2(z + (e2_*b)*np.sin(t)**3, r - (e2*a)*np.cos(t)**3)
    L = np.arctan2(y, x)
    n = a / np.sqrt(1 - e2*np.sin(B)**2)
    H = (r / np.cos(B)) - n
    return BLH(lat=B, lng=L, hgt=H)

ECEF_to_BLH = ECEF_to_BLH_approximate

def ECEF_to_ENU(pos, base):
    dx = pos.x - base.x
    dy = pos.y - base.y
    dz = pos.z - base.z
    base_blh = ECEF_to_BLH(base)
    sin_B = np.sin(base_blh.lat)
    cos_B = np.cos(base_blh.lat)
    sin_L = np.sin(base_blh.lng)
    cos_L = np.cos(base_blh.lng)
    e = -sin_L*dx + cos_L*dy
    n = -sin_B*cos_L*dx - sin_B*sin_L*dy + cos_B*dz
    u =  cos_B*cos_L*dx + cos_B*sin_L*dy + sin_B*dz
    return ENU(east=e, north=n, up=u)

def ENU_to_AZEL(enu):
    e = enu.east
    n = enu.north
    u = enu.up
    elevation = np.arctan2(u, np.sqrt(e**2 + n**2))
    azimuth   = np.arctan2(e, n)
    zenith    = (0.5 * np.pi) - elevation
    return AZEL(elevation=elevation,
                azimuth=azimuth,
                zenith=zenith)

def ECEF_to_AZEL(pos, base):
    return ENU_to_AZEL(ECEF_to_ENU(pos, base))

def haversine_distance(blh_1, blh_2):
    dlat = blh_2.lat - blh_1.lat
    dlng = blh_2.lng - blh_1.lng
    a = np.sin(dlat/2)**2 + np.cos(blh_1.lat) * np.cos(blh_2.lat) * np.sin(dlng/2)**2
    dist = 2 * C.HAVERSINE_RADIUS * np.arcsin(np.sqrt(a))
    return dist

def hubenys_distance(blh_1, blh_2):
    Rx = C.WGS84_SEMI_MAJOR_AXIS
    Ry = C.WGS84_SEMI_MINOR_AXIS
    E2 = C.WGS84_SQUARED_FIRST_ECCENTRICITY
    num_M = Rx * (1 - E2)
    Dy = blh_1.lat - blh_2.lat
    Dx = blh_1.lng - blh_2.lng
    P  = 0.5 * (blh_1.lat + blh_2.lat)
    W  = np.sqrt(1 - E2 * np.sin(P)**2)
    M  = num_M / W**3
    N  = Rx / W
    d2 = (Dy * M)**2 + (Dx * N * np.cos(P))**2
    d  = np.sqrt(d2)
    return d

def jacobian_BLH_to_ECEF(blh):
    a  = C.WGS84_SEMI_MAJOR_AXIS
    e2 = C.WGS84_SQUARED_FIRST_ECCENTRICITY
    B = blh.lat
    L = blh.lng
    H = blh.hgt
    cos_B = np.cos(B)
    sin_B = np.sin(B)
    cos_L = np.cos(L)
    sin_L = np.sin(L)
    N = a / np.sqrt(1 - e2*sin_B**2)
    dNdB = a * e2 * sin_B * cos_B * (1 - e2*sin_B**2)**(-3/2)
    N_plus_H = N + H
    cos_B_cos_L = cos_B * cos_L
    cos_B_sin_L = cos_B * sin_L
    sin_B_cos_L = sin_B * cos_L
    sin_B_sin_L = sin_B * sin_L

    dXdB = dNdB*cos_B_cos_L - N_plus_H*sin_B_cos_L
    dYdB = dNdB*cos_B_sin_L - N_plus_H*sin_B_sin_L
    dZdB = (1-e2)*dNdB*sin_B + (1-e2)*N_plus_H*cos_B

    dXdL = - N_plus_H * cos_B_sin_L
    dYdL =   N_plus_H * cos_B_cos_L
    dZdL = np.zeros_like(dXdL)

    dXdH = cos_B_cos_L
    dYdH = cos_B_sin_L
    dZdH = sin_B

    J = np.stack([[dXdB, dXdL, dXdH],
                  [dYdB, dYdL, dYdH],
                  [dZdB, dZdL, dZdH]], axis=0)
    axes = list(range(2, J.ndim)) + [0, 1]
    J = np.transpose(J, axes)
    return J

def jacobian_ECEF_to_ENU(blh):
    B = blh.lat
    L = blh.lng
    cos_B = np.cos(B)
    sin_B = np.sin(B)
    cos_L = np.cos(L)
    sin_L = np.sin(L)
    
    dEdX = -sin_L
    dEdY =  cos_L
    dEdZ = np.zeros_like(dEdX)
    
    dNdX = -sin_B*cos_L
    dNdY = -sin_B*sin_L
    dNdZ =  cos_B

    dUdX = cos_B*cos_L
    dUdY = cos_B*sin_L
    dUdZ = sin_B

    J = np.stack([[dEdX, dEdY, dEdZ],
                  [dNdX, dNdY, dNdZ],
                  [dUdX, dUdY, dUdZ]], axis=0)
    axes = list(range(2, J.ndim)) + [0, 1]
    J = np.transpose(J, axes)
    return J

def jacobian_BL_to_EN(BLH):
    J_ECEF_BLH = jacobian_BLH_to_ECEF(BLH)
    J_ENU_ECEF = jacobian_ECEF_to_ENU(BLH)
    J_EN_BL    = np.einsum('nij,njk->nik', J_ENU_ECEF[:, 0:2, :], J_ECEF_BLH[:, :, 0:2])
    return J_EN_BL

def pd_haversine_distance(df1, df2):
    blh1 = BLH(
        lat=np.deg2rad(df1['latDeg'].values),
        lng=np.deg2rad(df1['lngDeg'].values),
        hgt=0,
    )
    blh2 = BLH(
        lat=np.deg2rad(df2['latDeg'].values),
        lng=np.deg2rad(df2['lngDeg'].values),
        hgt=0,
    )
    return haversine_distance(blh1, blh2)