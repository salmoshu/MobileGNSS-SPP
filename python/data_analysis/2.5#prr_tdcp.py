import os
import sys
import shutil
import numpy as np
import matplotlib.pyplot as plt

# 对应文档的2.5节（TDCP）

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.join(current_dir, '..')
src_dir = os.path.join(current_dir, '../rtklibpy')
sys.path.append(parent_dir)
sys.path.append(src_dir)
cfgfile = os.path.join(src_dir, 'config_spp.py')
shutil.copyfile(cfgfile, '__ppk_config.py')

import __ppk_config as cfg
import rinex as rn
from rtkcmn import rCST, sat2prn, id2sat, gpst2utc, pos2ecef
from rtkpos import rtkinit

'''
    1. 通过rtklib-py获取观测数据(仅分析L1)
    2. 分析观测数据的特性：
        - Q1. 载波相位的有效性
            TDCP与通过多普勒计算的伪距率之间的差异应为零(零上下波动)。
            比较TDCP与通过多普勒计算的伪距率之间的差异，是否存在系统误差。
        - Q2. 伪距率是否由多普勒获取
            有的设备伪距率由多普勒获取，而有的设备则由载波相位获取；
            具体通过TDCP与通过多普勒计算的伪距率进行比较，如果两者完全相同说明是由TDCP获取的，反之则是由多普勒获取的。
        - Q3. 伪距和多普勒是否相关
            通常伪距计算的接收机钟差(单位时间变化率)和从伪距率计算的接收机时钟漂移应该一致
        - Q4. 伪距和载波相位是否相关
            TDCP中包含的接收机钟差变化也与从伪距计算的接收机钟差变化不一致；
            TDCP计算的钟差变化：历元间载波相位观测相减得到的钟差差值
            伪距计算的钟差变化：历元间伪距观测相减得到的钟差差值
'''

def time2unix(t):
    t = gpst2utc(t)
    return t.time + t.sec

def get_raw_obs(obs_path):
    rov = rn.rnx_decode(cfg)
    nav = rtkinit(cfg)

    rov.decode_obsfile(nav, obs_path, None)
    return rov.obslist

def get_sat_info(nav_path):
    # load rover obs
    rov = rn.rnx_decode(cfg)
    nav = rtkinit(cfg)
    rov.decode_obsfile(nav, obs_path, None)

    return rov.decode_nav(nav_path, nav)

# gt for staitc postion
def get_middle_obs(obss):
    nepoch = len(obss)

    last_obs = obss[0]
    for i in range(0, nepoch):
        nsat = len(obss[i].sat)
        obss[i].prr0 = [np.nan] * nsat
        obss[i].tdcp0 = [np.nan] * nsat
        obss[i].dpr0 = [np.nan] * nsat

        for j in range(nsat):
            sat = obss[i].sat[j]
            sys, _ = sat2prn(sat)
            freq_idx = cfg.freq_ix0[sys]
            freq = cfg.freq[freq_idx]

            # get prr
            prr0 = - obss[i].D[j, 0] * rCST.CLIGHT / freq
            obss[i].prr0[j] = prr0

        # get tdcp
        if i == 0:
            continue
        last_nsat = len(last_obs.sat)
        for j in range(0, nsat):
            for k in range(0, last_nsat):
                if obss[i].sat[j] == last_obs.sat[k]:
                    # tdcp = dADR / dt, dADR = dL * c / f
                    dADR0 = (obss[i].L[j, 0] - last_obs.L[k, 0]) * rCST.CLIGHT / freq
                    tdcp0 = dADR0 / (time2unix(obss[i].t) - time2unix(last_obs.t))

                    dpr0 = obss[i].P[j, 0] - last_obs.P[k, 0]
                    dpr0 = dpr0 / (time2unix(obss[i].t) - time2unix(last_obs.t))

                    obss[i].tdcp0[j] = tdcp0
                    obss[i].dpr0[j] = dpr0

                    break
            else:
                obss[i].tdcp0[j] = np.nan
        last_obs = obss[i]

    return obss

# Q1. 载波相位的有效性 & TDCP中是否存在系统误差
def get_res_Q1(obss, id, mode=0):
    sat = id2sat(id)
    sys, prn = sat2prn(sat)
    nepoch = len(obss)

    prr_tdcp_idx = []
    prr_tdcp_val = []
    slip_idx = []
    slip_val = []

    for i in range(0, nepoch):
        nsat = len(obss[i].sat)
        for j in range(0, nsat):
            if obss[i].sat[j] == sat:
                prr0 = obss[i].prr0[j]
                tdcp0 = obss[i].tdcp0[j]
                dpr0 = obss[i].dpr0[j]

                freq_idx = cfg.freq_ix0[sys]
                freq = cfg.freq[freq_idx]

                if prr0 is not None and tdcp0 is not None:
                    prr_tdcp_idx.append(i)

                    if mode == 0:
                        ppr_tdcp_val = dpr0 - tdcp0
                    elif mode == 1:
                        ppr_tdcp_val = dpr0 - prr0
                    elif mode == 2:
                        ppr_tdcp_val = tdcp0 - prr0
                        # ppr_tdcp_val = tdcp0

                    prr_tdcp_val.append(ppr_tdcp_val)
                    if obss[i].lli[j, 0] & 1 != 0:
                        slip_idx.append(i)
                        slip_val.append(ppr_tdcp_val)

    return np.array(prr_tdcp_idx), np.array(prr_tdcp_val), np.array(slip_idx), np.array(slip_val)

if __name__ == "__main__":
    MODE_P_L = 0
    MODE_P_D = 1
    MODE_L_D = 2

    mode = MODE_L_D
    sta_id = 'G02'
    obs_path = r'.\rover.obs'
    obss = get_raw_obs(obs_path)

    # nav_path = r'.\rover.nav'
    # nav = get_sat_info(nav_path)

    obss = get_middle_obs(obss)

    prr_tdcp_idx, prr_tdcp_val, slip_idx, slip_val= get_res_Q1(obss, sta_id, mode=mode)
    prr_tdcp_idx = prr_tdcp_idx[1:]
    prr_tdcp_val = prr_tdcp_val[1:]
    prr_tdcp_val_cpy = prr_tdcp_val.copy()

    plt_str = ''
    if mode == MODE_P_L:
        plt_str = 'P-L'
    elif mode == MODE_P_D:
        plt_str = 'P-D'
    elif mode == MODE_L_D:
        plt_str = 'L-D'

    prr_tdcp_val[np.abs(prr_tdcp_val) > 10] = np.nan
    prr_tdcp_mean = np.mean(prr_tdcp_val)
    prr_tdcp_val[np.abs(prr_tdcp_val - prr_tdcp_mean) > 10] = np.nan
    prr_tdcp_mean = np.nanmean(prr_tdcp_val)
    print(f'{plt_str}(mean): {prr_tdcp_mean:.3f}')

    plt.plot(prr_tdcp_idx, prr_tdcp_val_cpy, label=sta_id+f': {plt_str}:mean='+f'{prr_tdcp_mean:.3f}'+'m/s')
    # plt.plot(prr_tdcp_idx, prr_tdcp_val_cpy, label=sta_id+f': tdcp')
    plt.plot(slip_idx, slip_val, 'ro', label='cycle slip')
    plt.legend()
    plt.grid()
    plt.show()
    