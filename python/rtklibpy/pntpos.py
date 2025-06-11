"""
pntpos.py : module for standalone positioning

Copyright (c) 2021 Rui Hirokawa (from CSSRLIB)
Copyright (c) 2022 Tim Everett
"""
import numpy as np
from numpy.linalg import norm, lstsq
from rtkcmn import rCST, ecef2pos, geodist, satazel, ionmodel, tropmodel, \
     Sol, tropmapf, uGNSS, trace, timeadd, \
     xyz2enu
import rtkcmn as gn
from ephemeris import seleph, satposs
from rinex import rcvstds

NX =        7           # num of estimated parameters, pos + clock
MAXITR =    10          #  max number of iteration or point pos
ERR_ION =   5.0         #  ionospheric delay Std (m)
ERR_TROP =  3.0         #  tropspheric delay Std (m)
ERR_SAAS =  0.3         #  Saastamoinen model error Std (m)
ERR_BRDCI = 0.5         #  broadcast ionosphere model error factor
ERR_CBIAS = 0.3         #  code bias error Std (m)
REL_HUMI =  0.7         #  relative humidity for Saastamoinen model
MIN_EL = np.deg2rad(5)  #  min elevation for measurement

def varerr(nav, sys, el, snr, rcvstd):
    """ variation of measurement """
    # RTKLIB:
    # var1 = a^2
    # var2 = (b / sin(el))^2
    # var3 = c^2 * 10^(0.1*(snr_max-snr))
    # var4 = (d*rcv_std)^2
    # var = erate * fact * (var1+var2+var3) + var4
    s_el = np.sin(el)
    fact = nav.eratio[0] # only for L1
    if s_el <= 0.0:
        return 0.0
    a = 0.003 # use simple weighting, since only used for approx location
    b = nav.err[2] 
    c = nav.err[4]
    var1 = a**2
    var2 = (b / s_el)**2
    var3 = 0
    if c > 0:
        dsnr = max(nav.snrmax-snr, 0)
        var3 = c**2 * (10**(0.1*dsnr))
    var = var1 + var2 + var3
    var *= (fact**2)
    var *= nav.efact[sys]
    return  var

def gettgd(sat, eph, type=0):
    """ get tgd: 0=E5a, 1=E5b  """
    sys = gn.sat2prn(sat)[0]
    if sys == uGNSS.GLO:
        return eph.dtaun * rCST.CLIGHT
    else:
        return eph.tgd[type] * rCST.CLIGHT
    

def prange(nav, obs, i):
    eph = seleph(nav, obs.t, obs.sat[i])
    P1 = obs.P[i,0]
    if P1 == 0 or P1 < 1e7:
        return 0
    sys = gn.sat2prn(obs.sat[i])[0]
    if sys == uGNSS.GPS or sys == uGNSS.QZS:
        b1 = gettgd(obs.sat[i], eph, 0)
        return P1 - b1
    elif sys == uGNSS.GAL:
        b1 = gettgd(obs.sat[i], eph, 1)
        return P1 - b1
    elif sys == uGNSS.BDS:
        b1 = gettgd(obs.sat[i], eph, 0)
        return P1 - b1
    else:  # GLONASS
        # TODO:  abstract hard coded freqs
        gamma = nav.freq[4] / nav.freq[5]  # G1/G2
        b1 = gettgd(obs.sat[i], eph, 0)
        return P1 - b1 / (gamma - 1)

def rescode(iter, obs, nav, rs, dts, svh, x):
    """ calculate code residuals """
    ns = len(obs.sat)  # measurements
    trace(4, 'rescode : n=%d\n' % ns)
    v = np.zeros(ns + NX - 3)
    H = np.zeros((ns + NX - 3, NX))
    mask = np.zeros(NX - 3) # clk states 
    azv = np.zeros(ns)
    elv = np.zeros(ns)
    var = np.zeros(ns + NX - 3)
    
    rr = x[0:3].copy()
    dtr = x[3]
    pos = ecef2pos(rr)
    trace(3, 'rescode: rr=%.3f %.3f %.3f\n' % (rr[0], rr[1], rr[2]))
    rcvstds(nav, obs) # decode stdevs from receiver
    
    nv = 0
    vion = 0
    vtrp = 0
    for i in np.argsort(obs.sat):
        sys = nav.sysprn[obs.sat[i]][0]
        if norm(rs[i,:]) < rCST.RE_WGS84:
            continue
        if gn.satexclude(obs.sat[i], var[i], svh[i], nav):
            continue
        # geometric distance and elevation mask
        r, e = geodist(rs[i], rr)
        if r < 0:
            continue
        [az, el] = satazel(pos, e)
        if el < nav.elmin:
            continue
        if iter > 0:
            # test CNR
            if obs.S[i,[0]] < nav.cnr_min[0]:
                continue
            # ionospheric correction
            dion = ionmodel(obs.t, pos, az, el, nav.ion)
            freq = gn.sat2freq(obs.sat[i], 0, nav)
            dion *= (nav.freq[0] / freq)**2
            # tropospheric correction
            trop_hs, trop_wet, _ = tropmodel(obs.t, pos, el, REL_HUMI)
            mapfh, mapfw = tropmapf(obs.t, pos, el)
            dtrp = mapfh * trop_hs + mapfw * trop_wet

            # vion = (dion*1)**2 * (nav.freq[0] / freq)**4
            # vtrp = (0.3 / el + 0.1)**2
        else:
            dion = dtrp = 0
        # psendorange with code bias correction
        P = prange(nav, obs, i)
        if P == 0:
            continue
        # pseudorange residual
        v[nv] = P - (r + dtr - rCST.CLIGHT * dts[i][0] + dion + dtrp)
        trace(4, 'sat=%d: v=%.3f P=%.3f r=%.3f dtr=%.6f dts=%.6f dion=%.3f dtrp=%.3f\n' %
              (obs.sat[i],v[nv],P,r,dtr,dts[i][0],dion,dtrp))
        # design matrix 
        H[nv, 0:3] = -e
        H[nv, 3] = 1
        # time system offset and receiver bias correction
        if sys == uGNSS.GLO:
            v[nv] -= x[4]
            H[nv, 4] = 1.0
            mask[1] = 1
        elif sys == uGNSS.GAL:
            v[nv] -= x[5]
            H[nv, 5] = 1.0
            mask[2] = 1
        elif sys == uGNSS.BDS:
            v[nv] -= x[6]
            H[nv, 6] = 1.0
            mask[3] = 1
        else:
            mask[0] = 1
            
        azv[nv] = az
        elv[nv] = el
        snr = obs.S[i,0]
        # vmeas = 0.3**2
        vmeas = 0.0
        vmeas += varerr(nav, sys, el, snr, nav.rcvstd[obs.sat[i]-1,0])
        var[nv] = vmeas + vion + vtrp
        nv += 1

    # constraint to avoid rank-deficient
    for i in range(NX - 3):
        if mask[i] == 0:
            v[nv] = 0.0
            H[nv, i+3] = 1
            var[nv] = 0.01
            nv += 1
    v = v[0:nv]
    H = H[0:nv, :]
    azv = azv[0:nv]
    elv = elv[0:nv]
    var = var[0:nv]
    return v, H, azv, elv, var

def clkOffset(n, H, x, v):
    """接收机时钟偏移校正（系统间）"""
    # 识别各GNSS系统的索引位置
    idx_G = idx_E = idx_C = 0
    for i in range(n):
        if H[i,3] == 1 and H[i,4] == 0 and H[i,5] == 0:  # GPS
            idx_G = i
        elif H[i,4] == 1:  # Galileo
            idx_E = i
        elif H[i,5] == 1:  # BDS
            idx_C = i

    # 计算各系统时钟偏移中值
    offsetG = np.median(v[:idx_G+1]) if idx_G+1 > 2 else 0
    offsetE = np.median(v[idx_G+1:idx_E+1]) - offsetG if idx_E-idx_G > 2 else 0
    offsetC = np.median(v[idx_E+1:idx_C+1]) - offsetG if idx_C-idx_E > 2 else 0

    # 更新状态向量
    x[3] += offsetG  # GPS时钟
    x[4] += offsetE  # Galileo时钟偏移
    x[5] += offsetC  # BDS时钟偏移

    trace(3, f"clk offset: {' '*36} dtr={x[3]:.5f} {x[3]+x[4]:.5f} {x[3]+x[5]:.5f}\n")

    # 校正残差向量
    for i in range(n):
        if H[i,3] == 1:  # GPS
            v[i] -= offsetG
        if H[i,4] == 1:  # Galileo
            v[i] -= offsetE
        if H[i,5] == 1:  # BDS
            v[i] -= offsetC

def estpos(obs, nav, rs, dts, svh, sol):
    """ estimate position and clock errors with standard precision """
    x = np.zeros(NX)
    x[0:3] = nav.x[0:3]
    trace(3, 'estpos  : n=%d\n' % len(rs))
    for iter in range(MAXITR):
        v, H, az, el, var = rescode(iter, obs, nav, rs[:,0:3], dts, svh, x)
        nv = len(v)    

        # # 新增时钟偏移校正
        # if iter == 0:
        #     clkOffset(nv, H, x, v)

        if nv < NX:
            trace(3, 'estpos: lack of valid sats nsat=%d nv=%d\n' % 
                  (len(obs.sat), nv))
            return sol

        # weight by residuals
        v_fabs = np.abs(v)
        v_fabs[v_fabs < 0.001] = 0.001
        v_fabs[v_fabs > 10000.0] = 10000.0
        v /= v_fabs
        H /= v_fabs[:,None]

        # weight by variance (lsq uses sqrt of weight)
        std = np.sqrt(var)
        v /= std
        H /= std[:,None]
        # least square estimation
        dx = lstsq(H, v, rcond=None)[0]
        x += dx
        if norm(dx) < 1e-4:
            break
    else: # exceeded max iterations
        sol.stat = gn.SOLQ_NONE
        trace(3, 'estpos: solution did not converge\n')
    sol.stat = gn.SOLQ_SINGLE
    sol.t = timeadd(obs.t, -x[3] / rCST.CLIGHT )
    sol.dtr = x[3:5] / rCST.CLIGHT
    sol.rr[0:3] = x[0:3]
    sol.rr[3:6] = 0
    Q = np.linalg.inv(H.T @ H)
    sol.qr[0:3, 0:3] = Q[0:3, 0:3]
    sol.ns = nv

def resdop(obs, nav, rs, dts, rr, x, azel, vsat):
    """
    Calculate range rate residuals using Doppler observations.

    Args:
        obs: Observation data (obsd_t-like structure).
        nav: Navigation data (nav_t-like structure).
        rs: Satellite positions and velocities [x, y, z, vx, vy, vz] (n x 6).
        dts: Satellite clock errors [dts, dts_dot] (n x 2).
        rr: Receiver position [x, y, z] (3 x 1).
        x: Current velocity state [vx, vy, vz, dtr_drift] (4 x 1).
        azel: Azimuth and elevation angles [az, el] (n x 2).
        vsat: Satellite validity flags (n x 1).

    Returns:
        v: Range rate residuals (m/s).
        H: Design matrix.
        var: Variance of residuals.
    """
    ns = len(obs.sat)  # Number of satellites
    trace(3, 'resdop  : n=%d\n' % ns)

    # Initialize arrays
    v = np.zeros(ns)
    H = np.zeros((ns, 4))  # State: [vx, vy, vz, dtr_drift]
    var = np.zeros(ns)

    # Convert receiver position to geodetic coordinates
    pos = ecef2pos(rr)
    E = xyz2enu(pos)  # ECEF to ENU transformation matrix

    nv = 0
    for i in range(ns):
        # Check if Doppler observation is valid
        if obs.D[i, 0] == 0.0 or vsat[i] or norm(rs[i, 3:6]) <= 0.0:
            continue

        # Get signal frequency
        freq = gn.sat2freq(obs.sat[i], 0, nav)
        if freq == 0.0:
            continue

        # Line-of-sight (LOS) vector in ENU
        cosel = np.cos(azel[i, 1])
        a = np.array([
            np.sin(azel[i, 0]) * cosel,
            np.cos(azel[i, 0]) * cosel,
            np.sin(azel[i, 1])
        ])
        e = E.T @ a  # Transform to ECEF

        # Satellite velocity relative to receiver
        vs = rs[i, 3:6] - x[0:3]  # vs = v_sat - v_receiver

        # Range rate with Earth rotation correction
        rate = np.dot(vs, e) + rCST.OMGE / rCST.CLIGHT * (
            rs[i, 4] * rr[0] + rs[i, 1] * x[0] -
            rs[i, 3] * rr[1] - rs[i, 0] * x[1]
        )

        # Standard deviation of range rate error (m/s)
        err = nav.err[7] if hasattr(nav, 'err') else 1.0  # Doppler error (Hz)
        sig = err * rCST.CLIGHT / freq if err > 0.0 else 1.0

        # Range rate residual (m/s)
        v[nv] = (-obs.D[i, 0] * rCST.CLIGHT / freq - (rate + x[3] - rCST.CLIGHT * dts[i, 1])) / sig

        # Design matrix
        H[nv, 0:3] = -e / sig  # Velocity components
        H[nv, 3] = 1.0 / sig    # Clock drift

        # Variance (already scaled by sig)
        var[nv] = 1.0  # Since we divided by sig, variance is 1 after scaling

        nv += 1

    # Trim arrays
    v = v[:nv]
    H = H[:nv, :]
    var = var[:nv]

    return v, H, var

def estvel(obs, nav, rs, dts, vsat, sol):
    """
    Estimate receiver velocity using Doppler observations.

    Args:
        obs: Observation data (obsd_t-like structure).
        nav: Navigation data (nav_t-like structure).
        rs: Satellite positions and velocities [x, y, z, vx, vy, vz] (n x 6).
        dts: Satellite clock errors [dts, dts_dot] (n x 2).
        vsat: Satellite validity flags (n x 1).
        sol: Solution object to store results.

    Returns:
        sol: Updated solution object with velocity and covariance.
    """
    # Compute azimuth and elevation angles (needed for resdop)
    rr = sol.rr[0:3]  # Receiver position from position solution
    pos = ecef2pos(rr)
    azel = np.zeros((len(obs.sat), 2))
    for i in range(len(obs.sat)):
        r, e = geodist(rs[i, 0:3], rr)
        if r < 0:
            continue
        az, el = satazel(pos, e)
        azel[i, 0] = az
        azel[i, 1] = el

    # Initialize velocity state [vx, vy, vz, dtr_drift]
    x = np.zeros(4)
    trace(3, 'estvel  : n=%d\n' % len(rs))

    for iter in range(MAXITR):
        # Compute range rate residuals
        v, H, var = resdop(obs, nav, rs, dts, rr, x, azel, vsat)
        nv = len(v)

        # Check if enough observations
        if nv < 4:  # Need at least 4 observations for [vx, vy, vz, dtr_drift]
            trace(3, 'estvel: lack of valid sats nv=%d\n' % nv)
            return sol

        # Weight by variance
        std = np.sqrt(var)
        v /= std
        H /= std[:, None]

        # Least squares estimation
        dx = lstsq(H, v, rcond=None)[0]
        x += dx

        # Check convergence
        if norm(dx) < 1e-6:
            trace(3, 'estvel : vx=%.3f vy=%.3f vz=%.3f, n=%d\n' % (x[0], x[1], x[2], nv))
            break
    else:
        # Exceeded max iterations
        sol.stat = gn.SOLQ_NONE
        trace(3, 'estvel: solution did not converge\n')
        return sol

    # Update solution
    sol.rr[3:6] = x[0:3]  # Store velocity [vx, vy, vz]
    Q = np.linalg.inv(H.T @ H)  # Covariance matrix
    # Store velocity covariance (diagonal and off-diagonal elements)
    sol.qv[0:3, 0:3] = Q[0:3, 0:3]

def pntpos(obs, nav):
    """ single-point positioning ----------------------------------------------------
    * compute receiver position, velocity, clock bias by single-point positioning
    * with pseudorange and doppler observables
    * args   : obs      I   observation data
    *          nav      I   navigation data
    * return : sol      O   """
    sol = Sol()
    rs, _, dts, svh = satposs(obs, nav)
    estpos(obs, nav, rs, dts, svh, sol)
    estvel(obs, nav, rs, dts, svh, sol)
    return sol
    

