import numpy as np
import scipy.signal
from scipy.special import sinc

def make_sinc_filter(F_cutoff, dt, n_sinc=3, n_gauss=2):
    N_FLT = round(n_sinc / (2 * F_cutoff * dt))
    x = (n_sinc / (N_FLT + 1)) * (np.arange(2*N_FLT+1) - N_FLT)
    S = sinc(x)
    W = scipy.signal.windows.gaussian(2*N_FLT+1, std=N_FLT/n_gauss)
    flt = S * W
    flt = flt * (1 / np.sum(flt))
    return flt

def main():
    import scipy.signal
    import matplotlib.pyplot as plt

    dt = 2.5 * 1e-3
    F_cutoff = 2.0
    FLT = make_sinc_filter(F_cutoff=F_cutoff, dt=dt)

    f_nyquist = (1 / (2 * dt))
    print(f'dt = {1000*dt} [ms]')
    print(f'f_nyquist = {f_nyquist:.2f} [Hz]')
    
    num = FLT
    den = np.zeros_like(num)
    den[0] = 1
    sys = scipy.signal.TransferFunction(num, den, dt=dt)

    frange = np.logspace(-1, 2, 1000)
    wrange = (2 * np.pi * dt) * frange
    _, mag, phase = sys.bode(wrange)

    assert(len(num) % 2 == 1)
    N_FLT = (len(num) - 1) // 2
    it = np.linspace(-N_FLT, N_FLT, 2*N_FLT+1)
    print(f'N_FLT = {N_FLT}')

    plt.figure()
    plt.plot(it, num)

    plt.figure()
    plt.semilogx(frange, mag)
    plt.grid(True)
    plt.xlim([0.1, 100])
    plt.ylim([-80, 20])

    plt.show()
    return