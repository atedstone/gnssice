import math
import numpy as np

def gaussfiltcoef(SR, fco):
    """
    Return coefficients of Gaussian low-pass filter.

    SR=sampling rate, fco=cutoff (-3dB) freq, both in Hz. 
    Coeffs for FIR filter of length L (L always odd) are computed.
    This symmetric FIR filter of length L=2N+1 has delay N/SR seconds.
    Examples of use
        Compute Gaussian filter frequency response for SR=1000, fco=50 Hz:
        freqz(gaussfiltcoef(1000,50),1,256,1000);
        Filter signal X sampled at 5kHz with Gaussian filter with fco=500:
        y=filter(gaussfiltcoef(5000,500),1,X);
    SR, fco are not sanity-checked.  WCR 2006-10-11.
    https://www.mathworks.com/matlabcentral/fileexchange/12606-1d-gaussian-lowpass-filter

    Ported to Python by AJT, 2022-07-28. Works with scipy.signal.filtfilt.

    """

    a = 3.011 * fco
    N = math.ceil(0.398 * SR / fco)

    L = 2 * N + 1

    b = []
    for n in range(-N, N):
        b.append(3.011 * (fco/SR) * np.exp(-math.pi * (a*n/SR)**2))
    b = np.array(b)
    b = b / np.sum(b)

    return b