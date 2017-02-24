#Xk=∑n=0N−1(xn⋅e^(−i 2π k n / N))

#Else use...
#numpy.fft
#scipy.fftpack

import numpy as np
import timeit
import matplotlib.pyplot as plt

#FFT Implimentation

#matrix Multiplication
def DFT_slow(x):
    """Compute the discrete Fourier Transform of the 1D array x"""
    x = np.asarray(x, dtype=float)
    N = x.shape[0]

    n = np.arange(N)
    k = n.reshape((N, 1))
    M = np.exp(-2j * np.pi * k * n / N)
    return np.dot(M, x)


x = np.random.random(1024)
print(np.allclose(DFT_slow(x), np.fft.fft(x)))

def FFT(x):
    """A recursive implementation of the 1D Cooley-Tukey FFT"""
    x = np.asarray(x, dtype=float)
    N = x.shape[0]

    #Split evens and odds
    if N % 2 > 0:
        raise ValueError("size of x must be a power of 2")
    elif N <= 32:  # this cutoff should be optimized
        return DFT_slow(x)
    else:
        X_even = FFT(x[::2])
        X_odd = FFT(x[1::2])
        factor = np.exp(-2j * np.pi * np.arange(N) / N)
        return np.concatenate([X_even + factor[:N / 2] * X_odd,
                               X_even + factor[N / 2:] * X_odd])


x = np.random.random(1024 * 16)
print(np.fft.fft(x))

plt.plot(x)
#plt.plot(np.fft.fft(x))

plt.show()
