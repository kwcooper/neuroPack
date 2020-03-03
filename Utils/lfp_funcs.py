# More that just me (Credits):
# Voyteck lab 
# McNaughton lab

import numpy as np 
import scipy as sp
import scipy.io
import scipy.signal

# 1. Filters and utils
# 2. PAC analysis
# 3. Ripple utils



###########################################################################################
##  Filters and LFP utils
###########################################################################################

# Define butter bandpass filter
def butter_bandpass(lowcut, highcut, fs, order=4):
    #lowcut is the lower bound of the frequency that we want to isolate
    #hicut is the upper bound of the frequency that we want to isolate
    #fs is the sampling rate of our data
    nyq = 0.5 * fs #nyquist frequency - see http://www.dspguide.com/ if you want more info
    low = float(lowcut) / nyq
    high = float(highcut) / nyq
    b, a = sp.signal.butter(order, [low, high], btype='band')
    return b, a

# Use to call the OG goto filter (IMO)
def butter_bandpass_filter(mydata, lowcut, highcut, fs, order=4):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = sp.signal.filtfilt(b, a, mydata)
    return y

# Compute the RMS of the signal
def RMS(lfp):
	return np.sqrt(np.mean(lfp**2))

# https://stackoverflow.com/a/8260297/5311230
def RMS_window(a, window_size):
  a2 = np.power(a,2)
  window = np.ones(window_size)/float(window_size)
  return np.sqrt(np.convolve(a2, window, 'valid'))


###########################################################################################
##  Phase Coupling Analysis
###########################################################################################

# Used for PAC analysis below
# can't do a simple pearson correlation here because the phase of theta 
# is periodic (between -pi and pi), and you'll get incorrect 
# results if you run a pearson on periodic data...
# This function should correct for that. 
def circCorr(ang,line):
    n = len(ang)
    rxs = sp.stats.pearsonr(line,np.sin(ang))
    rxs = rxs[0]
    rxc = sp.stats.pearsonr(line,np.cos(ang))
    rxc = rxc[0]
    rcs = sp.stats.pearsonr(np.sin(ang),np.cos(ang))
    rcs = rcs[0]
    rho = np.sqrt((rxc**2 + rxs**2 - 2*rxc*rxs*rcs)/(1-rcs**2)) #r
    r_2 = rho**2 #r squared
    pval = 1- sp.stats.chi2.cdf(n*(rho**2),1)
    standard_error = np.sqrt((1-r_2)/(n-2))

    return rho, pval, r_2, standard_error


# Phase amplitude coupling analysis
# give it the phase providing band
# the amplitude providing band
# your data and the fs
def PAC(lfp, fs, ppb, apb, pp=1, fig=0): 
	# ppb = [4,8] for theta; apb = [80-125] for gamma
	#calculating phase of theta
	phase_data = butter_bandpass_filter(lfp, ppb[0], ppb[1], round(fs));
	phase_data = scipy.signal.hilbert(phase_data);
	phase_data = np.angle(phase_data);

	#calculating amplitude envelope of high gamma
	amp_data = butter_bandpass_filter(lfp, apb[0], apb[1], round(fs));
	amp_data = scipy.signal.hilbert(amp_data);
	amp_data = np.abs(amp_data);

	# Now run stats to see what we get:
	# todo: port this to a variable
	rho, pval, r_2, standard_error = circCorr(phase_data, amp_data)

	if fig: # TODO rename and split to seperate func
		#let's look at a small chunk of our data
		plt.figure(figsize = (15,6));
		dat_norm = (dat[1:int(fs)*2]-np.mean(dat[1:int(fs)*2]))/np.std(dat[1:int(fs)*2])
		#plt.plot(dat_norm,label= 'Raw Data'); 
		plt.plot(phase_data[1:int(fs)*2]*np.mean(amp_data),label= 'Phase of Theta');
		plt.plot(amp_data[1:int(fs)*2],label= 'Amplitude of High Gamma'); 
		plt.xlabel('Two Seconds of Theta Phase and High Gamma Amplitude')
		plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		plt.show()

	if pp:
		#print(rho, pval, r_2, standard_error)
		print(rho, pval, r_2, standard_error)
	else:
		return rho, pval, r_2, standard_error







