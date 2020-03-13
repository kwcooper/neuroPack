# More that just me (Credits):
# Newman lab, Voyteck lab, McNaughton lab, Fortin Lab

import numpy as np 
import scipy as sp # TODO: Clean this up... 
import scipy.io
import scipy.signal

# 1. Filters and utils
# 2. PAC analysis
# 3. Ripple utils


###########################################################################################
##  1. Filters and LFP utils
###########################################################################################

# TODO: Change filter bands to lists inds

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
def butter_bandpass_filter(lfp, lowcut, highcut, fs, order=4):
	b, a = butter_bandpass(lowcut, highcut, fs, order=order)
	y = sp.signal.filtfilt(b, a, lfp)
	return y

# FIR filter like fir1 in MATLAB
# Adapted from:
# https://github.com/voytekresearch/pacpy/blob/0cf3ed0c39656f1c512d1e2dcd9db5b7ba316e3f/pacpy/filt.py#L9
def fir1(lfp, lowcut, highcut, fs, w=3):
	# w = Length of the filter (by cycles)
	nyq = np.float(fs / 2)
	if np.any(np.array([lowcut, highcut]) > nyq):
		raise ValueError('Filter frequencies must be below nyquist rate.')
	Ntaps = int(np.floor(w * fs / lowcut))
	taps = sp.signal.firwin(Ntaps, [lowcut/nyq, highcut/nyq], pass_zero=False)
	return sp.signal.filtfilt(taps, [1], lfp)

# Compute the RMS of the signal
def RMS(lfp):
	return np.sqrt(np.mean(lfp**2))

# inspired by https://stackoverflow.com/a/8260297/5311230
def RMS_window(a, window_size):
  a2 = np.power(a, 2)
  window = np.ones(window_size)/float(window_size)
  return np.sqrt(np.convolve(a2, window, 'valid'))


# Func to pad list of lists (by max in list) to make them the same len:
# There's a cleaner way, I know, was in a rush...
def pad_items(listOlist):
	maxLen = max([len(i) for i in listOlist])
	padList = []
	for ll in listOlist:
		targLen = maxLen - len(ll)
		padl, padr = [0]*int((targLen / 2)), [0]*int((targLen / 2)) #np.zeros(int((targLen / 2)))
		if targLen % 2 == 0:
			padl.extend(ll)
			padl.extend(padr[:-1])
			padList.append(padl)
		else:
			padl.extend(ll)
			padl.extend(padr)
			padList.append(padl)
	return padList


###########################################################################################
##  2. Phase Amplitude Coupling Analysis (Rudimentary)
###########################################################################################

# TODO:
# 	Add separate arguments for low and hi bands for multi-region analysis...

def getPhaseAmp(lfp, fs, ppb, apb): 
	#calculating phase of theta
	phase_data = butter_bandpass_filter(lfp, ppb[0], ppb[1], round(fs));
	phase_data = scipy.signal.hilbert(phase_data);
	phase_data = np.angle(phase_data);

	#calculating amplitude envelope of high gamma
	amp_data = butter_bandpass_filter(lfp, apb[0], apb[1], round(fs));
	amp_data = scipy.signal.hilbert(amp_data);
	amp_data = np.abs(amp_data);

	# lo, hi
	return phase_data, amp_data


# Used for PAC analysis below
# can't do a simple pearson correlation here because the phase of theta 
# is periodic (between -pi and pi), and you'll get incorrect 
# results if you run a pearson on periodic data...
# This function should correct for that. 
def circCorr(ang, line):
	n = len(ang)
	rxs = sp.stats.pearsonr(line, np.sin(ang))
	rxc = sp.stats.pearsonr(line, np.cos(ang))
	rcs = sp.stats.pearsonr(np.sin(ang), np.cos(ang))
	rxs, rxc, rcs = rxs[0], rxc[0], rcs[0]
	rho = np.sqrt((rxc**2 + rxs**2 - 2*rxc*rxs*rcs)/(1-rcs**2)) #r
	r_2 = rho**2 #r squared
	pval = 1.0 - sp.stats.chi2.cdf(float(n)*(float(rho)**2.0), 1.0)
	stand_error = np.sqrt((1-r_2)/(n-2))

	return rho, pval, r_2, stand_error

# https://gist.github.com/kn1cht/89dc4f877a90ab3de4ddef84ad91124e
def corrcoef(x, y, deg=True, test=False):
	'''Circular correlation coefficient of two angle data(default to degree)
	Set `test=True` to perform a significance test.
	'''
	convert = np.pi / 180.0 if deg else 1
	sx = np.frompyfunc(np.sin, 1, 1)((x - mean(x, deg)) * convert)
	sy = np.frompyfunc(np.sin, 1, 1)((y - mean(y, deg)) * convert)
	r = (sx * sy).sum() / np.sqrt((sx ** 2).sum() * (sy ** 2).sum())

	if test:
		l20, l02, l22 = (sx ** 2).sum(),(sy ** 2).sum(), ((sx ** 2) * (sy ** 2)).sum()
		test_stat = r * np.sqrt(l20 * l02 / l22)
		p_value = 2 * (1 - sp.stats.norm.cdf(abs(test_stat)))
		return tuple(round(v, 7) for v in (r, test_stat, p_value))
	return round(r, 7)


# Phase amplitude coupling analysis
# give it the phase providing band
# the amplitude providing band
# your data and the fs
def PAC(lfp, fs, ppb, apb, pp=1, fig=0): 
	# ppb = [4,8] for theta; apb = [80-125] for gamma

	phase_data, amp_data = getPhaseAmp(lfp, fs, ppb, apb)

	# Now run stats to see what we get:
	# TODO: port this to a variable
	rho, pval, r_2, standard_error = circCorr(phase_data, amp_data)

	if fig: # TODO rename and split to seperate func
		#let's look at a small chunk of our data
		plt.figure(figsize = (15,6));
		#dat_norm = (dat[1:int(fs)*2]-np.mean(dat[1:int(fs)*2]))/np.std(dat[1:int(fs)*2])
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


# Calculate PAC via Ozkurt & Schnitzler, 2011
# TODO: Need to Validate!!! 
# Depreciated? 
def PAC_ozkurt_single(lfp, fs, ppb, apb):
	lo, hi = getPhaseAmp(lfp, fs, ppb, apb)
	return np.abs(np.sum(hi * np.exp(1j * lo))) / (np.sqrt(len(lo)) * np.sqrt(np.sum(hi**2)))

def PAC_ozkurt(lo, hi, f_lo, f_hi, fs, w_lo=3, w_hi=3):

	# TODO: move these gut's to a seperate function... 
	# Compute the phase of the lo signal
	ph = fir1(lo, f_lo[0], f_lo[1], fs, w=w_lo)
	lo = np.angle(sp.signal.hilbert(ph))

	# Compute the amplitude of the hi signal
	ap = fir1(hi, f_hi[0], f_hi[1], fs, w=w_hi)
	hi = np.abs(sp.signal.hilbert(ap))

	return np.abs(np.sum(hi * np.exp(1j * lo))) / (np.sqrt(len(lo)) * np.sqrt(np.sum(hi**2)))


# Calculate PAC using the phase-locking value (PLV) method
# TODO: Need to Validate!!!
def PAC_PLV(lfp, fs, ppb, apb):
	lo, hi = getPhaseAmp(lfp, fs, ppb, apb)
	return np.abs(np.mean(np.exp(1j * (lo - hi))))


# compute the phase of greatest coupling
def couplingPhase(phase_data, amp_data, bin_size=10):
	bins = range(-180, 180+bin_size, bin_size) 
	bins = np.dot(bins, 0.0174532925)

	#filling phase bins with amplitudes
	amps = []
	for x in range(len(bins)-1):
		amps_above_lo_bound = np.where(phase_data >= bins[x])[0]
		amps_below_hi_bound = np.where(phase_data < bins[x+1])[0]
		amps_below_hi_bound = set(amps_below_hi_bound)
		amp_inds_in_this_bin = [amp_val for amp_val in amps_above_lo_bound if amp_val in amps_below_hi_bound]
		amps_in_this_bin = amp_data[amp_inds_in_this_bin]
		amps.append(np.mean(amps_in_this_bin))

	bins = bins[:len(bins)-1]
	amps = (amps-np.mean(amps))/(np.std(amps)) #normalize for clarity

	return amps, bins

# Sweep the PAC algo to see all of the coupling frequencies 
def sweepFreqAmp(lfp, fs, maxFreq=200, maxAmp=200, freqWin=5, ampWin=5):
	# Allocate array
	rho_mat = np.zeros((int(maxFreq/freqWin), int(maxAmp/ampWin)))
	r_2_mat = np.zeros((int(maxFreq/freqWin), int(maxAmp/ampWin)))
	
	xc = 0 # slide a window across the image
	for freq in range(1, maxFreq, freqWin):
		yc = 0
		for amp in range(1, maxAmp, ampWin):
			
			# Compute the PAC...
			phsFreq, ampFreq = [freq, freq+freqWin], [amp, amp+ampWin]
			rho, pval, r_2, standard_error = PAC(lfp, fs, phsFreq, ampFreq, pp=0, fig=0)

			rho_mat[xc, yc] = np.log(rho)
			r_2_mat[xc, yc] = np.log(r_2)
			yc += 1
		xc += 1

	return rho_mat, r_2_mat


###########################################################################################
##  3. Ripple Utils
###########################################################################################


# Finds list of rpl segments above a threshold. 
# TODO: Make second threshold optional
def seg_inds(sig, thresh, thresh_low):
	listosegs, listoinds = [], []
	segsBuff, indsBuff = [], []
	inseg = 0
	for i in range(len(sig)):
		if sig[i] > thresh_low:
			inseg = 1
			segsBuff.append(sig[i])
			indsBuff.append(i)
		elif inseg and (sig[i] <= thresh_low):
			inseg = 0
			# Check for second threshold. Only append if passes
			if not all(i < thresh_low for i in segsBuff):
				listosegs.append(segsBuff)
				listoinds.append(indsBuff)
			segsBuff = []
			indsBuff = []

	return listosegs, listoinds


# Molle ... McNaughton et al 2006 method
def detect_ripples(lfp, freqs, fs, nstd):
	# freqs = [150, 250]
	# nstd = 3
	
	# Filter the signal for ripples and compute the hilbert transform
	lfp_rpl = butter_bandpass_filter(lfp,freqs[0],freqs[1],fs)
	lfp_amp = np.abs(scipy.signal.hilbert(lfp_rpl)); # TODO: Unused? 

	# The threshold for ripple detection was set to 3 SDs above the mean RMS signal. 
	win = 10 # ms;  TODO: write function to compute this from fs (not hardcoded)
	std_rmsW = RMS_window(lfp_rpl, win) # Compute the RMS of the data

	std_rms, mu_rms = np.std(std_rmsW), np.mean(std_rmsW) # Basic stats
	amp_thr = mu_rms + (nstd * std_rms) # threshold 
	amp_thr_low = mu_rms + std_rms # low threshold

	# Pad the signal? 
	padding = np.zeros(len(lfp_rpl) - len(std_rmsW)) # TODO: Add to RMS function? 
	std_rmsW = np.append(padding, std_rmsW) # TODO: Prbably should edit this to both sides... 

	# The beginning and end of a ripple were marked at points at which the RMS signal dropped below 1 SD
	los, loi = seg_inds(std_rmsW, amp_thr, amp_thr_low)

	# provided that these two points were separated by 25–75 ms. 
	# TODO: fix this for fs not at 1000!!!
	tme_thresh = 25 # ms
	los = [ind for ind in los if len(ind) > tme_thresh]
	loi = [ind for ind in loi if len(ind) > tme_thresh]

	# Grab the min val for the timestamp for ripple trough: 
	lot = []
	for seg, ind in zip(los, loi):
		tmpDex = seg.index(min(seg)) # grab the index of the trough. 
		lot.append(ind[tmpDex-1])

	return lfp_rpl, los, loi, lot

# Compute number of ripples per second
# TODO: need a more robhust method, maybe without lfp? 
def ripples_per_sec(lfp, lot, fs): 
	rps = []
	rc = 0
	for i in range(len(lfp)):
		if i % fs != 0:
			if i in lot:
				rc += 1
		else: 
			rps.append(rc)
			rc = 0
	return rps 


