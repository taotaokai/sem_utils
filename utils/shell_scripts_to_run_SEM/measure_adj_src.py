#!/usr/bin/env python

import numpy as np
from scipy.fftpack import shift 
from scipy.signal import butter, filtfilt 
#import matplotlib as mpl
#import matplotlib.pyplot as plt
from obspy import read

####################
def tukeywin(npts, r):
    win = np.ones(npts)
    low = int(r*(npts-1)/2+0.5)
    high = int((1-r/2)*(npts-1))
    for n in np.arange(0,low):
        Z = np.pi * (2*n/r/(npts-1)-1)
        win[n] = 0.5 + np.cos(Z)/2
    for n in np.arange(high,npts):
        Z = np.pi * (2*n/r/(npts-1)-2/r+1)
        win[n] = 0.5 + np.cos(Z)/2
    return win

####################
def butter_bandpass(data, lowcut, highcut, fs, order=2):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    data_filt = filtfilt(b, a, data) # use zero phase
    return data_filt 

####################
DEG2RAD = np.pi/180
def en2rt(en, baz):
    '''rotate e/n to r/t (r-t-z right-hand axes).
    Parameters:
        en: input E-/N- components, size(en)=(2,npts),
            E- first row, N- second row;
        baz: back azimuth
    Return value: R-/T- components .'''
    # rotation matrix
    te = np.cos(baz*DEG2RAD)
    tn = np.cos((baz+90)*DEG2RAD)
    re = tn
    rn = -te
    rot = np.array([[re,rn],
                    [te,tn]])
    return np.dot(rot,en)

#-------------------
def rt2en(rt, baz):
    '''rotate r/t to e/n (r-t-z right-hand axes).
    Parameters:
        rt: input R-/T- components, size(rt)=(2,npts),
            R- first row, T- second row;
        baz: back azimuth
    Return value: E-/N- components .'''
    # rotation matrix
    te = np.cos(baz*DEG2RAD)
    tn = np.cos((baz+90)*DEG2RAD)
    re = tn
    rn = -te
    rot = np.array([[re,te],
                    [rn,tn]])
    return np.dot(rot,rt)

####################
def measure_adj_src(st_obs, st_syn, freqmin, freqmax, phase_wins):
    '''measure adjoint source using correlation coefficient as
       the objective function. 
    
    Parameters:
        (obspy.stream) st_obs: seismograms [0:2] -> e/n/z
                       st_syn: st_synthetics  [0:2] -> e/n/z
        freqmin,freqmax: frequency range;
        phase_wins:({name, ttime, begin, end, comp},{},...) phase windows
    
    Return value:
        (obspy.stream) st_adj: adjoint sources [0:4] -> e/n/z/r/t.'''

    # band-pass filter
    for icomp in range(0, 3):
        data = st_syn[icomp].data
        fs = st_syn[icomp].stats.sampling_rate
        st_syn[icomp].data = butter_bandpass(data,freqmin,freqmax,fs)
        data = st_obs[icomp].data
        fs = st_obs[icomp].stats.sampling_rate
        st_obs[icomp].data = butter_bandpass(data,freqmin,freqmax,fs)

    # shift st_obs so the first sample lies at the same time as in st_syn
    for icomp in range(0, 3):
        t0_st_syn = st_syn[icomp].stats.starttime
        t0_st_obs = st_obs[icomp].stats.starttime
        t_shift = t0_st_syn - t0_st_obs
        #print "shift st_obs backward ", t_shift
        npts = st_obs[icomp].stats.npts
        delta = st_obs[icomp].stats.delta
        tlen = npts*delta
        data_shift = shift(st_obs[icomp].data,t_shift,tlen)
        st_obs[icomp].data = data_shift
    # resample st_obs as the same sampling rate as st_syn
    fs = st_syn[0].stats.sampling_rate
    st_obs.resample(fs)

    # copy seismogram data into matrix 
    npts = st_syn[0].stats.npts
    obs_enzrt = np.zeros((5,npts))
    syn_enzrt = np.zeros((5,npts))
    for icomp in range(0, 3):
        obs_enzrt[icomp,:] = st_obs[icomp].data[0:npts]
        syn_enzrt[icomp,:] = st_syn[icomp].data
    # get rotated seismograms
    baz = st_syn[0].stats.sac.baz
    obs_enzrt[3:5,:] = en2rt(obs_enzrt[0:2,:], baz)
    syn_enzrt[3:5,:] = en2rt(syn_enzrt[0:2,:], baz)
    # define the time samples of syn
    bt = st_syn[0].stats.sac.b
    ot = st_syn[0].stats.sac.o
    delta = st_syn[0].stats.delta
    npts = st_syn[0].stats.npts
    ts = bt - ot + delta*np.arange(0,npts) # zero corresponds to the origin time

    # measure adj_src in each time window
    adj_enzrt = np.zeros((5,npts))
    weight = np.zeros(npts)
    for ph in phase_wins:
        # set data window
        idx = (ts>=ph['ttime']+ph['begin']) & \
              (ts<=ph['ttime']+ph['end'])
        npts_win = sum(idx)
        window = tukeywin(npts_win,0.2)
        # record the weighting
        weight[idx] = weight[idx]+window
        # select components
        comp = ph['comp']
        obs_enzrt_win = np.zeros((5,npts_win))
        syn_enzrt_win = np.zeros((5,npts_win))
        if ('z' in comp):
            obs_enzrt_win[2,:] = obs_enzrt[2,idx]
            syn_enzrt_win[2,:] = syn_enzrt[2,idx]
        if (('r' in comp) and ('t' in comp) or
            ('e' in comp) and ('n' in comp)):
            obs_enzrt_win[0:2,:] = obs_enzrt[0:2,idx]
            syn_enzrt_win[0:2,:] = syn_enzrt[0:2,idx]
        elif (('t' in comp) and ('r' not in comp) and 
              ('e' not in comp) and ('n' not in comp)):
            obs_enzrt_win[4,:] = obs_enzrt[4,idx]
            syn_enzrt_win[4,:] = syn_enzrt[4,idx]
            obs_enzrt_win[0:2,:] = rt2en(obs_enzrt_win[3:5,:],baz)
            syn_enzrt_win[0:2,:] = rt2en(syn_enzrt_win[3:5,:],baz)
        else:
            print "Error: invalid 'comp' value!"
            print "phase = ", ph
            exit()
        # apply the window taper
        obs_enzrt_win = obs_enzrt_win*window
        syn_enzrt_win = syn_enzrt_win*window
        # calculate adj_src
        norm2_obs = np.sum(obs_enzrt_win[0:3,:]**2, axis=(0,1))
        norm2_syn = np.sum(syn_enzrt_win[0:3,:]**2, axis=(0,1))
        cc = np.sum(obs_enzrt_win[0:3,:]*syn_enzrt_win[0:3,:], axis=(0,1))
        amp_ratio = cc/norm2_syn
        norm = np.sqrt(norm2_syn*norm2_obs)
        cc = cc/norm
        #print "cc=", cc, "amp_ratio(obs/syn)=", amp_ratio
        #print "norm=", norm
        adj_enzrt_win = np.zeros((5,npts_win))
        #### KEY FORMULAR HERE ####
        # adj_src = (obs - A*syn)/||obs||/||syn||
        adj_enzrt_win[0:3,:] = (obs_enzrt_win[0:3,:] - 
                                amp_ratio*syn_enzrt_win[0:3,:])/norm * fs
        # source origin time sensitivity
        #syn_enzrt_win_fft = np.zeros((5,npts_win))
        #syn_enzrt_win_fft[0:3,:] = np.fft.fft(syn_enzrt_win[0:3,:])
        #freqs = np.fft.fftfreq(npts_win,d=delta)
        #r = sum(abs(freqs)>freqmax)/npts_win
        #freq_weight = np.fft.fftshift(tukeywin(npts_win,r))
        #plt.plot(freqs,freq_weight)
        #plt.show()
        #syn_enzrt_win_fft[0:3,:] = syn_enzrt_win_fft[0:3,:]*2j*np.pi*freqs
        #syn_enzrt_win_fft[0:3,:] = np.fft.ifft(syn_enzrt_win_fft[0:3,:]*freq_weight)
        #ker_orig = -1 * np.sum(syn_enzrt_win_fft[0:3,:]*adj_enzrt_win[0:3,:], axis=(0,1))
        #print "ker_orig=", ker_orig

        ###########################
        adj_enzrt_win[3:5,:] = en2rt(adj_enzrt_win[0:2,:],baz)

        # put adj_enzrt_win into adj_enzrt
        adj_enzrt[:,idx] = adj_enzrt[:,idx] + adj_enzrt_win

        # output info
        netwk = st_syn[0].stats.network 
        stnm  = st_syn[0].stats.station 
        stla  = st_syn[0].stats.sac.stla 
        stlo  = st_syn[0].stats.sac.stlo 
        stel  = st_syn[0].stats.sac.stel 
        gcarc = st_syn[0].stats.sac.gcarc 
        baz   = st_syn[0].stats.sac.baz
        print "{:5s} {:5s} {:7.3f} {:8.3f} {:6.1f} | {:5.1f} {:5.1f} | {} {} {} {:3s} {:4.2f} {:5.2f}".format(netwk,stnm,stla,stlo,stel,gcarc,baz,ph['name'],ph['begin'],ph['end'],comp,cc,amp_ratio)
    #end for ph

    idx = weight>1.
    weight[idx] = 1/weight[idx]
    adj_enzrt = adj_enzrt*weight

    # return st_adj
    st_adj = st_syn.copy()
    st_syn.rotate('NE->RT',baz)
    st_adj.append(st_syn[1])
    st_adj.append(st_syn[0])
    st_adj[0].data = adj_enzrt[0,:]
    st_adj[1].data = adj_enzrt[1,:]
    st_adj[2].data = adj_enzrt[2,:]
    st_adj[3].data = adj_enzrt[3,:]
    st_adj[4].data = adj_enzrt[4,:]

    return st_adj, ts
#END def measure_adj_src(st_syn, st_obs, freqlim, phase_wins):

####################
if __name__ == '__main__':

    import sys
    from obspy.taup.taup import getTravelTimes

    # read arguments
    argv = sys.argv
    if len(argv) != 8:
        print "Usage: measure_adj_src.py stalst phalst obsdir syndir adjdir freqmin freqmax"
        sys.exit()
    stalst = argv[1]
    phalst = argv[2]
    obsdir = argv[3]
    syndir = argv[4]
    adjdir = argv[5]
    freqmin = float(argv[6])
    freqmax = float(argv[7])
    #----------
    f_sta = open(stalst,"r")
    try: 
        for sta in f_sta: 
            # skip commented line
            if sta.strip().startswith('#'):
                continue
            #print sta.strip()
            sta = sta.split()
            stnm = sta[0]
            netwk = sta[1]
            # read obs
            fn_obs = obsdir + '/' + netwk + '.' + stnm + '.*[ENZ]'
            st_obs = read(fn_obs)
            # read syn
            fn_syn = syndir + '/' + netwk + '.' + stnm + '.*[ENZ].sem.sac'
            st_syn = read(fn_syn)
            # get travel time information 
            gcarc = st_syn[0].stats.sac.gcarc
            evdp = st_syn[0].stats.sac.evdp
            #print gcarc, evdp
            ttaup = getTravelTimes(gcarc, evdp, 'ak135')
            # calculate phase windows
            f_pha = open(phalst,"r")
            phase_wins = list({})
            try:
                for pha in f_pha:
                    # skip commented line
                    if pha.strip().startswith('#'):
                        continue
                    pha = pha.split()
                    pha_inf = dict(name=pha[0], ttime=0, 
                                   begin=float(pha[1]), 
                                   end=float(pha[2]), comp=pha[3])
                    # make lists of phase windows
                    for tt in ttaup:
                        if tt['phase_name'] == pha[0]:
                            pha_inf['ttime'] = tt['time']
                            phase_wins.append(pha_inf)
                            break
            finally:
                f_pha.close()
            # calculate adj_src
            st_adj, ts = measure_adj_src(st_obs, st_syn, 
                                            freqmin, freqmax, phase_wins)
            # output
            for tr in st_adj:
                tr.write(adjdir + '/' + netwk + '.' + stnm + '.' + 
                         tr.stats.channel + ".adj.sac", "sac")
                #ascii for SEM
                f = open(adjdir + '/' + netwk + '.' + stnm + '.' +
                         tr.stats.channel + ".adj", "w")
                for i in xrange(len(tr.data)):
                    f.write("{:16.9e}  {:16.9e}\n".format(ts[i],tr.data[i]))
                f.close
    finally:
        f_sta.close()
#END main
