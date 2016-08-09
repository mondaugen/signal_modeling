# Show and print spectrograms of the original files and the source separated
# ones
import numpy as np
import matplotlib.pyplot as plt
import os
import neplot as nep
# Color contrast config
# values further from 1, more contrast
clr_gamma=4.
clr_mapper=nep.PowerNormalize(clr_gamma)

plotoutpath=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/'

plt.rc('text',usetex=True)
plt.rc('font',family='serif')

fpath_ac_gtr_orig='tmp/ac_gtr_a3_op_sr16k.f64'
fpath_ac_gtr_ss='tmp/xylo_fs4_ac_gtr_a3_sr16k_p0.f64'
fpath_xylo_orig='tmp/xylo_fs4_sr16k.f64'
fpath_xylo_ss='tmp/xylo_fs4_ac_gtr_a3_sr16k_p1.f64'

with open(fpath_ac_gtr_orig,'r') as f:
    ac_gtr_orig=np.fromfile(f,'float64')
with open(fpath_ac_gtr_ss,'r') as f:
    ac_gtr_ss=np.fromfile(f,'float64')
with open(fpath_xylo_orig,'r') as f:
    xylo_orig=np.fromfile(f,'float64')
with open(fpath_xylo_ss,'r') as f:
    xylo_ss=np.fromfile(f,'float64')

# Hop size
H=512
# Sample rate of file
Fs=16000
# Analysis window size and FFT size
M=4096

plt.figure(1)
plt.specgram(ac_gtr_orig,M,Fs,noverlap=M-H,norm=clr_mapper,cmap='Greys')
tmp_title='Spectrogram of acoustic guitar'
plt.xlabel('Time (seconds)')
plt.ylabel('Frequency (Hz)')
plt.xlim(0,(len(ac_gtr_orig)-H)/Fs)
plt.savefig(plotoutpath+'ac_gtr_orig_spec.eps')
with open(plotoutpath+'ac_gtr_orig_spec.txt','w') as f:
    f.write(tmp_title+'%')
plt.figure(2)
plt.specgram(ac_gtr_ss,M,Fs,noverlap=M-H,norm=clr_mapper,cmap='Greys')
tmp_title='Spectrogram of source separated acoustic guitar'
plt.xlabel('Time (seconds)')
plt.ylabel('Frequency (Hz)')
plt.xlim(0,(len(ac_gtr_ss)-H)/Fs)
plt.savefig(plotoutpath+'ac_gtr_ss_spec.eps')
with open(plotoutpath+'ac_gtr_ss_spec.txt','w') as f:
    f.write(tmp_title+'%')
plt.figure(3)
plt.specgram(xylo_orig,M,Fs,noverlap=M-H,norm=clr_mapper,cmap='Greys')
tmp_title='Spectrogram of xylophone'
plt.xlabel('Time (seconds)')
plt.ylabel('Frequency (Hz)')
plt.xlim(0,(len(xylo_orig)-H)/Fs)
plt.savefig(plotoutpath+'xylo_orig_spec.eps')
with open(plotoutpath+'xylo_orig_spec.txt','w') as f:
    f.write(tmp_title+'%')
plt.figure(4)
plt.specgram(xylo_ss,M,Fs,noverlap=M-H,norm=clr_mapper,cmap='Greys')
tmp_title='Spectrogram of source separated xylophone'
plt.xlabel('Time (seconds)')
plt.ylabel('Frequency (Hz)')
plt.xlim(0,(len(xylo_ss)-H)/Fs)
plt.savefig(plotoutpath+'xylo_ss_spec.eps')
with open(plotoutpath+'xylo_ss_spec.txt','w') as f:
    f.write(tmp_title+'%')

plt.show()
