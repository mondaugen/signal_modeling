# Load in the frame-by-frame classifications as given by the octave/MATLAB
# script hsrp_test_9.m
#
# Make a set of nodes representing these groups.
#
# Prepare a linear program that will find a good set of paths through these
# nodes.
#
# Plot the results.
#
# PYTHONPATH=$PWD:$PWD/test python test/ptpath_pp_groups_test_1.py
import sys
import ptpath, ptpath_test
import math
import matplotlib.pyplot as plt
J=2 # number of paths
Z=ptpath_test.load_grph_harm_sines_rp(sys.argv[1])
S_=reduce(lambda x,y: x+y,Z['S'])
F_=reduce(lambda x,y: x+y,Z['F'])
D=dict()
for s,f in zip(S_,F_):
    D[f]=s
R=ptpath_test.make_pp_groups(Z,D)
print 'num frames'
print len(R['F'])
# make LP
cf_opt=dict()
cf_opt['H']=Z['opt'][0]['H']
cf_opt['sig']=.1

def P__(x,y):
    return math.sqrt(2.*math.pi)*cf_opt['sig']*math.exp(-0.5*((float(x['w'])
        -float(y['w'])/cf_opt['sig'])**2.))

def P__2(x,y):
    return 1.
    
cf_opt['P']=P__2

def tmp_cost_func(a,b):
    return ptpath_test.LPNode_ppg_dist(a,b,cf_opt)

d=ptpath.g_f_2lp(R['S'],R['F'],J,cost_func=tmp_cost_func,opt={'calc_mean':0,'min_mean_dev':0})
sol=ptpath_test.sol_test_1(d)
paths=ptpath_test.lp_sol_extract_paths(sol,R['S'],R['F'])
ptpath_test.pp_groups_plot_paths(R['S'],paths,cf_opt,show=False,fignum=0)
ptpath_test.plot_hsrpc_test(Z,D,show=False,fignum=1)
plt.figure(0)
plt.title('Source tracking using LP')
plt.xlabel('Sample Number')
plt.ylabel('Frequency in Radians / Sample')
plt.figure(1)
plt.title('Raw frame-by-frame classification')
plt.xlabel('Sample Number')
plt.ylabel('Frequency in Radians / Sample')
plt.show()
#ptpath_test.lp_sol_plot_paths(sol,R['S'],R['F'])
