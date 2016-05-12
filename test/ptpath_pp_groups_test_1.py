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
J=3 # number of paths
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
cf_opt['P']=lambda x,y: 1

def tmp_cost_func(a,b):
    return ptpath_test.LPNode_ppg_dist(a,b,cf_opt)

d=ptpath.g_f_2lp(R['S'],R['F'],J,cost_func=tmp_cost_func,opt={'calc_mean':0,'min_mean_dev':0})
sol=ptpath_test.sol_test_1(d)
paths=ptpath_test.lp_sol_extract_paths(sol,R['S'],R['F'])
#ptpath_test.pp_groups_plot_paths(R['S'],paths,cf_opt,show=True)
ptpath_test.lp_sol_plot_paths(sol,R['S'],R['F'])
