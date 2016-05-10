# Compare using the more descriptive parameters (frequency slope etc.) in the
# cost function with the simple Euclidean distance between frequency parameters
# invoke with
# PYTHONPATH=$PWD:$PWD/test python ptpath_hsrp_cmp_1.py
import ptpath, ptpath_test
Z=ptpath_test.load_grph_harm_sines_rp("../masters_thesis/data/hsrp_20160504T161736928.mat")
S_=reduce(lambda x,y: x+y,Z['S'])
F_=reduce(lambda x,y: x+y,Z['F'])
D=dict()
for s,f in zip(S_,F_):
    D[f]=s
print 'Calculating d1'
d1=ptpath.g_f_2lp(D,Z['F'],10,ptpath_test.LPNode_rp_dist,{"calc_mean":0,"min_mean_dev":0})
print 'Calculating d2'
d2=ptpath.g_f_2lp(D,Z['F'],10,ptpath_test.LPNode_rp_eu_dist,{"calc_mean":0,"min_mean_dev":0})
print 'Solving d1'
sol1=ptpath_test.sol_test_1(d1)
print 'Solving d2'
sol2=ptpath_test.sol_test_1(d2)
ptpath_test.plot_lp_hsrp_cmp(sol1,sol2,D,Z['trues'])
