# Compare using the more descriptive parameters (frequency slope etc.) in the
# cost function with the simple Euclidean distance between frequency parameters
# invoke with
# PYTHONPATH=$PWD:$PWD/test python ptpath_hsrp_cmp_1.py
import sys
import ptpath, ptpath_test
Z=ptpath_test.load_grph_harm_sines_rp(sys.argv[1])
S_=reduce(lambda x,y: x+y,Z['S'])
F_=reduce(lambda x,y: x+y,Z['F'])
D=dict()
for s,f in zip(S_,F_):
    D[f]=s
ptpath_test.plot_hsrpc_test(Z,D)
