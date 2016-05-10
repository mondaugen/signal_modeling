import ptpath, ptpath_test
g=ptpath_test.rand_graph_1(4,4,30)
print 'done gen graph'
d=ptpath.g_f_2lp(g["S"],g["F"],3,opt={'calc_mean':0,'min_mean_dev':0})
print 'done gen program'
sol=ptpath_test.sol_test_1(d)
print 'size A'
print d['A'].size
print 'size G'
print d['G'].size
ptpath_test.plot_test(sol,g["S"])
