import ptpath, ptpath_test
import matplotlib.pyplot as plt
J=4
g=ptpath_test.rand_graph_1(6,6,10)
(C,C_cxns)=ptpath_test.shortest_paths_cost_lattice(g['S'],g['F'],J)
print 'Solving via Viterbi'
q=ptpath_test.shortest_paths_viterbi(C,C_cxns)
print 'Done'
ptpath_test.plot_spv(g['S'],g['F'],q,C_cxns,show=False,fignum=0)
d=ptpath.g_f_2lp(g["S"],g["F"],J,opt={'calc_mean':0,'min_mean_dev':0})
print 'Solving via LP'
sol=ptpath_test.sol_test_1(d)
print 'Done'
ptpath_test.plot_test(sol,g["S"],show=False,fignum=1)
plt.show()
