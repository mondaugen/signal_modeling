import ptpath, ptpath_test
import matplotlib.pyplot as plt
J=3
g=ptpath_test.rand_graph_1(3,5,15)
(C,C_cxns)=ptpath_test.shortest_paths_cost_lattice(g['S'],g['F'],J)
print 'Solving via Viterbi'
q=ptpath_test.shortest_paths_viterbi(C,C_cxns)
print 'Done'
ptpath_test.plot_spv(g['S'],g['F'],q,C_cxns,show=False,fignum=0)
plt.title('Shortest paths using Viterbi algorithm')
d=ptpath.g_f_2lp(g["S"],g["F"],J,opt={'calc_mean':0,'min_mean_dev':0})
print 'Solving via LP'
sol=ptpath_test.sol_test_1(d)
print 'Done'
ptpath_test.plot_test(sol,g["S"],show=False,fignum=1)
plt.title('Shortest paths using LP algorithm')
plt.show()
