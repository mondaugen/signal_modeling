import ptpath, ptpath_test
g=ptpath_test.rand_graph_1(6,6,10)
(C,C_cxns)=ptpath_test.shortest_paths_cost_lattice(g['S'],g['F'],5)
q=ptpath_test.shortest_paths_viterbi(C,C_cxns)
ptpath_test.plot_spv(g['S'],g['F'],q,C_cxns)
