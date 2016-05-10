import ptpath, ptpath_test
g=ptpath_test.rand_graph_1(10,20,50)
q=ptpath_test.shortest_eu_path_viterbi(g['S'],g['F'])
ptpath_test.plot_sepv(g['S'],g['F'],q)
