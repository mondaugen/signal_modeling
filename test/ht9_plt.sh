#Plot recent data generated by hsrp_test_9.m
# Invoke from parent directory as
# bash ht9_plt.sh
fname=`octave hsrp_test_9.m`
echo ${fname}
PYTHONPATH=$PWD/..:$PWD python ptpath_hsrpc_plot_test_1.py ${fname}
