#Plot recent data generated by hsrp_test_9.m
# Invoke from parent directory as
# bash ht9_plt_grps.sh
fname=`octave -q hsrp_test_9.m`
echo ${fname}
PYTHONPATH=$PWD/..:$PWD python ptpath_pp_groups_test_2.py ${fname}
