# invoke from root of repo as
# bash test/plot_all_am_fm_sep.sh
for i in {'0.01','0.001'}
do
    for j in {'0.0001','0.00001'}
    do
        octave -qf test/hsrp_test_7_save_data_specvar.m ${i} ${j}
        PYTHONPATH=$PWD:$PWD/test python \
            test/plot_hsrp_test_7_dat_specfile_lp_path.py \
            hsrp_test_7_${i}_${j}.dat
    done
done
