./bin/SID_ERMS_Rescore.default.linuxgccdebug \
    -database ../database \
    -in::file::l ../inputs/inputs/in_list \
    -out:file:o ../inputs/test_output.out \
    -complex_type ../inputs/C2.sym \
    -ERMS ../inputs/exp_ERMS.csv \
    -RMSE true 
