function ms_error_N = objectivefcn1(s_in,U0,V0,N,N0)

S1 = [diag(s_in) zeros(N0,N-N0)];
dKin1 = U0*S1*V0';  

[tree_bin_ofs_N,rms_error_N,epct_N,error_N] = taylor_freq_ofs_forward0(N,N0,dKin1);

ms_error_N = rms_error_N.^2;