%
% run parameters
%

%class_name = 'double'; %plot_ID=1;
class_name = 'single'; %plot_ID=2;

Nt_pfb0 = 1024;

M0 = 2*Nt_pfb0;   % number of subfilters = FFT size
M0 = 4*Nt_pfb0;   % number of subfilters = FFT size
M0 = 16*Nt_pfb0;   % number of subfilters = FFT size
%M0 = 32*Nt_pfb0;   % number of subfilters = FFT size

fs_in = 1.0*M0;
n_in = Nt_pfb0*M0;

fs0 = fs_in/M0;
Ts0 = 1/fs0;
t_max = Nt_pfb0/fs0;
bin_bw = fs_in/M0;  % nominal bandwidth of PFB bins

df_dt_min = 0;
df_dt_max =  fs0.^2;

PlotLineWidth = 1;
AxisFontSize = 12;
AxisFontWeight = 'bold';
EnableFinalPlots = 0;

plot_taylor_det_peaks = 0;
plot_fastDD_det_peaks = 0;
plot_taylor_shifted_spgram = 0;
plot_fastDD_shifted_spgram = 0;
plot_fastDD_det_plane = 0;
plot_fastDD_det_profile = 0;
plot_taylor_det_plane = 0;
plot_taylor_det_profile = 0;
plot_multichirp_sg_rfi = 1;  % combined input and output SG
plot_multichirp_sg = 0;
plot_multichirp_cat_sg1 = 0;
plot_multichirp_cat_sg2 = 0;  % zoomed version
plot_ID = 1;

%ver_DD = 3;
ver_DD = 4; % fast mex DD versions
use_mex = 1; % =1 to use mex pfb and DD

use_weighted_mean = 1;  % use Gaussian-weighted mean
%use_weighted_mean = 0;  % use uniformly-weighted mean

%N_tap_list = [3 4 8];
%N_tap_list = 1;
%N_tap_list = 4;
N_tap_list = 3;

Lf_list = [1 2 4];
N_pre_DD_nom = 1;    % number of SG integrations prior to DD, must be power of 2
rfi_proc_enable = 0;
rfi_alg_ID = 1;
f1_rfi = []; A_rfi = [];

use_max_Lr = 1;  % =1 to set Lr=Lf, =0 use Lr_list
Lr_list = 2;

%window1_list = {'RectWin','Hann'};
window1_list = {'RectWin'};
%window2_list = {'lpfxxx','Hann'};
window2_list = {'lpfxxx'};
%Mx_list = [1 2];  % multiplier to initial M0 value
Mx_list = [1];  % nominal multiplier to initial M0 value

alg_baseline = 'Taylor';
N0_baseline = 4;
N_tap_baseline = 1;
Lf_baseline = 1;
Lr_baseline = 1;
N_pre_DD_baseline = 1; % number of SG integrations prior to DD, must be power of 2
%fastDD_alg_ID = 1;
fastDD_alg_ID = 3;

snr_db_list = [Inf];  % zero noise analysis mode
DT_std_dev = 5;

N0_list = [16];

%N0_taylor0 = [];  % if empty, use same value as N0_list
N0_taylor0 = [8]; % use one value for taylor

m0_chirp_list = [10:5:1010];
m0_offset_list = NaN;  % random offsets
f_offset_list = NaN;  % random offsets

%m0_offset_list = 0;

%f_offset_list = 0;
%f_offset_list = 0.25*bin_bw;
%f_offset_list = 0.5*bin_bw;
%f_offset_list = [0 .25 .5]*bin_bw;
%f_offset_list = [0:.125:.375]*bin_bw;
%f_offset_list = [.125 .375]*bin_bw;

mod_type = []; n_mod_ID = 1;
%mod_type='BPSK'; n_mod_ID = 10;
%mod_type='AM';   n_mod_ID = 10;

pos_neg = 1;  % =1 for positive drift rates
              % =-1 for negative drift rates
              % =0 for pos & neg drift rates

max_n_diff = 2;

run_string = '';
run_string = run_id;

if strcmp(run_id,'test')
  snr_db_list = [-3:3];  % statistical detection mode
  snr_db_list = Inf;  % statistical detection mode
  m0_chirp_list = [10:5:1010];
  plot_ID = 1;
  if (1)
    % single run for mex testing
    N_tap_list = [3];
    use_max_Lr = 1;  % =1 to set Lr=Lf, =0 use Lr_list
    Lf_list = [2];
    N0_taylor0 = [4]; % use one value for taylor
  elseif (1)
    % repeated runs for profiling
    N_tap_list = [3];
    use_max_Lr = 1;  % =1 to set Lr=Lf, =0 use Lr_list
    Lf_list = [1]*ones(1,10); N_tap_list = [1];
    %Lf_list = [1]*ones(1,10);
    %Lf_list = [3/2]*ones(1,10);
    %Lf_list = [7/4]*ones(1,10);
    %Lf_list = [2]*ones(1,10);
    %Lf_list = [3]*ones(1,10);
    N0_taylor0 = [4]; % use one value for taylor
  elseif (1)
    % Mx=2 N_tap=1 case - works poorly
    N_tap_list = [1];
    use_max_Lr = 1;  % =1 to set Lr=Lf, =0 use Lr_list
    Lf_list = [1 5/4 3/2 7/4 2 3 4];
    window1_list = {'RectWin','Hann'};
    Mx_list = 2;
  end
  run_string = '_test';
elseif contains(run_id,'baseline')
  % 1x-2x-4x Hann N_tap=4-8 baseline run - Inf SNR
  plot_ID = 1;
  %plot_ID = 2;  % TEMP
  %plot_ID = 10*ver_DD;
  Lf_list = [1 5/4 3/2 7/4 2 3 4];
  Lr_list = [1 5/4 3/2 7/4 2 3 4];
  N0_list = [8];
  N0_taylor0 = [8]; % use one value for taylor
  pos_neg = 1;
  snr_db_list = [Inf];  % zero noise analysis mode
  m0_chirp_list = [10:5:1010];
  if contains(run_id,'baseline-inf1') % baseline-inf1a <-> baseline-inf1d
    Mx_list = 1; plot_id = 10; % multiplier to initial M0 value
    if strcmp(run_id(end),'a'), Mx_list = 1; N_tap_list = [1]; plot_id = 10;
    elseif strcmp(run_id(end),'b'), Mx_list = 1; N_tap_list = [3]; plot_id = 20;
    elseif strcmp(run_id(end),'c'), Mx_list = 2; N_tap_list = [1]; plot_id = 30;
    elseif strcmp(run_id(end),'d'), Mx_list = 2; N_tap_list = [3]; plot_id = 40;
    end
    use_max_Lr = 1;  % =1 to set Lr=Lf, =0 use Lr_list
    Lf_list = [1 5/4 3/2 7/4 2 3 4];
    window2_list = {'lpfxxx'};
    run_string = 'baseline-inf1';
  elseif contains(run_id,'baseline-inf2') % baseline-inf2a <-> baseline-inf2c
    % sweep N_tap
    Mx_list = 1; plot_id = 10;
    if strcmp(run_id(end),'a'), Mx_list = 1; plot_id = 10;
    elseif strcmp(run_id(end),'b'), Mx_list = 2; plot_id = 20;
    elseif strcmp(run_id(end),'c'), Mx_list = 0.5; plot_id = 30;
    end
    use_max_Lr = 1;  % =1 to set Lr=Lf, =0 use Lr_list
    Lf_list = [1 3/2 2 3 4];
    N_tap_list = [3 4 6 8];
    window2_list = {'lpfxxx'};
    run_string = 'baseline-inf2';
  elseif contains(run_id,'baseline-inf3')
    % compare Hann lpfxxx and RectWin and lpfxxx
    use_max_Lr = 1;  % =1 to set Lr=Lf, =0 use Lr_list
    Lf_list = [1 3/2 2 3 4];
    if strcmp(run_id(end),'a'), window2_list = {'lpfxxx','Hann'}; N_tap_list = [3]; plot_id = 10;
    elseif strcmp(run_id(end),'b'), window2_list = {'lpfxxx'}; N_tap_list = [1 3]; plot_id = 20;
    end
    run_string = 'baseline-inf3';
  elseif strcmp(run_id,'baseline-inf4')
    % examine cases for Lf=Lr and Lf~=Lr
    use_max_Lr = 0;  % =1 to set Lr=Lf, =0 use Lr_list
    Lf_list = [1 3/2 7/4 2 2 3 4 4];
    Lr_list = [1 3/2 7/4 1 2 3 2 4];
    N_tap_list = [3];
    window2_list = {'lpfxxx'};
  elseif strcmp(run_id,'baseline-inf5')
    % examine fs=1 and fs=0.5
    use_max_Lr = 1;  % =1 to set Lr=Lf, =0 use Lr_list
    Mx_list = [1 2];
    Lf_list = [1 3/2 2];
    N_tap_list = [1 3];
    window2_list = {'lpfxxx'};
    N0_list = [8];
    N0_taylor0 = [];  % if empty, use same value as N0_list
  elseif strcmp(run_id,'baseline-inf6')
    use_max_Lr = 1;  % =1 to set Lr=Lf, =0 use Lr_list
    Lf_list = [1 3/2 2 3 4];
    N_tap_list = [3];
    m0_chirp_list = [10:5:1010];
    m0_offset_list = 0;  % no random offsets
    f_offset_list = 0; 
    snr_db_list = [Inf];  % zero noise analysis mode
    plot_taylor_det_peaks = 0;
    plot_fastDD_det_peaks = 1;
    plot_multichirp_sg = 1;
    plot_multichirp_cat_sg2 = 1;  % zoomed version
  elseif contains(run_id,'baseline-snr1') % baseline-snr1a <-> baseline-snr1b
    % sweep N_tap
    Mx_list = 1; plot_id = 10;
    if strcmp(run_id(end),'a'), Mx_list = 1; plot_id = 10;
    elseif strcmp(run_id(end),'b'), Mx_list = 2; plot_id = 20;
    end
    use_max_Lr = 1;  % =1 to set Lr=Lf, =0 use Lr_list
    Lf_list = [1 3/2 2 3 4];
    N_tap_list = [1 3];
    m0_chirp_list = [10:5:1010];
    window2_list = {'lpfxxx'};
    snr_db_list = [-5:5];  % statistical detection mode
    run_string = 'baseline-snr1';
  elseif strcmp(run_id,'baseline-snr2')
    use_max_Lr = 1;  % =1 to set Lr=Lf, =0 use Lr_list
    Lf_list = [4];
    N_tap_list = [3];
    %m0_chirp_list = [10:17:1010];
    m0_chirp_list = [10:33:1010];
    snr_db_list = [-5:5];  % statistical detection mode
    window2_list = {'lpfxxx'};
    plot_fastDD_det_plane = 1;
    plot_taylor_det_plane = 0;
    plot_fastDD_det_profile = 1;
    plot_taylor_det_profile = 0;
  end
elseif contains(run_id,'N0-sweep32-Mx')
  Mx_list = str2num(run_id(14:end));
  snr_db_list = [Inf];  % zero noise analysis mode
  N0_list = [4 8 16 32];
  N0_taylor0 = [];  % if empty, use same value as N0_list
  use_max_Lr = 1;  % =1 to set Lr=Lf, =0 use Lr_list
  N_tap_list = [1 3];
  Lf_list = [1 3/2 2 3 4];
  run_string = 'N0-sweep32-Mx';
elseif strcmp(run_id,'N0-sweep')
  snr_db_list = [Inf];  % zero noise analysis mode
  N0_list = [4 8 16 32 64 128 256 512 1024];
  N0_taylor0 = [];  % if empty, use same value as N0_list
  use_max_Lr = 1;  % =1 to set Lr=Lf, =0 use Lr_list
  N_tap_list = [1 3];
  plot_ID = 1;
  if (1)
    N_tap_list = [1 4 8];
    Lf_list = [1 3/2 2 4];
  else
    N_tap_list = [4];
    Lf_list = [2];
  end
elseif strcmp(run_id(1:end-1),'alg-test')
  fastDD_alg_ID = str2num(run_id(end));
  plot_ID = 10*fastDD_alg_ID;
  snr_db_list = [-10:10];  % statistical detection mode
  snr_db_list = Inf;  % statistical detection mode
  use_max_Lr = 1;  % =1 to set Lr=Lf, =0 use Lr_list
  if (1)
    use_max_Lr = 1;  % =1 to set Lr=Lf, =0 use Lr_list
    N_tap_list = [3];
    Lf_list = [1 3/2 2 3 4];
    Mx_list = [1 2];  % multiplier to initial M0 value
  end
  pos_neg = 1;
  run_string = '_alg-test';
elseif strcmp(run_id(1:4),'rfi-')
  rfi_inr_db = str2num(run_id(end-7:end-6));
  N_avg1 = str2num(run_id(end-5:end-3));
  if (N_avg1==999)||(N_avg1==0), N_avg1 = Nt_pfb0; end
  rfi_duty_idx = str2num(run_id(end-2:end-1));
  rfi_proc_enable = str2num(run_id(end));
  plot_ID = 10*str2num(run_id(end-7:end));

  mod_params_rfi.T_start = 0;
  if rfi_duty_idx==0
    rfi_duty_pct = 100;
    mod_params_rfi.T_on = t_max;
    mod_params_rfi.T_off = 0;
  else
    rfi_duty_pct = rfi_duty_idx;
    mod_params_rfi.T_on = .125*t_max*rfi_duty_pct/100;
    mod_params_rfi.T_off = .125*t_max*(1-rfi_duty_pct/100);
  end
  if (rfi_inr_db>0)
    A_rfi = 10.^(.1*rfi_inr_db)*100/rfi_duty_pct;
  else
    rfi_inr_db = -Inf;
    A_rfi = 0;
  end

  rfi_alg_ID = rfi_proc_enable;
  rfi_bin_limit = 12;
  rfi_bin_replace = 2.0;
  
  N_rfi_avg_t = N_avg1; N_rfi_avg_f = 1;
  if (1)
    if (rfi_alg_ID==0)||(rfi_alg_ID==1), rfi_weight_exp = 1;
    elseif (rfi_alg_ID==2),  rfi_weight_exp = 1.25;
    elseif (rfi_alg_ID==3),  rfi_weight_exp = 1.5;
    elseif (rfi_alg_ID==4),  rfi_weight_exp = 1.75;
    elseif (rfi_alg_ID==5),  rfi_weight_exp = 2;
    elseif (rfi_alg_ID==6),  rfi_weight_exp = 3;
    elseif (rfi_alg_ID==7),  rfi_weight_exp = 4;
    elseif (rfi_alg_ID==8),  rfi_weight_exp = 1; rfi_bin_limit = 0; rfi_bin_replace = 2.0;
    elseif (rfi_alg_ID==9),  rfi_weight_exp = 1; rfi_bin_limit = 12; rfi_bin_replace = 12;
    end
  elseif (1)
    rfi_weight_exp = 1.0;
    if (rfi_alg_ID==0)||(rfi_alg_ID==1),rfi_bin_limit = 13; rfi_bin_replace = 2.0;
    elseif (rfi_alg_ID==2), rfi_bin_limit = 12; rfi_bin_replace = 2.0;
    elseif (rfi_alg_ID==3), rfi_bin_limit = 11; rfi_bin_replace = 2.0;
    elseif (rfi_alg_ID==4), rfi_bin_limit = 10; rfi_bin_replace = 2.0;
    end
  elseif (1)
    rfi_weight_exp = 1.0;
    if (rfi_alg_ID==0)||(rfi_alg_ID==1),rfi_bin_limit = 12; rfi_bin_replace = 3.0;
    elseif (rfi_alg_ID==2), rfi_bin_limit = 12; rfi_bin_replace = 2.5;
    elseif (rfi_alg_ID==3), rfi_bin_limit = 12; rfi_bin_replace = 2.0;
    elseif (rfi_alg_ID==4), rfi_bin_limit = 12; rfi_bin_replace = 1.5;
    end
  end
  fprintf('INR=%.0f, N_avg_t=%.0f, N_avg_f=%.0f, duty=%.0f%%, alg=%.0f\n',...
          rfi_inr_db,N_rfi_avg_t,N_rfi_avg_f,rfi_duty_pct,rfi_alg_ID);
  fprintf('limit=%.0f, replace=%.2f, exp=%.2f, plot_ID=%.0f\n',...
          rfi_bin_limit,rfi_bin_replace,rfi_weight_exp,plot_ID);
  f1_rfi = [-750:100:500];
  df_dt_rfi = 0;
  mod_params_rfi.mod_type = 'BPSK';
  mod_params_rfi.f_sym = 4*fs0;
  mod_params_rfi.T_sym=1/mod_params_rfi.f_sym;
  
  snr_db_list = [1];  % statistical detection mode
  %snr_db_list = [-2:5];  % statistical detection mode
  m0_chirp_list = [10:17:1010];
  %m0_chirp_list = [10:5:300];
  N_tap_list = [4];
  Lf_list = [2];
  plot_multichirp_sg_rfi = 1;  % combined input and output SG
  plot_multichirp_sg = 0;   
  plot_multichirp_cat_sg1 = 1;
  plot_multichirp_cat_sg2 = 0;  % zoomed version
  plot_fastDD_det_plane = 0;
  plot_fastDD_det_profile = 0;
  plot_taylor_det_profile = 1;
  plot_taylor_det_plane = 1;
  pos_neg = 1;
  run_string = 'rfi_test';
elseif strcmp(run_id,'wide')
  % wide drift rate range run
  df_dt_min =  0;
  df_dt_max =  2*fs0.^2;
  m0_chirp_list = [10:17:2*Nt_pfb0];
  snr_db_list = [Inf];  % zero noise analysis mode
  N_tap_list = [1 3 8];
  Lf_list = [1 3/2 2 4];
%elseif ~isempty(strfind(run_id,'preDD-test'))
elseif contains(run_id,'preDD-test')
  N_pre_DD_nom = str2num(run_id(11:end));
  plot_ID = 10*N_pre_DD_nom;
  snr_db_list = Inf;  % statistical detection mode
  if (0)
    use_max_Lr = 0;  % =1 to set Lr=Lf, =0 use Lr_list
    Lr_list = [1 2 4];
  else
    use_max_Lr = 1;  % =1 to set Lr=Lf, =0 use Lr_list
  end
    N_tap_list = [3];
    Lf_list = [1 1.5 2 4];
  pos_neg = 1;
  run_string = '_preDD-test';
elseif strcmp(run_id,'C-BPSK1')
  % AM envelope characteristics
  snr_db_list = [Inf];  % zero noise analysis mode
  plot_ID = 1;
  Mx_list = [1 2];  % multiplier to initial M0 value
  N_tap_list = [1 3];
  Lf_list = [1.0 1.5 2];
  mod_type='C-BPSK1';   n_mod_ID = 10;
elseif strcmp(run_id,'C-BPSK2')
  % AM envelope characteristics
  snr_db_list = [Inf];  % zero noise analysis mode
  plot_ID = 1;
  Mx_list = [1 2];  % multiplier to initial M0 value
  N_tap_list = [1 3];
  Lf_list = [1.0 1.5 2];
  mod_type='C-BPSK2';   n_mod_ID = 10;
elseif strcmp(run_id,'C-QPSK')
  % AM envelope characteristics
  snr_db_list = [Inf];  % zero noise analysis mode
  plot_ID = 1;
  Mx_list = [1 2];  % multiplier to initial M0 value
  N_tap_list = [1 3];
  Lf_list = [1.0 1.5 2];
  mod_type='C-QPSK';   n_mod_ID = 10;
elseif strcmp(run_id,'BPSK')
  % BPSK envelope characteristics
  snr_db_list = [Inf];  % zero noise analysis mode
  Mx_list = [1 2];  % multiplier to initial M0 value
  plot_ID = 1;
  N_tap_list = [1 3];
  Lf_list = [1.0 1.5 2];
  mod_type='BPSK';   n_mod_ID = 10;
elseif strcmp(run_id,'QPSK')
  % BPSK envelope characteristics
  snr_db_list = [Inf];  % zero noise analysis mode
  Mx_list = [1 2];  % multiplier to initial M0 value
  plot_ID = 1;
  N_tap_list = [1 3];
  Lf_list = [1.0 1.5 2];
  mod_type='BPSK';   n_mod_ID = 10;
end


