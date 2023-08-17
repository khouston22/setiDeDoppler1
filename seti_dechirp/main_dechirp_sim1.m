%
% SETI De-Deoppler evaluation script
%

clear;

addpath('..\misc_fns');
addpath('..\pfb');
addpath('..\pfb_mex');
addpath('..\seti_dechirp_mex');

run_id_list = {'baseline-inf1a','baseline-inf1b','baseline-inf1c','baseline-inf1d',...
               'baseline-inf2a','baseline-inf2b',...
               'baseline-inf3','baseline-inf3b',...
               'baseline-inf4','baseline-inf5','baseline-inf6',...
               'baseline-snr1a','baseline-snr1b','baseline-snr2',...
               'alg-test1','alg-test2','alg-test3','alg-test4','alg-test5','alg-test6',...
               'N0-sweep32-Mx1','N0-sweep32-Mx2'};
run_id_list = {'baseline-inf1a','baseline-inf1b','baseline-inf1c','baseline-inf1d'};
%run_id_list = {'baseline-inf1b','baseline-inf1d'};
%run_id_list = {'baseline-inf2a','baseline-inf2b'};
%run_id_list = {'baseline-inf3a','baseline-inf3b'};
% runid_list = {'baseline-inf1a','baseline-inf1b','baseline-inf1c','baseline-inf1d',...
%                'baseline-inf2a','baseline-inf2b',...
%                'baseline-inf3','baseline-inf5'};
%run_id_list = {'baseline-snr1a','baseline-snr1b','baseline-snr2','N0-sweep32-Mx1','N0-sweep32-Mx2'};
%run_id_list = {'baseline-snr1a','baseline-snr1b','baseline-snr2'};
%run_id_list = {'baseline-inf6'};
%run_id_list = {'test'};
%run_id_list = {'alg-test1','alg-test2','alg-test3','alg-test4','alg-test5','alg-test6'};
%run_id_list = {'N0-sweep32-Mx1','N0-sweep32-Mx2'};
%run_id_list = {'wide'};
%run_id_list = {'BPSK','QPSK'};
%run_id_list = {'C-BPSK1','C-BPSK2','C-QPSK'};
%run_id_list = {'preDD-test4'};
%run_id_list = {'preDD-test1','preDD-test2','preDD-test4','preDD-test8','preDD-test16'};
%run_id_list = {'baseline-inf','alg-test1','alg-test2','alg-test3','alg-test4','alg-test5','alg-test6'};
% nominal without/with RFI proc
% run_id_list = {'rfi-00999000','rfi-00999001','rfi-10999000','rfi-10999001',...
%                'rfi-10999050','rfi-10999051','rfi-10999100''rfi-10999101',...
%                'rfi-10999500','rfi-10999501'};
%run_id_list = {'rfi-00999000','rfi-00999001','rfi-10999051','rfi-10999101','rfi-10999501'};
% time window size sweep
% run_id_list = {'rfi-00999000','rfi-00999001','rfi-00125001','rfi-00251001''rfi-00501001',...
%                'rfi-10999000','rfi-10999001','rfi-10125001','rfi-10251001''rfi-10501001',...
%                'rfi-10999050','rfi-10999051','rfi-10125051','rfi-10251051''rfi-10501051',...
%                'rfi-10999100','rfi-10999101','rfi-10125101','rfi-10251101''rfi-10501101',...
%                'rfi-10999500','rfi-10999501','rfi-10125501','rfi-10251501''rfi-10501501'};
% INR sweep
% run_id_list = {'rfi-00999000','rfi-00999001','rfi-00999002','rfi-00999003','rfi-00999004','rfi-00999005',...
%                'rfi-10999000','rfi-10999001','rfi-10999002','rfi-10999003','rfi-10999004','rfi-10999005',...
%                'rfi-15999000','rfi-15999001','rfi-15999002','rfi-15999003','rfi-15999004','rfi-15999005',...
%                'rfi-20999000','rfi-20999001','rfi-20999002','rfi-20999003','rfi-20999004','rfi-20999005'};
% run_id_list = {'rfi-10999000','rfi-10999001','rfi-10999050','rfi-10999051',...
%                'rfi-15999000','rfi-15999001','rfi-15999050','rfi-15999051',...
%                'rfi-20999000','rfi-20999001','rfi-20999050','rfi-20999051'};
% limit/replace options
% run_id_list = {'rfi-00999000','rfi-00999001','rfi-00999008','rfi-00999009',...
%                'rfi-10999000','rfi-10999001','rfi-10999008','rfi-10999009',...
%                'rfi-10999050','rfi-10999051','rfi-10999058','rfi-10999059',...
%                'rfi-10999100','rfi-10999101','rfi-10999108','rfi-10999109',...
%                'rfi-10999500','rfi-10999501','rfi-10999508','rfi-10999509'};

n_run = length(run_id_list);

for i_run = 1:n_run
  
  run_id = run_id_list{i_run};
  fprintf('\n\nRun %.0f of %.0f, Run ID = %s\n',i_run,n_run,run_id);
  
  %
  % set run parameters
  %

  set_sim_parameters1;

  if use_mex && (i_run==1)
    fprintf('Testing using mex functions\n');
    mex ../pfb_mex/pfb_Lfx_mex.c ../pfb_mex/pfb_fns.c ../pfb_mex/kiss_fft.c
    mex ../seti_dechirp_mex/taylor_DD_mex.c ../seti_dechirp_mex/taylor_DD_fns.c
    mex ../seti_dechirp_mex/fast_DD_mex.c ../seti_dechirp_mex/taylor_DD_fns.c ../seti_dechirp_mex/fast_DD_fns.c
  end
  
  %
  % set up for plotting
  %

  if ~exist('./plots-dechirp/','dir'), mkdir('./plots-dechirp/'); end
  if isempty(run_string)
    output_directory=sprintf('./plots-dechirp/dechirp');
  else
    output_directory=sprintf('./plots-dechirp/%s/',run_string);
  end
  if ~exist(output_directory,'dir'), mkdir(output_directory); end

  default_figure_position = [775 43 758 591];	% set figure size for png import to pptx
  figure(1);
  set(1,'Position',default_figure_position);
  commandwindow;

  png = print_fig;
  png.output_directory = output_directory;

  % create/empty log file
  % fp = fopen(fullfile(output_directory,'_log_taylor.txt'),'w');
  % fclose(fp);
  % fp = fopen(fullfile(output_directory,'_log_fastDD.txt'),'w');
  % fclose(fp);
  
  %
  % nested loops
  %

  for mod_ID=1:n_mod_ID
    mod_params = [];
    if isempty(mod_type)
      mod_str = '';
    else
      mod_params.mod_type = mod_type;
      mod_params.f_sym = (.25*(mod_ID-1)+.01)*fs0;
      mod_params.T_sym=1/mod_params.f_sym;
      if strcmp(mod_params.mod_type,'BPSK')
        mod_str = sprintf('BPSK, f-sym=%.2f, T-sym=%.2f',mod_params.f_sym,...
                           mod_params.T_sym);
        mod_str_fn = sprintf('BPSK-%.0f',mod_ID);
      else
        mod_params.Am=.25;
        mod_str = sprintf('%s %.2f, f-sym=%.2f, T-sym=%.2f',...
                           mod_params.mod_type,mod_params.Am,...
                           mod_params.f_sym,mod_params.T_sym);
        mod_str_fn = sprintf('%s-%.0f',mod_params.mod_type,mod_ID);
      end
      if (mod_ID==1)
        mod_str = 'No Modulation';
        mod_params.mod_type = 'NONE';
      end
    end
    
    %
    % determine amplitude and noise parameters
    %

    n_snr = length(snr_db_list);

    for i_snr=1:n_snr
      if isinf(snr_db_list(i_snr))
        sigma_list(i_snr) = 0; A_list(i_snr) = 1;
      else
        sigma_list(i_snr) = 1; A_list(i_snr) = 10.^(.05*snr_db_list(i_snr));
      end
    end

    %
    % generate parameters for multiplexed chirps
    %

    m0_chirp_list = sort(m0_chirp_list);
    n_m0_chirp = length(m0_chirp_list);
    n_m0_offset = length(m0_offset_list);
    n_f_offset = length(f_offset_list);
    n_chirp = n_m0_chirp*n_m0_offset*n_f_offset;

    f_incr_offset = 10*bin_bw;
    [f1_all,f2_all,df_dt_all,m0_chirp_all,f_offset_all,f_total] = ...
        gen_multi_chirp_parameters(...
             fs0,Nt_pfb0,M0,f_incr_offset,m0_chirp_list,m0_offset_list,f_offset_list);

    if any(m0_chirp_all<0)
      m0_chirp_all
      error(sprintf('Error in multichirp_sweep, m0_chirp_all must be non-negative'));
    end

    %
    % generate input waveform for multiplexed chirps
    %

    fprintf('SigGen %s %.0f       ',class_name,n_chirp);
    tstart = tic;
    sigma = 0; A = 1;
    s_in = gen_chirp_wf_multi(sigma,A,f1_all,df_dt_all,fs0,Nt_pfb0,M0,mod_params,class_name);
    if ~isempty(f1_rfi)
      fprintf('RFI %.0f signals ',length(f1_rfi));
      s_in = s_in + ...
        gen_chirp_wf_multi(0,A_rfi,f1_rfi,df_dt_rfi,fs0,Nt_pfb0,M0,mod_params_rfi,class_name);
    end
    fprintf(' %.1f sec\n',toc(tstart));
    
    %
    % set up list of N_tap, Lf and Lr values to process
    % note: i_case=1 is baseline
    %
    
    N_tap_case  = N_tap_baseline;   
    Lf_case = Lf_baseline; 
    Lr_case = Lr_baseline; 
    wind_case = {'RectWin'};
    M_case = M0;
    fs_case = fs_in/M0;
    i_case = 1;
    
    for Mx = Mx_list
      for i_Lf = 1:length(Lf_list)
        Lf = Lf_list(i_Lf);
        for i_N_tap = 1:length(N_tap_list)
          if (N_tap_list(i_N_tap)==1)
            w_list = window1_list;  % 'RectWin' & others
          else
            w_list = window2_list;
          end
          for i_w_list = 1:length(w_list)
            if use_max_Lr
              i_case = i_case+1;
              N_tap_case(i_case)  = N_tap_list(i_N_tap);
              Lf_case(i_case) = Lf;
              Lr_case(i_case) = Lf;
              M_case(i_case) = M0*Mx;
              fs_case(i_case) = fs_in/(M0*Mx);
              wind_case(i_case) = w_list(i_w_list);
            elseif (length(Lr_list)==length(Lf_list))
              i_case = i_case+1;
              N_tap_case(i_case)  = N_tap_list(i_N_tap);
              Lf_case(i_case) = Lf;
              Lr_case(i_case) = Lr_list(i_Lf);
              M_case(i_case) = M0*Mx;
              fs_case(i_case) = fs_in/(M0*Mx);
              wind_case(i_case) = w_list(i_w_list);
            else
              for Lr = Lr_list
                if (Lr<=Lf) && (Lr>=Lf/2)
                  i_case = i_case+1;
                  N_tap_case(i_case)  = N_tap_list(i_N_tap);
                  Lf_case(i_case) = Lf;
                  Lr_case(i_case) = Lr;
                  M_case(i_case) = M0*Mx;
                  fs_case(i_case) = fs_in/(M0*Mx);
                  wind_case(i_case) = w_list(i_w_list);
                end
              end
            end
          end
        end
      end
    end
    n_case = i_case;
    fprintf('M:    '); fprintf('%4.0f ',M_case); fprintf('\n');
    fprintf('fs:   '); fprintf('%4.1f ',fs_case); fprintf('\n');
    fprintf('N_tap:'); fprintf('%4.0f ',N_tap_case); fprintf('\n');
    fprintf('Lf:   '); fprintf('%4.1f ',Lf_case); fprintf('\n');
    fprintf('Lr:   '); fprintf('%4.1f ',Lr_case); fprintf('\n');
    fprintf('Wind: ');
    for i_case=1:n_case, fprintf('%s ',wind_case{i_case}(1:4)); end;
    fprintf('\n');

    n_N0 = length(N0_list);

    for i_N0 = 1:n_N0
      
      N0 = N0_list(i_N0);

      i_case = 0;
      run_legend_str = {};
      
      SG_size = NaN(n_case);
      DP_size  = NaN(n_case);
      SG_ratio = NaN(n_case);
      DP_ratio = NaN(n_case);

      for i_case = 1:n_case

        N_tap = N_tap_case(i_case);
        Lf = Lf_case(i_case);
        Lr = Lr_case(i_case);
        M = M_case(i_case);
        Nt_pfb = Nt_pfb0*M0/M;
        fs = fs_in/M;
        Ts = 1/fs;
        T_avg = Nt_pfb0/fs0;
        window_name = wind_case{i_case};
        if strcmp(window_name,'lpfxxx')
          window_name = choose_lpf_window(Lf,N_tap);
        end
        
        if (i_case==1)    % baseline
          N0_taylor = N0_baseline;
        elseif isempty(N0_taylor0)
          N0_taylor = N0;
        else
          N0_taylor = N0_taylor0;
        end
        
        if (i_case==1)    % baseline
          N_pre_DD = N_pre_DD_baseline;
        else
          N_pre_DD = N_pre_DD_nom;
        end
        N_DD = Nt_pfb/N_pre_DD;
        T_line = N_pre_DD*Ts;

        freq = [-Lf*M/2:Lf*M/2-1]/(Lf*M)*fs_in;
        i_freq0 = Lf*M/2+1;

        if use_max_Lr
          run_legend_str{i_case} = sprintf('Ntap=%.0f %s %.2fx',N_tap,window_name,Lf);
        else
          run_legend_str{i_case} = sprintf('Ntap=%.0f %s %.1fx Lr=%.2f',N_tap,window_name,Lf,Lr);
        end
        
        if (N_pre_DD>1)
          run_legend_str{i_case} = [run_legend_str{i_case} sprintf(' NpreDD=%.0f',N_pre_DD)];
        end
        
        %if (length(Mx_list)>1)
        %  run_legend_str{i_case} = [sprintf('M=%.0f ',M) run_legend_str{i_case}];
        %end
        
        fprintf('\nRun %.0f of %.0f, Run ID = %s ',i_run,n_run,run_id);
        fprintf('Case %.0f of %.0f, %.0f Chirps, %s\n',...
          i_case,n_case,n_chirp,run_legend_str{i_case});
        
        n_freq = length(freq);
        df_bin = freq(2)-freq(1);

        if (length(Mx_list)>1)
          config_str0 = sprintf('fs=%.1f-%.1f, Tavg=%.0f',min(fs_case),max(fs_case),T_avg);
          config_str0_fn = sprintf('fs-%.1f-%.1f-Tavg-%.0f',min(fs_case),max(fs_case),T_avg);
          run_legend_str{i_case} = [run_legend_str{i_case} sprintf(' fs=%.1f',fs)];
        else
          config_str0 = sprintf('fs=%.1f, Tavg=%.0f',fs,T_avg);
          config_str0_fn = sprintf('fs-%.1f-Tavg-%.0f',fs,T_avg); 
        end
        config_str1 = sprintf('Ntap=%.0f %s %.2fx Lr=%.2f, fs=%.1f, N=%.0f',...
          N_tap,window_name,Lf,Lr,fs,Nt_pfb);
        config_str2 = sprintf('Ntap=%.0f %s %.2fx Lr=%.2f, fs=%.1f Hz, N=%.0f, N0=%.0f',...
          N_tap,window_name,Lf,Lr,fs,Nt_pfb,N0);

        config_str1_fn = sprintf('%.2fx-Lr-%.1f-M-%.0f-Ntap-%.0f-%s-fs-%.1f-N-%.0f-%.0f',...
          Lf,Lr,M,N_tap,window_name,fs,Nt_pfb,N_pre_DD);
        config_str2_fn = sprintf('%.2fx-Lr-%.1f-M-%.0f-Ntap-%.0f-%s-fs-%.1f-N-%.0f-N0-%.0f-%.0f',...
          Lf,Lr,M,N_tap,window_name,fs,Nt_pfb,N0,N_pre_DD);

        if (N_pre_DD>1)
          config_str0 = [config_str0 sprintf(' NpreDD=%.0f',N_pre_DD)];
          config_str1 = [config_str1 sprintf(' NpreDD=%.0f',N_pre_DD)];
          config_str2 = [config_str2 sprintf(' NpreDD=%.0f',N_pre_DD)];
        end
        
        if (i_case==1)
          init_run_arrays;
        end

        %
        % run pfb with noise only
        %

        fprintf('PFB %.2fx Noise      ',Lf);

        tstart = tic;
        sigma = sqrt(M0); A = 0;
        [noise_out,freq1,t_out,coef] = gen_chirp_pfb_matrix(...
                  sigma,A,[],[],fs,Nt_pfb,M,N_tap,window_name,Lf,mod_params,class_name);
        t_cpu_pfb(1:n_chirp,i_case) = toc(tstart);
        fprintf(' %.1f sec\n',t_cpu_pfb(1,i_case));
        
        if (N_pre_DD==1)
          xx = abs(noise_out).^2;
        else
          [xx,t_out] = pre_DD_avg(abs(noise_out).^2,t_out,N_pre_DD);
        end

        bin_noise_magsq0 = mean(xx(:));

        if rfi_proc_enable
          [xx,bin_power_in,bin_power_out] = rfi_agc1(xx,N_rfi_avg_f,N_rfi_avg_t,rfi_weight_exp,...
                      rfi_bin_limit,rfi_bin_replace,bin_noise_magsq0);
        end
      
        %bin_noise_magsq = rms(noise_out(:)).^2;
        bin_noise_magsq = mean(xx(:));
        bin_noise_db0 = 10*log10(bin_noise_magsq);

        % run noise through taylor & determine detection threshold

        fprintf('Taylor-%2.0f v%.0f Noise  ',N0_taylor,ver_DD);
        tstart = tic;
        if (ver_DD==4)
          [det_noise_taylor] = taylor_DD_mex(xx,freq,T_line,Lf,Lr,N_DD,N0_taylor,df_dt_min,df_dt_max);
        elseif (ver_DD==3)
          [det_noise_taylor] = taylorDD3(xx,freq,T_line,Lf,Lr,N_DD,N0_taylor,df_dt_min,df_dt_max);
        end
        fprintf(' %.1f sec\n',toc(tstart));
        
        DT_taylor = mean(det_noise_taylor(:))+DT_std_dev*std(det_noise_taylor(:));
        mdl_taylor = DT_taylor - bin_noise_magsq;
        DT_taylor_db = 10*log10(DT_taylor);
        delta_DT_taylor_db = 0;

        % run noise through fastDD & determine detection threshold
        
        fprintf('fastDD%.0f-%2.0f v%.0f Noise  ',fastDD_alg_ID,N0,ver_DD);
        tstart = tic;
        if (ver_DD==4)
          [det_noise_fastDD] = fast_DD_mex(xx,freq,T_line,Lf,Lr,N_DD,N0,df_dt_min,df_dt_max,fastDD_alg_ID);
        elseif (ver_DD==3)
          [det_noise_fastDD] = fastDD3(xx,freq,T_line,Lf,Lr,N_DD,N0,df_dt_min,df_dt_max,fastDD_alg_ID);
        end
        fprintf(' %.1f sec\n',toc(tstart));        
        
        DT_fastDD = mean(det_noise_fastDD(:))+DT_std_dev*std(det_noise_fastDD(:));
        mdl_fastDD = DT_fastDD - bin_noise_magsq;
        DT_fastDD_db = 10*log10(DT_fastDD);
        %delta_DT_fastDD_db = DT_fastDD_db - DT_taylor_db;
        delta_DT_fastDD_db = 10*log10(mdl_fastDD/mdl_taylor);
        
        SG_size(i_case) = prod(size(xx));
        DP_size(i_case) = prod(size(det_noise_fastDD));
        SG_ratio(i_case) = SG_size(i_case)/SG_size(1);
        DP_ratio(i_case) = DP_size(i_case)/DP_size(1);
        
        if (1)
        fprintf('Case %.0f fs %.1f Lf %.2f: SG size %.0f x %.0f (%.2fx), fastDD DP size %.0f x %.0f (%.2fx)\n',...
          i_case,fs,Lf,size(xx),SG_ratio(i_case),size(det_noise_fastDD),DP_ratio(i_case));
        end

%         fprintf('\n*** %s\n',config_str1);
%         fprintf('DT_taylor = %.2f, bin_noise_magsq = %.2f, mdl_taylor = %.2f, DT_taylor_db = %.2f\n',...
%             DT_taylor,bin_noise_magsq,mdl_taylor,DT_taylor_db);
%         fprintf('mean, std = %.2f %.2f\n',mean(det_noise_taylor(:)),std(det_noise_taylor(:)));
%         fprintf('DT_fastDD = %.2f, bin_noise_magsq = %.2f, mdl_fastDD = %.2f, DT_fastDD_db = %.2f\n',...
%             DT_fastDD,bin_noise_magsq,mdl_fastDD,DT_fastDD_db);
%         fprintf('mean, std = %.2f %.2f, delta_DT_fastDD_db = %.2f\n',mean(det_noise_fastDD(:)),std(det_noise_fastDD(:)),delta_DT_fastDD_db);
%         fprintf('***\n\n');

        %
        % sweep drift rates and offsets
        %

        multichirp_sweep;

        if (n_chirp==1)
          continue;
        end

        %
        % compute means for each case
        %

        for i_snr=1:n_snr

          snr_db = snr_db_list(i_snr);

          if (i_snr==n_snr)
            if use_weighted_mean
              mean_best_db(i_case) = GWmean(taylor_best_db(:,i_snr,i_case),df_dt_all);
            else
              mean_best_db(i_case) = mean(taylor_best_db(:,i_snr,i_case));
            end
            mean_rel_snr_best_db(i_case) = mean_best_db(i_case) - bin_noise_db0;
          end

          if (~isinf(snr_db) & (i_snr==n_snr))
            if (n_snr>1)
              det_snr_taylor_db(i_N0,i_case) = interp1(mean_SE_taylor_db(:,i_case),snr_db_list,0,'linear','extrap');
              det_snr_fastDD_db(i_N0,i_case) = interp1(mean_SE_fastDD_db(:,i_case),snr_db_list,0,'linear','extrap');
            else
              det_snr_taylor_db(i_N0,i_case)  = snr_db_list(i_snr)-mean_SE_taylor_db(i_snr,i_case);
              det_snr_fastDD_db(i_N0,i_case) = snr_db_list(i_snr)-mean_SE_fastDD_db(i_snr,i_case);
            end
          end

          if (i_snr==n_snr)
            mean_t_cpu_taylor(i_N0,i_case)  = mean(mean(t_cpu_taylor(:,:,i_case)));
            mean_t_cpu_fastDD(i_N0,i_case) = mean(mean(t_cpu_fastDD(:,:,i_case)));
            mean_t_cpu_taylor_total(i_N0,i_case)  = mean_t_cpu_pfb(i_N0,i_case) + mean_t_cpu_taylor(i_N0,i_case);
            mean_t_cpu_fastDD_total(i_N0,i_case) = mean_t_cpu_pfb(i_N0,i_case) + mean_t_cpu_fastDD(i_N0,i_case);
            t_cpu_baseline = mean_t_cpu_taylor_total(i_N0,1);
            t_cpu_str = sprintf('CPU: Taylor %.2fx, fastDD %.2fx Baseline',...
                                mean_t_cpu_taylor_total(i_N0,i_case)/t_cpu_baseline,...
                                mean_t_cpu_fastDD_total(i_N0,i_case)/t_cpu_baseline);
          end
          
          %
          % store baseline separately
          %
          
          if (i_case==1)
            if strcmp(lower(alg_baseline),'taylor')
              baseline_db = taylor_db(:,:,1);
              baseline_SE_db = taylor_SE_db(:,:,1);
              rel_snr_baseline_db = rel_snr_taylor_db(:,1);
              SE_baseline_db = SE_taylor_db(:,:,1);
              mean_baseline_db = mean_taylor_db(:,1);
              mean_rel_snr_baseline_db = mean_rel_snr_taylor_db(i_N0,1);
              mean_SE_baseline_db = mean_SE_taylor_db(:,1);
              det_snr_baseline_db = det_snr_taylor_db(1);
            elseif strcmp(lower(alg_baseline),'fastDD')
              baseline_db = fastDD_db(:,:,1);
              baseline_SE_db = fastDD_SE_db(:,:,1);
              rel_snr_baseline_db = rel_snr_fastDD_db(:,1);
              SE_baseline_db = SE_fastDD_db(:,:,1);
              mean_baseline_db = mean_fastDD_db(:,1);
              mean_rel_snr_baseline_db = mean_rel_snr_fastDD_db(i_N0,1);
              mean_SE_baseline_db = mean_SE_fastDD_db(:,1);
              det_snr_baseline_db = det_snr_fastDD_db(1);
            end
          end

          % SNR improvement over baseline
          if (isinf(snr_db) & (i_snr==n_snr))
            delta_snr_taylor_db(i_N0,i_case) = ...
                     mean_rel_snr_taylor_db(i_N0,i_case) - mean_rel_snr_baseline_db;
            delta_snr_fastDD_db(i_N0,i_case) = ...
                     mean_rel_snr_fastDD_db(i_N0,i_case) - mean_rel_snr_baseline_db;
            value_metric_taylor(i_N0,i_case) = delta_snr_taylor_db(i_N0,i_case)...
                       /mean_t_cpu_taylor_total(i_N0,i_case);
            value_metric_fastDD(i_N0,i_case) = delta_snr_fastDD_db(i_N0,i_case)...
                       /mean_t_cpu_fastDD_total(i_N0,i_case);
          end
          
          if (i_snr==n_snr)
            delta_det_snr_taylor_db(i_N0,i_case) = ...
                     det_snr_baseline_db - det_snr_taylor_db(i_N0,i_case);
            delta_det_snr_fastDD_db(i_N0,i_case) = ...
                     det_snr_baseline_db - det_snr_fastDD_db(i_N0,i_case);
          end
          
          %   mean_scp_db0 = mean(scp_db(:,:,i_case));
          %   mean_scp_db1 = mean_scp_db0 - mean_scp_db0(1);
          %   mean_scp_db(:,i_case) = mean_scp_db0';
          %   fprintf('Mean SCP dB     = '); fprintf('%6.2f ',mean_scp_db0); fprintf('\n');
          %   fprintf('Mean wrt n-diff=1 = '); fprintf('%6.2f ',mean_scp_db1); fprintf('\n');

          
          best_str    = sprintf('Best:        relSNR=%6.2f dB',...
                                mean_rel_snr_best_db(i_case));
          taylor_str0 = sprintf('Taylor-%.0f',N0_taylor);
          fastDD_str0 = sprintf('fastDD%.0f-%.0f',fastDD_alg_ID,N0);
          baseline_str0 = sprintf('Baseline Taylor-%.0f',N0_baseline);
                              
          taylor_str1 = sprintf('Taylor: relSNR=%6.2f dB, %6.2f dB Over Baseline',...
                                mean_rel_snr_taylor_db(i_N0,i_case),delta_snr_taylor_db(i_N0,i_case));
          fastDD_str1 = sprintf('fastDD: relSNR=%6.2f dB, %6.2f dB Over Baseline',...
                                mean_rel_snr_fastDD_db(i_N0,i_case),delta_snr_fastDD_db(i_N0,i_case));

          taylor_str2 = sprintf('Taylor: Mean Signal Excess=%6.2f dB',...
                                mean_SE_taylor_db(i_snr,i_case));
          fastDD_str2 = sprintf('fastDD: Mean Signal Excess=%6.2f dB',...
                                mean_SE_fastDD_db(i_snr,i_case));
                              
          baseline_str1  = sprintf(...
            'Baseline:   relSNR=%6.2f dB  (%s-%.0f Ntap=%.0f %.2fx)',...
                              mean_rel_snr_baseline_db,...
                              alg_baseline, N0_baseline,N_tap_baseline,Lf_baseline);
          baseline_str2  = sprintf(...
            'Baseline:  Mean Signal Excess=%6.2f dB (%s-%.0f Ntap=%.0f %.2fx)',...
                              mean_SE_baseline_db(i_snr),...
                              alg_baseline, N0_baseline,N_tap_baseline,Lf_baseline);

          noise_str  = ['Bin Noise   ' sprintf('%.2f ',bin_noise_db0) 'dB'];
          DT_str  = sprintf('DT  Taylor %.2f  FastDD %.2f delta %.2f dB',...
                            DT_taylor_db,DT_fastDD_db,delta_DT_fastDD_db);

          if (n_snr==1)
            fprintf('%s\n',best_str);
            fprintf('%s\n',taylor_str1);
            fprintf('%s\n',fastDD_str1);
            fprintf('%s\n',t_cpu_str);
          end

          if (isinf(snr_db) & (n_snr==1))
            % SNR=Inf analysis
            fp = fopen(fullfile(output_directory,'_log_taylor.txt'),'a');
            fprintf(fp,...
              '%.2fx, Lr=%.2f, fs=%.2f\t%s Ntap=%.0f\tTaylor %2.0f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',...
                    Lf,Lr,fs,window_name,N_tap,N0_taylor,mean_taylor_db(i_snr,i_case),-bin_noise_db0,...
                    -delta_DT_taylor_db,mean_rel_snr_taylor_db(i_N0,i_case),...
                    delta_snr_taylor_db(i_N0,i_case),0,...
                    mean_t_cpu_pfb(i_N0,i_case),mean_t_cpu_taylor(i_N0,i_case),...
                    mean_t_cpu_taylor_total(i_N0,i_case),value_metric_taylor(i_N0,i_case));
            fclose(fp);
            fp = fopen(fullfile(output_directory,'_log_fastDD.txt'),'a');
            fprintf(fp,...
              '%.2fx, Lr=%.2f, fs=%.2f\t%s Ntap=%.0f\tfastDD %2.0f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',...
                    Lf,Lr,fs,window_name,N_tap,N0,mean_fastDD_db(i_snr,i_case),-bin_noise_db0,...
                    -delta_DT_fastDD_db,mean_rel_snr_fastDD_db(i_N0,i_case),...
                    delta_snr_fastDD_db(i_N0,i_case),0,...
                    mean_t_cpu_pfb(i_case),mean_t_cpu_fastDD(i_case),...
                    mean_t_cpu_fastDD_total(i_case),value_metric_fastDD(i_N0,i_case));
            fclose(fp);
          elseif (~isinf(snr_db) & (i_snr==n_snr))
            % SNR sweep
            fp = fopen(fullfile(output_directory,'_log_taylor.txt'),'a');
            fprintf(fp,...
              '%.2fx, Lr=%.2f, fs=%.2f\t%s Ntap=%.0f\tTaylor %2.0f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',...
                    Lf,Lr,fs,window_name,N_tap,N0_taylor,det_snr_taylor_db(i_N0,i_case),...
                    delta_det_snr_taylor_db(i_N0,i_case),...
                    mean_t_cpu_pfb(i_N0,i_case),mean_t_cpu_taylor(i_N0,i_case),...
                    mean_t_cpu_taylor_total(i_N0,i_case),value_metric_taylor(i_N0,i_case));
            fclose(fp);
            fp = fopen(fullfile(output_directory,'_log_fastDD.txt'),'a');
            fprintf(fp,...
              '%.2fx, Lr=%.2f, fs=%.2f\t%s Ntap=%.0f\tfastDD %2.0f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',...
                    Lf,Lr,fs,window_name,N_tap,N0,det_snr_fastDD_db(i_N0,i_case),...
                    delta_det_snr_fastDD_db(i_N0,i_case),...
                    mean_t_cpu_pfb(i_N0,i_case),mean_t_cpu_fastDD(i_N0,i_case),...
                    mean_t_cpu_fastDD_total(i_N0,i_case),value_metric_taylor(i_N0,i_case));
            fclose(fp);
          end

          %
          % plot outputs
          % 

          if (i_case>1)
            plot_run_outputs;
          end

        end % snr
      end  % i_case

      if (0)    
        taylor_best_db, fastDD_db, taylor_db
        mean_best_db, mean_taylor_db, mean_fastDD_db
      end

      %
      % plot comparison of ideal deDoppler curves
      % 

      plot_meta_run;

    end % N0
  end % mod_ID

end % run_id_list


