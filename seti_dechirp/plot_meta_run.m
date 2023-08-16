%
% plot comparison of ideal deDoppler curves
%

if (1 && isinf(snr_db))
  if isempty(mod_str)
    png.file_name = sprintf('00-DD-best-all-db-%s',config_str0_fn);
  else
    png.file_name = sprintf('00-DD-best-all-db-%s-%s',config_str0_fn,mod_str_fn);
  end
  png.file_count = plot_ID-1;
  figure(1); clf;
  n_Lf = length(Lf_list);
  n_N_tap = length(N_tap_list);
  legend_str = [];
  for i_case = n_case:-1:2
    Lf1 = Lf_case(i_case);
    N_tap1 = N_tap_case(i_case);
    if (N_tap1==1), ctl='-*'; elseif (N_tap1==4), ctl='-^'; else ctl='-o'; end
    plot(df_dt_all,rel_snr_best_db(:,i_case),ctl); hold on;
    legend_str = [legend_str run_legend_str(i_case)];
    if (Lf1==1), set(gca,'ColorOrderIndex',1); end
  end
  v=axis; axis([v(1:2) -3 4]); v=axis;
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  xlabel('Drift Rate (Hz/sec)');
  ylabel('Relative SNR dB')
  title(sprintf('Relative SNR for DeDoppler with Perfect Indexing: %s',config_str0));
  legend(legend_str,'Location','NorthEast')
  if ~isempty(mod_str)
    text_sc(.05,.11,mod_str);
  end
  grid;
  print_fig(png);
end

if (1 && isinf(snr_db))

  n_N_tap = length(N_tap_list);

  for i=1:n_N_tap
    N_tap1 = N_tap_list(i);
    if isempty(mod_str)
      png.file_name = sprintf('01-DD-best-db-N-tap-%.0f-%s',N_tap1,config_str0_fn);
    else
      png.file_name = sprintf('01-DD-best-db-N-tap-%.0f-%s-%s',N_tap1,config_str0_fn,mod_str_fn);
    end
    png.file_count = plot_ID-1;
    figure(1); clf;
    legend_str = [];
    ii = find(N_tap_case(2:end)==N_tap1) + 1;
    plot(df_dt_all,squeeze(rel_snr_best_db(:,ii)),'-*'); hold on;
    if (max(m0_chirp_all)>Nt_pfb)
      plot([1 1],[-10 10],'--k','LineWidth',1);
    end
    v=axis; axis([v(1:2) -3 4]); v=axis;
    set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
    xlabel('Drift Rate (Hz/sec)');
    ylabel('Relative SNR dB')
    str = run_legend_str{ii(1)};
    title(sprintf('Relative SNR for DeDoppler with Perfect Indexing: %s',config_str0));
    if (max(m0_chirp_all)>Nt_pfb)
      legend([run_legend_str(ii) {'1 bin/time sample @1x'}],'Location','NorthEast')
    else
      legend(run_legend_str(ii),'Location','NorthEast');
    end
    text_sc(.05,.05,config_str0);
    text_sc(.05,.08,['Bin Noise = ' sprintf('%.2f ',bin_noise_db(1,1,ii(1))) 'dB']);
    if ~isempty(mod_str)
      text_sc(.05,.11,mod_str);
    end
    grid;
    png = print_fig(png);
  end
end


if (1 && isinf(snr_db) && n_N0>1 && i_N0==n_N0)
  for i_plt = 1:2
    if i_plt==1
      png.file_name = sprintf('02a-fastDD-dsnr-vs-comp-%s',config_str0_fn);
      mean_t_cpu_total = mean_t_cpu_fastDD_total;
      delta_snr_db = delta_snr_fastDD_db;
      alg_name = sprintf('fastDD%.0f',fastDD_alg_ID);
    else
      png.file_name = sprintf('02b-taylor-dsnr-vs-comp-%s',config_str0_fn);
      mean_t_cpu_total = mean_t_cpu_taylor_total;
      delta_snr_db = delta_snr_taylor_db;
      alg_name = 'Taylor';
   end
    t_cpu_baseline = mean(mean_t_cpu_taylor_total(:,1));
    if ~isempty(mod_str)
      png.file_name = [png.file_name sprintf('-%s',mod_str_fn)];
    end
    png.file_count = plot_ID-1;
    figure(1); clf;
    n_Lf = length(Lf_list);
    n_N_tap = length(N_tap_list);
    legend_str = [];
    for i_case = 2:n_case
      Lf1 = Lf_case(i_case);
      N_tap1 = N_tap_case(i_case);
      ctl = set_symbol(i_case);
      semilogx(mean_t_cpu_total(:,i_case)/t_cpu_baseline,...
               delta_snr_db(:,i_case),ctl,'LineWidth',2); hold on;
      legend_str = [legend_str run_legend_str(i_case)];
      %if (Lf1==1), set(gca,'ColorOrderIndex',1); end
    end
    semilogx(1,0,'^k','LineWidth',2); hold on;
    semilogx([1 21.5],[0 4],':k','LineWidth',2); hold on;
    semilogx([1 6.35],[0 4],':k','LineWidth',2); hold on;
    legend_str = [legend_str ['Baseline: ' run_legend_str{1}]];
    v=axis; axis([1 500 0 4]); v=axis;
    set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
    xlabel('Relative Compute Time')
    ylabel('SNR Over Baseline dB')
    title(sprintf('%s SNR Gain vs Compute Cost: %s, N0=%.0f-%.0f',alg_name,...
              config_str0,min(N0_list),max(N0_list)));
    legend(legend_str,'Location','SouthEast')
    text(4.0,.10,['N0=' sprintf('%.0f ',N0_list)],'FontSize',12);
    text(1.4,3.6,'Reference Line 1','FontSize',12);
    text(18,3.6,'Reference Line 2','FontSize',12);
    if ~isempty(mod_str)
      text(4.0,.25,mod_str,'FontSize',12);
    end
    grid;
    print_fig(png);
  end
end


if (1 && isinf(snr_db))
  if isempty(mod_str)
    png.file_name = sprintf('03a-fastDD-dsnr-vs-comp-%s-N0-%.0f',config_str0_fn,N0);
  else
    png.file_name = sprintf('03a-fastDD-dsnr-vs-comp-%s-N0-%.0f-%s',config_str0_fn,N0,mod_str_fn);
  end
  png.file_count = plot_ID-1;
  figure(1); clf;
  t_cpu_baseline = mean_t_cpu_taylor_total(i_N0,1);
  n_Lf = length(Lf_list);
  n_N_tap = length(N_tap_list);
  legend_str = [];
  for i_case = 2:n_case
    Lf1 = Lf_case(i_case);
    N_tap1 = N_tap_case(i_case);
    ctl = set_symbol(i_case);
    semilogx(mean_t_cpu_fastDD_total(i_N0,i_case)/t_cpu_baseline,...
             delta_snr_fastDD_db(i_N0,i_case),ctl,'LineWidth',2); hold on;
    legend_str = [legend_str run_legend_str(i_case)];
    %if (Lf1==1), set(gca,'ColorOrderIndex',1); end
  end
  semilogx(1,0,'^k','LineWidth',2); hold on;
  semilogx([1 21.5],[0 4],':k','LineWidth',2); hold on;
  semilogx([1 6.35],[0 4],':k','LineWidth',2); hold on;
  legend_str = [legend_str ['Baseline: ' run_legend_str{1}]];
  v=axis; axis([1 200 0 4]); v=axis;
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  xlabel('Relative Compute Time')
  ylabel('SNR Over Baseline dB')
  title(sprintf('fastDD%.0f SNR Gain vs Compute Cost: %s',fastDD_alg_ID,config_str0));
  legend(legend_str,'Location','SouthEast')
  text(4.0,.10,['N0=' sprintf('%.0f ',N0_list)],'FontSize',12);
  text(1.4,3.6,'Reference Line 1','FontSize',12);
  text(18,3.6,'Reference Line 2','FontSize',12);
  if ~isempty(mod_str)
    text(4.0,.25,mod_str,'FontSize',12);
  end
  grid;
  print_fig(png);
end

if (1 && isinf(snr_db))
  if isempty(mod_str)
    png.file_name = sprintf('03b-taylor-dsnr-vs-comp-%s-N0-%.0f',config_str0_fn,N0);
  else
    png.file_name = sprintf('03b-taylor-dsnr-vs-comp-%s-N0-%.0f-%s',config_str0_fn,N0,mod_str_fn);
  end
  png.file_count = plot_ID-1;
  figure(1); clf;
  t_cpu_baseline = mean_t_cpu_taylor_total(i_N0,1);
  n_Lf = length(Lf_list);
  n_N_tap = length(N_tap_list);
  legend_str = [];
  for i_case = 2:n_case
    Lf1 = Lf_case(i_case);
    N_tap1 = N_tap_case(i_case);
    ctl = set_symbol(i_case);
    semilogx(mean_t_cpu_taylor_total(i_N0,i_case)/t_cpu_baseline,...
             delta_snr_taylor_db(i_N0,i_case),ctl,'LineWidth',2); hold on;
    legend_str = [legend_str run_legend_str(i_case)];
    %if (Lf1==1), set(gca,'ColorOrderIndex',1); end
  end
  semilogx(1,0,'^k','LineWidth',2); hold on;
  semilogx([1 21.5],[0 4],':k','LineWidth',2); hold on;
  semilogx([1 6.35],[0 4],':k','LineWidth',2); hold on;
  legend_str = [legend_str ['Baseline: ' run_legend_str{1}]];
  v=axis; axis([1 200 0 4]); v=axis;
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  xlabel('Relative Compute Time')
  ylabel('SNR Over Baseline dB')
  title(sprintf('Taylor SNR Gain vs Compute Cost: %s',config_str0));
  legend(legend_str,'Location','SouthEast')
  if isempty(N0_taylor0)
    text(4.0,.10,['N0=' sprintf('%.0f ',N0_list)],'FontSize',12);
  else
    text(4.0,.10,['N0=' sprintf('%.0f ',N0_taylor0)],'FontSize',12);
  end
  text(1.4,3.6,'Reference Line 1','FontSize',12);
  text(18,3.6,'Reference Line 2','FontSize',12);
  if ~isempty(mod_str)
    text(4.0,.25,mod_str,'FontSize',12);
  end
  grid;
  print_fig(png);
end

if (1 && isinf(snr_db) && (n_case>2))
  if isempty(mod_str)
    png.file_name = sprintf('04a-taylor-fastDD-dSNR-db-%s-N0-%.0f',config_str0_fn,N0);
  else
    png.file_name = sprintf('04a-taylor-fastDD-dSNR-db-%s-N0-%.0f-%s',config_str0_fn,N0,mod_str_fn);
  end
  png.file_count = plot_ID-1;
  figure(1); clf;
  bar([1:n_case-1]',[delta_snr_taylor_db(i_N0,2:n_case) ; delta_snr_fastDD_db(i_N0,2:n_case)]');
  xlim([0 n_case]);
  set(gca,'XTick',[1:n_case-1]);
  set(gca,'XTickLabel',run_legend_str(2:n_case));
  ylim([0 5]);
  xtickangle(45);
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  %xlabel('Normalized SNR dB')
  ylabel('SNR Over Baseline dB')
  title(sprintf('SNR Gain vs PFB Configuration: %s, N0=%.0f',config_str0,N0));
  legend(sprintf('Taylor-%.0f',N0_taylor),sprintf('fastDD%.0f-%.0f',fastDD_alg_ID,N0),...
         'Location','NorthEast');
  text_sc(.05,.95,sprintf('Baseline: %s: %s-%.0f',...
          run_legend_str{1},alg_baseline,N0_baseline),11);
  grid;
  print_fig(png);
end

if (1 && isinf(snr_db) && (n_case>2))
  if isempty(mod_str)
    png.file_name = sprintf('04b-taylor-fastDD-dSNR-db-%s-N0-%.0f',config_str0_fn,N0);
  else
    png.file_name = sprintf('04b-taylor-fastDD-dSNR-db-%s-N0-%.0f-%s',config_str0_fn,N0,mod_str_fn);
  end
  png.file_count = plot_ID-1;
  figure(1); clf;
  plot([1:n_case-1]',[delta_snr_taylor_db(i_N0,2:n_case) ; delta_snr_fastDD_db(i_N0,2:n_case)]',...
       '-o','LineWidth',2);
  xlim([0 n_case]);
  set(gca,'XTick',[1:n_case-1]);
  set(gca,'XTickLabel',run_legend_str(2:n_case));
  ylim([0 5]);
  xtickangle(45);
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  %xlabel('Normalized SNR dB')
  ylabel('SNR Over Baseline dB')
  title(sprintf('SNR Gain vs PFB Configuration: %s, N0=%.0f',config_str0,N0));
  legend(sprintf('Taylor-%.0f',N0_taylor),sprintf('fastDD%.0f-%.0f',fastDD_alg_ID,N0),...
         'Location','NorthEast');
  text_sc(.05,.95,sprintf('Baseline: %s: %s-%.0f',...
          run_legend_str{1},alg_baseline,N0_baseline),11);
  grid;
  print_fig(png);
end

if (1 && isinf(snr_db) && n_N0>1 && i_N0==n_N0)
  if isempty(mod_str)
    png.file_name = sprintf('05a-fastDD-relSNR-db-%s',config_str0_fn);
  else
    png.file_name = sprintf('05a-fastDD-relSNR-db-%s-%s',config_str0_fn,mod_str_fn);
  end
  png.file_count = plot_ID-1;
  figure(1); clf;
  plot(1,mean_rel_snr_baseline_db,'^','LineWidth',2); hold on;
  legend_alg_str{1} = 'Baseline';
  for i_N0_ = 1:n_N0
    N0_ = N0_list(i_N0_);
    plot([2:n_case]',mean_rel_snr_fastDD_db(i_N0_,2:n_case),'-*','LineWidth',2); hold on;
    %legend_alg_str{1+i_N0_} = sprintf('fastDD-%.0f',N0_);
    legend_alg_str{1+i_N0_} = sprintf('fastDD%.0f-%.0f',fastDD_alg_ID,N0_);
  end
  xlim([0 n_case+1]);
  set(gca,'XTick',[1:n_case]);
  set(gca,'XTickLabel',{'Baseline',run_legend_str{2:n_case}});
  ylim([-3 4]);
  xtickangle(45);
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  %xlabel('Normalized SNR dB')
  ylabel('Relative SNR dB')
  title(sprintf('fastDD Relative SNR vs PFB Configuration: %s, N0=%.0f-%.0f',config_str0,...
        min(N0_list),max(N0_list)));
  legend(legend_alg_str,'Location','NorthWest');
  text_sc(.30,.95,sprintf('Baseline: %s: %s-%.0f',...
          run_legend_str{1},alg_baseline,N0_baseline),11);
  text_sc(.30,.90,['N0=' sprintf('%.0f ',N0_list)],11);
  grid;
  print_fig(png);
end

if (1 && isinf(snr_db) && n_N0>1 && i_N0==n_N0)
  if isempty(mod_str)
    png.file_name = sprintf('05b-taylor-relSNR-db-%s',config_str0_fn);
  else
    png.file_name = sprintf('05b-taylor-relSNR-db-%s-%s',config_str0_fn,mod_str_fn);
  end
  png.file_count = plot_ID-1;
  figure(1); clf;
  plot(1,mean_rel_snr_baseline_db,'^','LineWidth',2); hold on;
  legend_alg_str{1} = 'Baseline';
  for i_N0_ = 1:n_N0
    N0_ = N0_list(i_N0_);
    plot([2:n_case]',mean_rel_snr_taylor_db(i_N0_,2:n_case),'-*','LineWidth',2); hold on;
    legend_alg_str{1+i_N0_} = sprintf('Taylor-%.0f',N0_);
  end
  xlim([0 n_case+1]);
  set(gca,'XTick',[1:n_case]);
  set(gca,'XTickLabel',{'Baseline',run_legend_str{2:n_case}});
  ylim([-3 4]);
  xtickangle(45);
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  %xlabel('Normalized SNR dB')
  ylabel('Relative SNR dB')
  title(sprintf('Taylor Relative SNR vs PFB Configuration: %s, N0=%.0f-%.0f',config_str0,...
        min(N0_list),max(N0_list)));
  legend(legend_alg_str,'Location','NorthWest');
  text_sc(.30,.95,sprintf('Baseline: %s: %s-%.0f',...
          run_legend_str{1},alg_baseline,N0_baseline),12);
  grid;
  print_fig(png);
end

if (1 && isinf(snr_db))
  if isempty(mod_str)
    png.file_name = sprintf('05c-taylor-fastDD-relSNR-db-%s-N0-%.0f',config_str0_fn,N0);
  else
    png.file_name = sprintf('05c-taylor-fastDD-relSNR-db-%s-N0-%.0f-%s',config_str0_fn,N0,mod_str_fn);
  end
  png.file_count = plot_ID-1;
  figure(1); clf;
  plot(1,mean_rel_snr_baseline_db,'^','LineWidth',2); hold on;
  plot([2:n_case]',[mean_rel_snr_taylor_db(i_N0,2:n_case) ; ...
       [mean_rel_snr_fastDD_db(i_N0,2:n_case)]]','-*','LineWidth',2);
  xlim([0 n_case+1]);
  set(gca,'XTick',[1:n_case]);
  set(gca,'XTickLabel',{'Baseline',run_legend_str{2:n_case}});
  ylim([-3 4]);
  xtickangle(45);
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  %xlabel('Normalized SNR dB')
  ylabel('Relative SNR dB')
  title(sprintf('Relative SNR vs PFB Configuration: %s, N0=%.0f',config_str0,N0));
  legend('Baseline',sprintf('Taylor-%.0f',N0_taylor),sprintf('fastDD%.0f-%.0f',fastDD_alg_ID,N0),...
         'Location','NorthWest');
  text_sc(.30,.95,sprintf('Baseline: %s: %s-%.0f',...
          run_legend_str{1},alg_baseline,N0_baseline),12);
  grid;
  print_fig(png);
end

if (1 & ~any(isinf(snr_db)) & (n_snr>1) & (n_case>2))
  if isempty(mod_str)
    png.file_name = sprintf('11-taylor-fastDD-dMDL-db-%s-N0-%.0f',config_str1_fn,N0);
  else
    png.file_name = sprintf('11-taylor-fastDD-dMDL-db-%s-N0-%.0f-%s',config_str1_fn,N0,mod_str_fn);
  end
  png.file_count = plot_ID-1;
  figure(1); clf;
  bar([1:n_case-1]',[delta_det_snr_taylor_db(i_N0,2:n_case) ; ...
      delta_det_snr_fastDD_db(i_N0,2:n_case)]');
  xlim([0 n_case]);
  set(gca,'XTick',[1:n_case-1]);
  set(gca,'XTickLabel',run_legend_str(2:n_case));
  ylim([0 5]);
  xtickangle(45);
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  %xlabel('Normalized SNR dB')
  ylabel('MDL Gain Over Baseline dB')
  title(sprintf('Detection Gain vs PFB Configuration: %s, N0=%.0f',config_str0,N0));
  legend(sprintf('Taylor-%.0f',N0_taylor),sprintf('fastDD%.0f-%.0f',fastDD_alg_ID,N0),...
         'Location','NorthEast');
  text_sc(.05,.95,sprintf('Baseline: %s: %s-%.0f',...
          run_legend_str{1},alg_baseline,N0_baseline),12);
  grid;
  print_fig(png);
end

if (1 & ~any(isinf(snr_db)) & (n_snr>1))
  if isempty(mod_str)
    png.file_name = sprintf('12-taylor-fastDD-MDL-db-%s-N0-%.0f',config_str1_fn,N0);
  else
    png.file_name = sprintf('12-taylor-fastDD-MDL-db-%s-N0-%.0f-%s',config_str1_fn,N0,mod_str_fn);
  end
  png.file_count = plot_ID-1;
  figure(1); clf;
  plot(1,det_snr_baseline_db,'^','LineWidth',2); hold on;
  plot([2:n_case]',[det_snr_taylor_db(i_N0,2:n_case) ; ...
       [det_snr_fastDD_db(i_N0,2:n_case)]]','-*','LineWidth',2);
  xlim([0 n_case+1]);
  set(gca,'XTick',[1:n_case]);
  set(gca,'XTickLabel',{'Baseline',run_legend_str{2:n_case}});
  ylim([-6 3]);
  xtickangle(45);
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  %xlabel('Normalized SNR dB')
  ylabel('Minimum Normalized SNR dB')
  title(sprintf('Minimum Detectable Level vs PFB Configuration: %s, N0=%.0f',config_str1,N0));
  legend('Baseline',sprintf('Taylor-%.0f',N0_taylor),sprintf('fastDD%.0f-%.0f',fastDD_alg_ID,N0),...
         'Location','NorthEast');
  text_sc(.05,.95,sprintf('Baseline: %s: %s-%.0f',...
          run_legend_str{1},alg_baseline,N0_baseline),12);
  grid;
  print_fig(png);
end

% ,'LineWidth',2
function ctl = set_symbol(i_case)

n_sym = 6;
i_ctl = mod(i_case-2,n_sym);

if (i_ctl==0)
  ctl = '-o';
elseif (i_ctl==1)
  ctl = '-^';
elseif (i_ctl==2)
  ctl = '-v';
elseif (i_ctl==3)
  ctl = '-d';
elseif (i_ctl==4)
  ctl = '-s';
else
  ctl = '->';
end

end