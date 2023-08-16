%
% script to check Taylor (TurboSETI) and fastDD recursion tree error characteristics
% Taylor delays are exact, high confidence
% fastDD delays are a model, low confidence (don't publish)
%

clear;

addpath('..\misc_fns');

output_directory=sprintf('./plots-taylor-fastDD-diagram/');
if ~exist(output_directory,'dir'), mkdir(output_directory); end

%default_figure_position = [775 43 758 591];	% set figure size for png import to pptx
default_figure_position = [1001 43 758 591];	% set figure size for png import to pptx

figure(1);
set(1,'Position',default_figure_position);
commandwindow;

png = print_fig;
png.output_directory = output_directory;

PlotLineWidth = .6;
AxisFontSize = 12;
AxisFontWeight = 'bold';
EnableFinalPlots = 1;

for Lf = [1 2]
%for Lf = [1]

%  for N1 = [4:16]
  for N1 = [4:8]

    %for alg_idx = [1,2,3,4]
    for alg_idx = [1,2]
      if alg_idx==1
        alg_str = 'Taylor';
      elseif alg_idx==2
        alg_str = 'fastDD';
      elseif alg_idx==3
        alg_str = 'Alt3';
      elseif alg_idx==4
        alg_str = 'Alt4';
      end

      Lf_str = sprintf('%.0fx',Lf);
      fprintf('\n%s PFB-%s:\n',alg_str,Lf_str)

      N2 = 2*N1;
      Nf = N2;

      if strcmp(alg_str,'Taylor')
        drift_grid1 = [0:N1-1]/(N1-1);
        drift_grid2 = [0:N2-1]/(N2-1);
      elseif strcmp(alg_str,'fastDD')
        drift_grid1 = [0:N1]/(N1);
        drift_grid2 = [0:N2]/(N2);
      elseif strcmp(alg_str,'Alt3')
        drift_grid1 = [0:N1+1]/(N1+1);
        drift_grid2 = [0:N2+1]/(N2+1);
      elseif strcmp(alg_str,'Alt4')
        drift_grid1 = [0:N1-2]/(N1-2);
        drift_grid2 = [0:N2-2]/(N2-2);
      end
      f_grid_ctr = [0:1/Lf:Nf-1]+.5/Lf;
      t_grid_ctr = [0:N2-1]+.5;
      f_grid_edge = [0:1/Lf:Nf];
      t_grid_edge = [0:N2];
      Nr1 = length(drift_grid1);
      Nr2 = length(drift_grid2);
      DG1 = [zeros(1,Nr1) ; drift_grid1*N1];
      DG2 = [zeros(1,Nr2) ; drift_grid2*N2];

      if (1)
        png.file_name = sprintf('00-%s-grid_diagram-Lf-%.0f-N1-%.0f',alg_str,Lf,N1);
        png.file_count = 0;
        figure(1); clf;
        subplot(1,2,1)   % Even Drift Rates
        plot_grid;
        % N1 stage drift lines
        if strcmp(alg_str,'Taylor')
          plot(DG1,[0 N1]','-r','LineWidth',1); hold on;
          plot(DG1+[0:Nr1-1],[N1 N2]','-r','LineWidth',1); hold on;
        elseif strcmp(alg_str,'fastDD')
          plot(DG1,[0 N1]','-r','LineWidth',1); hold on;
          plot(DG1+[0:Nr1-1],[N1 N2]','-r','LineWidth',1); hold on;
        elseif strcmp(alg_str,'Alt3')
          plot(DG1,[0 N1]','-r','LineWidth',1); hold on;
          plot(DG1(:,1:end-1)+[0:Nr1-2],[N1 N2]','-r','LineWidth',1); hold on;
        elseif strcmp(alg_str,'Alt4')
          plot(DG1,[0 N1]','-r','LineWidth',1); hold on;
          plot([DG1 DG1(:,Nr1) DG1(:,Nr1)]+[0:Nr1+1],[N1 N2]','-r','LineWidth',1); hold on;
        end
        % N2 stage drift lines
        plot(DG2(:,1:2:end),[0 N2]','--b','LineWidth',1); hold on;
        v=axis; axis([-1 Nf+1 -1 N2+1]); axis ij; %axis equal;
        set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
        xlabel('Frequency/fs');
        ylabel('Time/Ts');
        title(sprintf('%s DeDoppler, Lf=%.0f, N1=%.0f\n"Even" Drift Rates',alg_str,Lf,N1));
        subplot(1,2,2)   % Odd Drift Rates
        plot_grid;
        % N1 stage drift lines
        if strcmp(alg_str,'Taylor')
          plot(DG1,[0 N1]','-r','LineWidth',1); hold on;
          plot(DG1+[1:Nr1],[N1 N2]','-r','LineWidth',1); hold on;
        elseif strcmp(alg_str,'fastDD')
          plot(DG1,[0 N1]','-r','LineWidth',1); hold on;
          plot(DG1(:,1:end-1)+[1:Nr1-1],[N1 N2]','-r','LineWidth',1); hold on;
          plot(DG1(:,2:end)+[0:Nr1-2],[N1 N2]','-r','LineWidth',1); hold on;
        elseif strcmp(alg_str,'Alt3')
          plot(DG1,[0 N1]','-r','LineWidth',1); hold on;
          plot(DG1(:,2:end)+[0:Nr1-2],[N1 N2]','-r','LineWidth',1); hold on;
        elseif strcmp(alg_str,'Alt4')
          plot(DG1,[0 N1]','-r','LineWidth',1); hold on;
          plot(DG1+[1:Nr1],[N1 N2]','-r','LineWidth',1); hold on;
        end
        % N2 stage drift lines
        plot(DG2(:,2:2:end),[0 N2]','--b','LineWidth',1); hold on;
        v=axis; axis([-1 Nf+1 -1 N2+1]); axis ij; %axis equal;
        set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
        xlabel('Frequency/fs');
        ylabel('Time/Ts');
        title(sprintf('"Odd" Drift Rates'));
        %grid;
        png = print_fig(png);
      end

    end % alg_str

  end % N1
end  % Lf

