%
% script to check Taylor (TurboSETI) and fastDD recursion tree error characteristics
% Taylor delays are exact, high confidence
% fastDD delays are a model, low confidence (don't publish)
%

clear;

addpath('..\misc_fns');

output_directory=sprintf('./plots-taylor-fastDD-error/');
if ~exist(output_directory,'dir'), mkdir(output_directory); end

%default_figure_position = [775 43 758 591];	% set figure size for png import to pptx
default_figure_position = [1001 43 758 591];	% set figure size for png import to pptx

figure(1);
set(1,'Position',default_figure_position);
commandwindow;

PlotLineWidth = .6;
AxisFontSize = 12;
AxisFontWeight = 'bold';
EnableFinalPlots = 1;


png = print_fig;
png.output_directory = output_directory;

exp_min1 = 1;
exp_min2 = 5;
exp_max = 10;

n_error_all = NaN(exp_min2,exp_max,2,2);
pct_error_all = NaN(exp_min2,exp_max,2,2);

for Lf = [1 2]
%for Lf = [1]

  for alg_idx = [1,2]
    if alg_idx==1
      alg_str = 'Taylor';
    else
      alg_str = 'fastDD';
    end
    
    Lf_str = sprintf('%.0fx',Lf);
    fprintf('\n%s PFB-%s:\n',alg_str,Lf_str)

    for exp_min = exp_min1:exp_min2

      N0 = 2.^exp_min;
      legend_str{exp_min} = sprintf('N0=%.0f',N0);

      for exp_i = exp_min:exp_max

        N=2.^exp_i;
        
        if strcmp(alg_str,'Taylor')
          freq_bin_offset_N = freq_bin_offset(N);
          freq_bin_offset_N2 = freq_bin_offset(N/2);
          FBO_N = freq_bin_offset_N*(N-1);
        else
          freq_bin_offset_N = freq_bin_offset2(N);
          freq_bin_offset_N2 = freq_bin_offset2(N/2);
          FBO_N = freq_bin_offset_N*N;
        end
        int_bin_offset_N  = round(freq_bin_offset_N);
        int_bin_offset_N2 = round(freq_bin_offset_N2);

        if (N==N0)
          if (Lf==1)
            tree_bin_ofs_N = int_bin_offset_N;
          else
            tree_bin_ofs_N = round(Lf*freq_bin_offset_N)/Lf;
          end
        else
          if strcmp(alg_str,'Taylor')
            tree_bin_ofs_N2 = tree_bin_ofs_N;
            tree_bin_ofs_N = zeros(N,N);
            for i=1:N/2
              tree_bin_ofs_N(2*i-1,1:N/2) = tree_bin_ofs_N2(i,:);
              tree_bin_ofs_N(2*i,1:N/2)   = tree_bin_ofs_N2(i,:);
              ii = N/2+1:N;
              tree_bin_ofs_N(2*i-1,ii) = tree_bin_ofs_N2(i,:)+i-1;
              tree_bin_ofs_N(2*i,ii)   = tree_bin_ofs_N2(i,:)+i;
            end
          else
            tree_bin_ofs_N2 = tree_bin_ofs_N;
            tree_bin_ofs_N = zeros(N+1,N);
            tree_bin_ofs_N(1:2:N+1,:) = [tree_bin_ofs_N2 [0:N/2]' + tree_bin_ofs_N2];
            tree_bin_ofs_lower = [tree_bin_ofs_N2(1:N/2,:) [0:N/2-1]' + tree_bin_ofs_N2(2:N/2+1,:)];
            tree_bin_ofs_upper = [tree_bin_ofs_N2(2:N/2+1,:) [1:N/2]' + tree_bin_ofs_N2(1:N/2,:)];
            tree_bin_ofs_mid = (tree_bin_ofs_lower + tree_bin_ofs_upper)/2;
            tree_bin_ofs_delta = (tree_bin_ofs_upper-tree_bin_ofs_lower)/2;
            tree_bin_ofs_N(2:2:N,:) = tree_bin_ofs_lower;
            %tree_bin_ofs_N(2:2:N,:) = tree_bin_ofs_upper;
            %rand_draw = .5*ones(N/2,N);
            %rand_draw = randi([0 1],N/2,1)*ones(1,N);
            %tree_bin_ofs_N(2:2:N,:) = rand_draw.*tree_bin_ofs_lower + (1-rand_draw).*tree_bin_ofs_upper;
            %sigma_offset = .50;
            %rand_offset = sigma_offset*randn(N/2,1)*ones(1,N);
            %rand_offset = (2*rand(N/2,1)-1)*ones(1,N);      % +/- 100% of upper-lower span
            %rand_offset = .5*(2*rand(N/2,1)-1)*ones(1,N);  % +/- 50% of upper-lower span
            %tree_bin_ofs_N(2:2:N,:) = tree_bin_ofs_mid + rand_offset.*tree_bin_ofs_delta;
          end
        end

        bin_error_N = tree_bin_ofs_N - freq_bin_offset_N;
        error_N = round(abs(bin_error_N)-.001);
        n_error_tree = length(find(abs(bin_error_N)>0.5));
        pct_error_tree = n_error_tree/length(error_N(:)) * 100;

        n_error_all(exp_min,exp_i,Lf,alg_idx) = n_error_tree;
        pct_error_all(exp_min,exp_i,Lf,alg_idx) = pct_error_tree;
        
        if (0)
          if (N<=32)
            print_compact_matrix('tree_bin_ofs_N*N',tree_bin_ofs_N*N,4,0);
            print_compact_matrix('bin_error_N',bin_error_N,4,2);
            print_compact_matrix('error_N',error_N,2,0);
          end
        end

        results_string = ...
          sprintf('N=%5.0f, RMS error=%.2f, min, max, mean, std error=%6.2f %6.2f %6.2f %6.2f, pct Q error=%.0f',...
            N,rms(bin_error_N(:)),min(bin_error_N(:)),max(bin_error_N(:)),...
            mean(bin_error_N(:)),std(bin_error_N(:)),pct_error_tree);
        fprintf('%s-PFB-%s %s\n',alg_str,Lf_str,results_string);

        if (1)
          png.file_name = sprintf('01-%s-error-hist-%s-N-%.0f-N0-%.0f',alg_str,Lf_str,N,N0);
          png.file_count = 0;
          figure(1); clf;
          hist(bin_error_N(:),101);
          v=axis; axis([-1 1 v(3:4)]);
          set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
          xlabel('Fractional Error (Frequency Bins)');
          ylabel('Relative Frequency');
          title(sprintf('%s DeDoppler %s Fractional Error, N=%.0f, N0=%.0f',alg_str,Lf_str,N,N0));
          text_sc(.05,.95,results_string);
          grid;
          png = print_fig(png);
          %pause(2);
        end

        if (1)
          png.file_name = sprintf('02-%s-error-spy-%s-N-%.0f-N0-%.0f',alg_str,Lf_str,N,N0);
          png.file_count = 0;
          figure(1); clf;
          spy(error_N);
          set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
          xlabel('Drift Rate Bin');
          ylabel('Frequency Bin');
          title(sprintf('%s DeDoppler %s Quantization Error, N=%.0f, N0=%.0f\n%s',alg_str,Lf_str,N,N0,results_string));
          grid;
          png = print_fig(png);
          %pause(2);
        end

        if (0)
          png.file_name = sprintf('03-%s-error-spy-%s-N-%.0f-N0-%.0f',alg_str,Lf_str,N,N0);
          png.file_count = 0;
          figure(1); clf;
          %imagesc([0:N-1],[0:N-1],tree_bin_ofs_N,[0 N]); % no "axis xy"
          imagesc([0:N-1],[0:N-1],abs(bin_error_N),[0 1.5]); % no "axis xy"
          colorbar; colormap(jet);
          colorbar_label([],'Error Freq Bins');
          xlabel('Time (Samples)')
          ylabel('Drift Rate Index')
          title(sprintf('%s DeDoppler %s Quantization Error, N=%.0f, N0=%.0f\n%s',alg_str,Lf_str,N,N0,results_string));
          png = print_fig(png);
          pause(2);
        end


      end  % N loop
      fprintf('\n');
    end  % exp_min loop
    
    if (1)
      png.file_name = sprintf('00-%s-pct-error-Lf-%.0f',alg_str,Lf);
      png.file_count = 0;
      figure(1); clf;
      plot(exp_min1:exp_max,pct_error_all(:,:,Lf,alg_idx),'-*','LineWidth',2);
      v=axis; axis([v(1:2) 0 30]);
      set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
      if (0)
        xlabel('log2(N_T)');
      else
        xlabel('Number of Integrations N_T');
        for i=exp_min1:exp_max, N_T_label{i}=sprintf('%.0f',2.^i); end
        set(gca,'XTick',exp_min1:exp_max);
        set(gca,'XTickLabel',N_T_label);
      end
      ylabel('Percentage Indexing Error');
      title(sprintf('%s DeDoppler Percent Indexing Error vs. N_T, %s',alg_str,Lf_str));
      grid;
      legend(legend_str,'Location','NorthWest');
      png = print_fig(png);
    end

  end % alg_str
end  % Lf

function FBO = freq_bin_offset(N)

FBO = [0:N-1]'*[0:N-1]/(N-1);  % Nominal Taylor TurboSETI alg

end

function FBO2 = freq_bin_offset2(N)

FBO2 = [0:N]'*[0:N-1]/N;    % fastDD alg

end


function print_compact_matrix(name_string,A,n_left,n_right);

  fprintf('%s\n',name_string);
  [m,n]=size(A);
  if (n_right>0)
    p_spec = sprintf('%%%.0f.%.0ff ',n_left+n_right+1,n_right);
  else
    p_spec = sprintf('%%%.0f.0f ',n_left);
  end
  for i=1:m
    fprintf(p_spec,A(i,:));
    fprintf('\n');
  end
end