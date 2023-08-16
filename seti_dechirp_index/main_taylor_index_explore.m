%
% script to check Taylor (TurboSETI) and fastDD recursion tree error characteristics
% Taylor delays are exact, high confidence
%

clear;

addpath('..\misc_fns');

output_directory=sprintf('./plots-taylor-fastDD-error/');
if ~exist(output_directory,'dir'), mkdir(output_directory); end

% %default_figure_position = [775 43 758 591];	% set figure size for png import to pptx
% default_figure_position = [1001 43 758 591];	% set figure size for png import to pptx
% 
% figure(1);
% set(1,'Position',default_figure_position);
% commandwindow;
% 
% png = print_fig;
% png.output_directory = output_directory;

if (1)
  Lf_eq_Lr = 1;
  Lf_string = 'Lf=Lr'; 
elseif (0)
  Lf = 2;
  Lf_eq_Lr = 0;
  Lf_string = 'Lf=2'; 
elseif (1)
  Lf = 1;
  Lf_eq_Lr = 0;
  Lf_string = 'Lf=1'; 
end
fprintf('\nExtended Taylor DD cases for %s\n',Lf_string);

for Nt0 = [4 5 6 8]

%for Nt_final=[16 32 64 128 256 512 1024]
%for Nt_final=[16 32 64 128 256]
%for Nt_final=[32]
%for Nt_final=[1024]
for Nt_final=[640]
  n_stage = round(log2(Nt_final/Nt0));
  Nt = Nt0*2^n_stage;
  for Nr0 = [Nt0 floor(1.4*Nt0):ceil(1.6*Nt0) 2*Nt0]
  %for Nr0 = Nt0 + [0:1];
    Lr = Nr0/Nt0;
    if Lf_eq_Lr
      Lf = Lr;
    end
    n_stage = round(log2(Nt/Nt0));
    Nr = (2.^n_stage)*Nr0;
    Nt2 = Nt/2;
    %fprintf('\nNt = %.0f, Nt0 = %.0f, Nr0 = %.0f, Lr = %.2f\n',Nt,Nt0,Nr0,Lr);

    freq_bin_ofs_Nt = freq_bin_ofs(Nt,Nr);
    freq_bin_ofs_Nt2 = freq_bin_ofs(Nt/2,Nr/2);
    freq_bin_ofs_Nt0 = freq_bin_ofs(Nt0,Nr0);
    int_bin_ofs_Nt  = round(freq_bin_ofs_Nt);
    int_bin_ofs_Nt2 = round(freq_bin_ofs_Nt2);
    int_bin_ofs_Nt0 = round(freq_bin_ofs_Nt0);

    % quantized error without fast algorithm

    rmse_direct = rms(int_bin_ofs_Nt(:)-freq_bin_ofs_Nt(:));
    epct_direct = find_epct(int_bin_ofs_Nt(:)-freq_bin_ofs_Nt(:));

    % unquantized Taylor error with fast algorithm

    [tree_ofs_taylorR,rmse_taylorR,error_taylorR] = ...
                    taylor_freq_ofs_forward(Nt,Nt0,Nr0,Lf,freq_bin_ofs_Nt0);
    epct_taylorR = find_epct(error_taylorR);
    mean_errR = mean(error_taylorR(:));

     % nominal Taylor Nt from int_bin_ofs_Nt0

    [tree_ofs_taylorQ,rmse_taylorQ,error_taylorQ] = ...
                    taylor_freq_ofs_forward(Nt,Nt0,Nr0,Lf,int_bin_ofs_Nt0);
    epct_taylorQ = find_epct(error_taylorQ);
    mean_errQ = mean(error_taylorQ(:));

    i12 = [1:2:Nr 2:2:Nr];
    [temp,i21] = sort(i12);
    P12 = eye(Nr); 
    P12 = P12(i12,:); % normal to [even rows ; odd rows]N
    P21 = eye(Nr); 
    P21 = P21(i21,:); % [even rows ; odd rows] to normal
% 
%     error_Nt0_Q = int_bin_ofs_Nt0 - freq_bin_ofs_Nt0;
    error_taylorR_even_odd = P12*error_taylorR;
    %mean(error_taylorR)

    fprintf('Nt=%4.0f Nr=%4.0f Nt0=%2.0f Nr0=%2.0f Lf=%.2f Lr=%.2f %.1f rmse R %.2f Q %.2f mean R %5.2f Q %5.2f Pct R %.0f Q %.0f\n',...
               Nt,Nr,Nt0,Nr0,Lf,Lr,Lf*Lr,rmse_taylorR,rmse_taylorQ,mean_errR,mean_errQ,epct_taylorR,epct_taylorQ);
    if (0)
      if (mod(log2(Nr0),1)==0)
        N_fft = 16*1024;
        M = N_fft/Lf;
        N_zp = (Lf-1)*M;
      else
        M = 16*1024;
        N_fft = M*Lf;
        N_zp = (Lf-1)*M;
      end
      fprintf('  Lf=%.2f N_fft==%.0f M=%.2f N_zp=%.2f\n',Lf,N_fft,M,N_zp);
    end
  end
  fprintf('\n');
end
end