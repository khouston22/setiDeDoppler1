function [FBO,QFBO,errFBO,rms_err,mean_err] = freq_bin_ofs(Nt,Nr)

FBO = [0:Nr-1]'*[0:Nt-1]/(Nr-1);  % Nominal Taylor TurboSETI alg
QFBO = round(FBO);
errFBO = FBO-QFBO;
rms_err = rms(errFBO(:));
mean_err = mean(errFBO(:));
end

