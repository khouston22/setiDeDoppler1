function w = hanning(N);


w = .5*(1 - cos(2*pi/(N+1)*[1:N]'));

