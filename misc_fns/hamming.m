function w = hamming(N);



w = .54 - .46*cos(2*pi/(N-1)*[0:N-1]');