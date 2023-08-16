function epct = find_epct(error_N)

error_N = error_N(:);
ii = find(abs(error_N)>.5);

epct = length(ii)/length(error_N)*100;

end
