function window_name = choose_lpf_window(Lf,N_tap);
%
% function to choose lpf window based on Lf and subfilter length
% attempts to keep crossover at ~-4 dB or higher
%

if (Lf==1)
  if (N_tap<=4)
    window_name = 'lpf114';
  elseif (N_tap<=6)
    window_name = 'lpf109';
  elseif (N_tap<=8)
    window_name = 'lpf107';
  end
elseif(Lf<1.3)
  % Lf ~ 1.25
  if (N_tap<=4)
    window_name = 'lpf094';
  elseif (N_tap<=6)
    window_name = 'lpf090';
  elseif (N_tap<=8)
    window_name = 'lpf088';
  end
elseif(Lf<1.5)
  % Lf ~ 1.33
  if (N_tap<=4)
    window_name = 'lpf092';
  elseif (N_tap<=6)
    window_name = 'lpf088';
  elseif (N_tap<=8)
    window_name = 'lpf086';
  end
elseif(Lf<1.7)
  % Lf ~ 1.5
  if (N_tap<=4)
    window_name = 'lpf080';
  elseif (N_tap<=6)
    window_name = 'lpf076';
  elseif (N_tap<=8)
    window_name = 'lpf074';
  end
elseif(Lf<2)
  % Lf ~ 1.75
  if (N_tap<=4)
    window_name = 'lpf070';
  elseif (N_tap<=6)
    window_name = 'lpf067';
  elseif (N_tap<=8)
    window_name = 'lpf064';
  end
elseif(Lf<2.5)
  % Lf ~ 2.0
  if (N_tap<=4)
    window_name = 'lpf057';
  elseif (N_tap<=6)
    window_name = 'lpf060';
  elseif (N_tap<=8)
    window_name = 'lpf057';
  end
else
  window_name = 'lpf050';
end
