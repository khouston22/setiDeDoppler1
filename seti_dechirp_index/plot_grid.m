%
% for diagram, plot grid lines and points
%

% horiz grid line
for i_t=1:length(t_grid_edge)
  plot(f_grid_edge([1 end]),t_grid_edge(i_t)*[1 1],'-c','LineWidth',.1); hold on;
end
% vert grid line
for i_f=1:length(f_grid_edge)
  plot((f_grid_edge([i_f]))*[1 1],t_grid_edge([1 end]),'-c','LineWidth',.1); hold on;
end
% grid center points
for i_t=1:length(t_grid_ctr)
  plot(f_grid_ctr,t_grid_ctr(i_t),'.k','LineWidth',1); hold on;
end
