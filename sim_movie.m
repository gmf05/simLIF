function sim_movie(V,Cax)

if nargin<2
  Cax = [min(min(V)), max(max(V))];
%   Cax = [-75 -55]; % good for voltage
end

global dt
global NT
global E_cell_dim
global I_cell_dim

N_E_cells = prod(E_cell_dim);
time = (1:NT)*dt;
t_delay = .001; % pause used in plotting [seconds]
% define handles sp1, sp2 which will be subplots?
% sp1 = [2 1 1];
% sp2 = [2 1 2];

for t = 1:NT
  subplot(2,1,1)
  imagesc(reshape(V(1:N_E_cells,t),E_cell_dim));
  caxis(Cax);
  title(num2str(time(t)));

  subplot(2,1,2)
  imagesc(reshape(V(N_E_cells+1:end,t),I_cell_dim));
  caxis(Cax);
  
  pause(t_delay);
end