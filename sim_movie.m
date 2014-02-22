function sim_movie(V)

global dt
global NT
global E_cell_dim
global I_cell_dim

N_E_cells = prod(E_cell_dim);
N_I_cells = prod(I_cell_dim);
time = (1:NT)*dt;
% define handles sp1, sp2 which will be subplots?
% sp1 = [2 1 1];
% sp2 = [2 1 2];
t_delay = .001;

% minV = -75; maxV = -55;
minV = min(min(V)); maxV = max(max(V));

for t = 1:NT
  subplot(2,1,1)
  imagesc(reshape(V(1:N_E_cells,t),E_cell_dim));
  caxis([minV,maxV]);
  title(num2str(time(t)));

  subplot(2,1,2)
  imagesc(reshape(V(N_E_cells+1:end,t),I_cell_dim));
  caxis([minV,maxV]);

  pause(t_delay);
end