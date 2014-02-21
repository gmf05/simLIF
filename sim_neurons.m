function [V,spikes,time] = sim_neurons()

% get parameters--------------------------------------------------
global dt
global NT
global E_cell_dim
global I_cell_dim
global network_topology
global C_e
global C_i
global gL_e
global gL_i
global EL_e
global EL_i
global I0
global Vmax
global Vth
global V0
global T
global Trefract_e
global Trefract_i
global g_e_max
global g_i_max
global E_e
global E_i
global beta_e
global beta_i
global alpha_max_e
global alpha_max_i
global k_e
global k_i
global sigma_e
global sigma_i
global W_ee
global W_ei
global W_ie
global W_ii
global g_tr_high
global tau_tr_e
global tau_tr_i
global Etr
global P0
global tau_rel_e
global tau_rel_i
global fd_e
global fd_i
get_sim_params();

N_E_cells = prod(E_cell_dim);
N_I_cells = prod(I_cell_dim);
N_cells = N_E_cells + N_I_cells;

% concatenate certain E/I parameters into lists for ease later on
concat_EI = @(x,y)[x*ones(N_E_cells,1); y*ones(N_I_cells,1)];
C = concat_EI(C_e, C_i);
gL = concat_EI(gL_e, gL_i);
EL = concat_EI(EL_e, EL_i);
tau_tr = concat_EI(tau_tr_e, tau_tr_i);
refs = concat_EI(Trefract_e, Trefract_i);
sigma = [sigma_e sigma_i];
W = [W_ee W_ie; W_ei W_ii];


% setup coordinates
coord = zeros(N_cells, 2);

count = 0;
for i = 1:E_cell_dim(1)
  coord(count + (1:E_cell_dim(2)), 1) = i;
  coord(count + (1:E_cell_dim(2)), 2) = 1:E_cell_dim(2);
  count = count + E_cell_dim(2);
end

for i = 1:I_cell_dim(1)
  coord(count + (1:I_cell_dim(2)), 1) = i;
  coord(count + (1:I_cell_dim(2)), 2) = 1:I_cell_dim(2);
  count = count + I_cell_dim(2);
end

% calculate distances
rho = unifrnd(0.5,1.5,[1,N_cells]); % randomness of connection strength
D = zeros(N_cells);

% distance function depends on network topology
switch network_topology
  case 'sheet'
    network_dist = @(X1,X2) norm(X1-X2,2);
    
  case' tube'
    []; % NEED TO FILL THIS IN
    network_dist = @(X1,X2) norm(X1-X2,2);
end


% loop over cell pairs
isIcell1 = 0;
for c1 = 1:N_cells
  if c1==N_E_cells+1, isIcell1 = 1; end
  isIcell2 = 0;
  i = coord(c1,1);  j = coord(c1,2);
  
  for c2 = 1:N_cells
    if c2==N_E_cells+1, isIcell2 = 1; end
    l = coord(c2,1);  m = coord(c2,2);
    Wc = W(1+isIcell1, 1+isIcell2);
    d = network_dist([i,j],[l,m]);
    D(c1,c2) = Wc*rho(c1)*exp(-d^2/sigma(1+isIcell2)^2);
  end
  
end

% initialize arrays-------------------------------------------
V = zeros(N_cells, NT); % all voltages
spikes = zeros(N_cells, NT); % raster plot
Vt = concat_EI(-73.6, -81.6); % current voltage initialized to val1 for E, val2 for I
last_spike = zeros(N_cells, 1); % when was each cell's last spike?
g_tr = zeros(N_cells, 1); % after-hyperpol. conductance
g_e = zeros(N_cells, 1); % synaptic conductance from E cells
g_i = zeros(N_cells, 1); % synaptic conductance from I cells
y = zeros(N_cells, 1); % spiking at each time step
% arrays used to compute gE, gI:
P_post_e = zeros(N_cells, 1); % post-syn. channel open prob.
P_post_i = zeros(N_cells, 1); % post-syn. channel open prob.
P_rel_e = P0 * ones(N_cells, 1); % synaptic release prob.
P_rel_i = P0 * ones(N_cells, 1); % synaptic release prob.
I = zeros(N_cells, 1);

% target input current to middle of grid
arrayMap = reshape(1:N_E_cells, E_cell_dim);
mid1 = 0.5*E_cell_dim(1);
mid2 = 0.5*E_cell_dim(2);
I(reshape(arrayMap(mid1-1:mid1+1,mid2-1:mid2+1),1,[])) = I0;

% main time loop----------------------------------------------
for t = 1:NT
  if mod(t,100)==0, fprintf(['t = ' num2str(t) '\n']), end
  % compute voltages (Eqn 1)
  
% % % %     alph=50*ones(N_cells,1);        
% % % %     alpha1=(alph+g_e+g_i);
% % % %     beta1=(14/3)*g_e-(2/3)*g_i;
% % % % 
% % % %     local_g2=exp(-dt/2/tau_ex)*g+(dt/2/tau_ex)*exp(-dt/2/tau_ex)*w;
% % % %     local_g_e2=g_ext+Ex*local_g2;
% % % %     local_g_inhibt2=exp(-dt/2/tau_in)*g_inhibt+(dt/2/tau_in)*exp(-dt/2/tau_in)*w_inhibt;
% % % %     local_g_i2=In*local_g_inhibt2+noise;
% % % %     alpha2=(alph+local_g_e2+local_g_i2);
% % % %     beta2=(14/3)*local_g_e2-(2/3)*local_g_i2;
% % % %     
% % % %     local_g3=exp(-dt/tau_ex)*g+(dt/tau_ex)*exp(-dt/tau_ex)*w;
% % % %     local_g_e3=g_ext+Ex*local_g3;
% % % %     local_g_inhibt3=exp(-dt/tau_in)*g_inhibt+(dt/tau_in)*exp(-dt/tau_in)*w_inhibt;
% % % %     local_g_i3=In*local_g_inhibt3+noise;
% % % %     alpha3=(alph+local_g_e3+local_g_i3);
% % % %     beta3=(14/3)*local_g_e3-(2/3)*local_g_i3;
  
  % Fourth order Runge-Kutta method to integrate V
  % A's and B's come from rewriting dV/dt
  % and evaluating at t, t+dt/2, and t+dt, resp.
  A1 = 1./C.*(gL + g_e + g_i + g_tr);
  B1 = 1./C.*(I + EL.*gL + Etr*g_tr + E_e*g_e + E_i*g_i);
  rk1 = -A1.*Vt+B1;
  
  g_e2 = g_e;
  g_i2 = g_i;
  g_tr2 = g_tr - dt/2*g_tr./tau_tr;
  A2 = 1./C.*(gL + g_e2 + g_i2 + g_tr2);
  B2 = 1./C.*(I + EL.*gL + Etr*g_tr2 + E_e*g_e2 + E_i*g_i2);
  rk2 = -A2.*(Vt + rk1*dt/2)+B2;
  rk3 = -A2.*(Vt + rk2*dt/2)+B2;
  
  g_e3 = g_e;
  g_i3 = g_i;
  g_tr3 = g_tr - dt*g_tr./tau_tr;
  A3 = 1./C.*(gL + g_e3 + g_i3 + g_tr3);
  B3 = 1./C.*(I + EL.*gL + Etr*g_tr3 + E_e*g_e3 + E_i*g_i3);
  rk4 = -A3.*(Vt + rk3*dt)+B3;
  
  Vt = Vt + dt/6*(rk1+2*rk2+2*rk3+rk4);
  
  % correct membrane potential where Vt > Vth
  Vth_ind = find(Vt > Vth);
  for i = Vth_ind'
    isIcell = (i > N_E_cells);
    if time(t) < last_spike(i) + T
      Vt(i) = Vmax;
    elseif time(t) < last_spike(i) + T + refract(1+isIcell)
      Vt(i) = V0;
    else
      last_spike(i) = time(t);
      Vt(i) = Vmax;
    end
  end
  
  % save voltage, spikes
  V(:,t) = Vt;
  y = y*0;
  y(Vt == Vmax) = 1;
  spikes(last_spike==time(t), t) = 1;
  
  % compute after-hyperpolarization conductances (Eqn 2)
  % first find which cells are refractory
  if time(t) > Trefract_e + T
    ghigh_ind = find(last_spike + T + refs > time(t));
  else
    ghigh_ind = [];
  end
  g_ind = setdiff(1:N_cells, ghigh_ind);
  g_tr(ghigh_ind) = g_tr_high;
  % for non-refractory cells, integrate g_tr
  g_tr(g_ind) = g_tr(g_ind) - dt*g_tr(g_ind)./tau_tr(g_ind); % SWITCH TO RUNGE-KUTTA??
  
  % compute new synaptic conductances (Eqns 3-7)
  % compute S_s (Eqn 7)
  S_e = D(:,1:N_E_cells)*y(1:N_E_cells);
  S_i = D(:,N_E_cells+1:end)*y(N_E_cells+1:end);
  
  % compute alpha_s (Eqn 6)
  alpha_e = alpha_max_e * (1 - exp(-S_e ./ k_e));
  alpha_i = alpha_max_i * (1 - exp(-S_i ./ k_i));
  
  % compute P_post_s (Eqn 5)
  P_post_e = P_post_e + dt*(alpha_e .* (1-P_post_e) - beta_e * P_post_e);
  P_post_i = P_post_i + dt*(alpha_i .* (1-P_post_i) - beta_i * P_post_i);
  
  % compute P_rel_s (Eqn 4)
  % NOTE: SEE BEFORE EQN 4 - do we need to multiply for fd_s at spiking cells?
  P_rel_e = P_rel_e + dt/tau_rel_e*(P0 - P_rel_e);
  P_rel_i = P_rel_i + dt/tau_rel_i*(P0 - P_rel_i);
  
  % compute g_s (Eqn 3)
  g_e = g_e_max * P_post_e .* P_rel_e;
  g_i = g_i_max * P_post_i .* P_rel_i;
  
end

end