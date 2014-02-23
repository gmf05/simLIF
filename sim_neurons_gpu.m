function [V,spikes,gs,time] = sim_neurons()

fprintf('Simulating network\n');
CLOCK = tic();

% get parameters--------------------------------------------------
fprintf('Getting parameters\n');
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
global ImVs
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
global E_tr
global P0
global tau_rel_e
global tau_rel_i
global fd_e
global fd_i
get_sim_params();
time = (1:NT)*dt;
N_E_cells = prod(E_cell_dim);
N_I_cells = prod(I_cell_dim);
N_cells = N_E_cells + N_I_cells;

toc(CLOCK);

% concatenate certain E/I parameters into lists for ease later on
concat_EI = @(x,y)[x*ones(N_E_cells,1); y*ones(N_I_cells,1)];
C = concat_EI(C_e, C_i);
gL = concat_EI(gL_e, gL_i);
EL = concat_EI(EL_e, EL_i);
tau_tr = concat_EI(tau_tr_e, tau_tr_i);
refract = [Trefract_e Trefract_i];
refvec = concat_EI(Trefract_e, Trefract_i);
% sigma = [sigma_e sigma_i];
% W = [W_ee W_ie; W_ei W_ii];

% setup coordinates
fprintf('Calculating cell geometry\n');
coord = zeros(N_cells, 2);

% NOTE: Assumes for now that there are more E cells than I cells
% and remaps coord of I cells to approximate pos. relative to E cells
% if there are more I cells, this mapping should be reversed
count = 0;
for i = 1:E_cell_dim(1)
  coord(count + (1:E_cell_dim(2)), 1) = i;
  coord(count + (1:E_cell_dim(2)), 2) = 1:E_cell_dim(2);
  count = count + E_cell_dim(2);
end

EIratio1 = E_cell_dim(1)/I_cell_dim(1);
EIratio2 = E_cell_dim(2)/I_cell_dim(2);
for i = 1:I_cell_dim(1)
  % no coordinate remapping (NOT advised)
%   coord(count + (1:I_cell_dim(2)), 1) = i;
%   coord(count + (1:I_cell_dim(2)), 2) = 1:I_cell_dim(2);
  
  % coordinate remapping
  coord(count + (1:I_cell_dim(2)), 1) = round(i*EIratio1);
  coord(count + (1:I_cell_dim(2)), 2) = round((1:I_cell_dim(2))*EIratio2);
  
  count = count + I_cell_dim(2);
end

% calculate distances
% NOTE: breaking distance & spike arrays
% into I & E subgroups speeds up computation by ~5x

% distance function depends on network topology
% switch network_topology
%   case 'sheet'
%     network_dist = @(X1,X2) norm(X1-X2,2);
%     
%   case' tube'
%     []; % NEED TO FILL THIS IN
%     network_dist = @(X1,X2) norm(X1-X2,2);
% end

% loop over cell pairs
% network_dist = @(i,j,l,m) sqrt((i-l)^2+(j-m)^2);
% load NN_dist_matrix;

% rho = unifrnd(0.5,1.5,[1,N_cells]);
load('rho_instance','rho');
tic
for c1 = 1:N_cells
  i = coord(c1,1);  j = coord(c1,2);  r = rho(c1);
  
  if c1 < N_E_cells+1
    for c2 = 1:N_E_cells
      l = coord(c2,1);  m = coord(c2,2);
      d = sqrt((i-l)^2+(j-m)^2); %FASTEST
      D_e(c1,c2) = W_ee*r*exp(-d^2/sigma_e^2);
    end
    for c2 = 1:N_I_cells
      l = coord(c2+N_E_cells,1);  m = coord(c2+N_E_cells,2);
      d = sqrt((i-l)^2+(j-m)^2); %FASTEST
      D_i(c1,c2) = W_ie*r*exp(-d^2/sigma_i^2);
    end
  else
    for c2 = 1:N_E_cells
      l = coord(c2,1);  m = coord(c2,2);
      d = sqrt((i-l)^2+(j-m)^2); %FASTEST
      D_e(c1,c2) = W_ei*r*exp(-d^2/sigma_e^2);
    end
    for c2 = 1:N_I_cells
      l = coord(c2+N_E_cells,1);  m = coord(c2+N_E_cells,2);
      d = sqrt((i-l)^2+(j-m)^2); %FASTEST
      D_i(c1,c2) = W_ii*r*exp(-d^2/sigma_i^2);
    end
  end
end


toc(CLOCK);
% save NN_dist_matrix D_e D_i
% error('asdf')

% initialize arrays-------------------------------------------
gs = zeros(N_cells, NT, 3); % all conductances
V = zeros(N_cells, NT); % all voltages
spikes = zeros(N_cells, NT); % raster plot
Vt = concat_EI(-73.6, -81.6); % current voltage initialized to val1 for E, val2 for I
I = zeros(N_cells,1); % input current [mV/s]
last_spike = zeros(N_cells, 1); % when was each cell's last spike?
g_tr = zeros(N_cells, 1); % after-hyperpol. conductance
g_e = zeros(N_cells, 1); % synaptic conductance from E cells
g_i = zeros(N_cells, 1); % synaptic conductance from I cells
P_post_e = zeros(N_cells, 1); % post-syn. channel open prob.
P_post_i = zeros(N_cells, 1); % post-syn. channel open prob.
P_rel_e = P0 * ones(N_cells, 1); % synaptic release prob.
P_rel_i = P0 * ones(N_cells, 1); % synaptic release prob.
y_e = zeros(N_E_cells, 1); % spiking at each time step
y_i = zeros(N_I_cells, 1); % spiking at each time step

% additional arrays used for Runge-Kutta
g_tr2 = zeros(N_cells, 1); % after-hyperpol. conductance
g_e2 = zeros(N_cells, 1); % synaptic conductance from E cells
g_i2 = zeros(N_cells, 1); % synaptic conductance from I cells
g_tr3 = zeros(N_cells, 1); % after-hyperpol. conductance
g_e3 = zeros(N_cells, 1); % synaptic conductance from E cells
g_i3 = zeros(N_cells, 1); % synaptic conductance from I cells
P_post_e2 = zeros(N_cells, 1); % post-syn. channel open prob.
P_post_i2 = zeros(N_cells, 1); % post-syn. channel open prob.
P_rel_e2 = P0 * ones(N_cells, 1); % synaptic release prob.
P_rel_i2 = P0 * ones(N_cells, 1); % synaptic release prob.

% arrays for GPU
D_e_gpu = gpuArray(D_e);
D_i_gpu = gpuArray(D_i);
y_e_gpu = gpuArray(y_e);
y_i_gpu = gpuArray(y_i);

% target input current to middle of grid
% NOTE: We also multiply I by 1000 so that I/C will be in units
% millivolt/sec, matching the units of V, rather than Volt/sec
arrayMap = reshape(1:N_E_cells, E_cell_dim);
mid1 = 0.5*E_cell_dim(1);
mid2 = 0.5*E_cell_dim(2);
I(reshape(arrayMap(mid1-1:mid1+1,mid2-1:mid2+1),1,[])) = ImVs;

% main time loop----------------------------------------------
fprintf('Computing voltages\n');
for t = 1:NT
  
  if mod(t,100)==0, fprintf(['t = ' num2str(t) '\n']), end
  
  % we compute values and plug into equations in reverse (Eqn 7 to Eqn 1)
  % for each conductance we evaluate at t, t+dt/2, t+dt for RK4 integration
  
  % first, compute new synaptic conductances g_s (Eqn 7 to Eqn 3)
  % compute S_s (Eqn 7)
%   S_e = D_e*y_e;
%   S_i = D_i*y_i;
  
  S_e = gather(D_e_gpu*y_e_gpu);
  S_i = gather(D_i_gpu*y_i_gpu);

  % compute alpha_s (Eqn 6)
  alpha_e = alpha_max_e * (1 - exp(-S_e ./ k_e));
  alpha_i = alpha_max_i * (1 - exp(-S_i ./ k_i));
  
  % compute P_post_s (Eqn 5)
  % probability of postsynaptic opening
  P_post_e2 = P_post_e + dt/2*(alpha_e .* (1-P_post_e) - beta_e * P_post_e);
  P_post_e = P_post_e + dt*(alpha_e .* (1-P_post_e) - beta_e * P_post_e);
  P_post_i2 = P_post_i + dt/2*(alpha_i .* (1-P_post_i) - beta_i * P_post_i);
  P_post_i = P_post_i + dt*(alpha_i .* (1-P_post_i) - beta_i * P_post_i);
  
  % compute P_rel_s (Eqn 4)
  % probability of presynaptic release
  P_rel_e2 = P_rel_e + dt/2/tau_rel_e*(P0 - P_rel_e);
  P_rel_e = P_rel_e + dt/tau_rel_e*(P0 - P_rel_e);
  P_rel_i2 = P_rel_i + dt/2/tau_rel_i*(P0 - P_rel_i);
  P_rel_i = P_rel_i + dt/tau_rel_i*(P0 - P_rel_i);
  
  % compute g_s (Eqn 3)
  g_e2 = g_e_max * P_post_e2 .* P_rel_e2;
  g_e3 = g_e_max * P_post_e .* P_rel_e;
  g_i2 = g_i_max * P_post_i2 .* P_rel_i2;
  g_i3 = g_i_max * P_post_i .* P_rel_i;
  
  % next, compute after-hyperpolarization conductances (Eqn 2)
  % first find which cells are refractory
  if time(t) > Trefract_e + T
    ghigh_ind = find(last_spike + T + refvec > time(t));
  else
    ghigh_ind = [];
  end
  g_ind = setdiff(1:N_cells,ghigh_ind); 
  g_tr(ghigh_ind) = g_tr_high;
  % for non-refractory cells, integrate g_tr
  g_tr2 = g_tr; g_tr3 = g_tr;
  g_tr2(g_ind) = g_tr(g_ind) - dt/2*g_tr(g_ind)./tau_tr(g_ind);
  g_tr3(g_ind) = g_tr(g_ind) - dt*g_tr(g_ind)./tau_tr(g_ind);
  
  % compute voltages (Eqn 1)
  % Fourth order Runge-Kutta method to integrate V
  % A's and B's come from rewriting dV/dt
  % and evaluating at t, t+dt/2, and t+dt, resp.
  % See Shelley & Tao, 2001 - Efficient and Accurate Time-Stepping
  % Schemes for Inte rate-and-Fire NNs
  A1 = 1./C.*(gL + g_e + g_i + g_tr);
  B1 = 1./C.*(I + EL.*gL + E_tr*g_tr + E_e*g_e + E_i*g_i);
  rk1 = -A1.*Vt+B1;
  
  A2 = 1./C.*(gL + g_e2 + g_i2 + g_tr2);
  B2 = 1./C.*(I + EL.*gL + E_tr*g_tr2 + E_e*g_e2 + E_i*g_i2);
  rk2 = -A2.*(Vt + rk1*dt/2)+B2;
  rk3 = -A2.*(Vt + rk2*dt/2)+B2;
  
  A3 = 1./C.*(gL + g_e3 + g_i3 + g_tr3);
  B3 = 1./C.*(I + EL.*gL + E_tr*g_tr3 + E_e*g_e3 + E_i*g_i3);
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
  
  % update g_s now that voltage is computed
  % CHECK INTEGRATION METHOD HERE
  g_e = g_e3;
  g_i = g_i3;
  g_tr = g_tr3;
%   g_e = (g_e+4*g_e2+g_e3)/6;
%   g_i = (g_i+4*g_i2+g_i3)/6;
%   g_tr = (g_tr+4*g_tr2+g_tr3)/6;
  
  % save voltage, spikes
  V(:,t) = Vt;
  spikes(last_spike==time(t), t) = 1;
  
  % local
%   y_e = y_e*0; y_i = y_i*0;
%   y_e(Vt(1:N_E_cells)==Vmax) = 1;
%   y_i(Vt(N_E_cells+1:end)==Vmax) = 1;
  
  % GPU
  y_e_gpu = y_e_gpu*0; y_i_gpu = y_i_gpu*0;
  y_e_gpu(Vt(1:N_E_cells)==Vmax) = 1;
  y_i_gpu(Vt(N_E_cells+1:end)==Vmax) = 1;
  
  % CHECK THIS??
  % presynaptic depression at spiking cells
  % release probability drops by factor fd_s < 1
  P_rel_e(spikes(:,t)==1) = fd_e * P_rel_e(spikes(:,t)==1);
  P_rel_i(spikes(:,t)==1) = fd_i * P_rel_i(spikes(:,t)==1);
  
  % temporary storage -- checking this
  gs(:,t,1) = g_e;
  gs(:,t,2) = g_i;
  gs(:,t,3) = g_tr;
  
end

toc(CLOCK);

end
