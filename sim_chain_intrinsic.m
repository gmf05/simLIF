function [V,spikes,gs,time] = sim_chain_intrinsic()
%{

Enhanced LIF model, arranged in a 1-dim chain.

Model:
  dV[i]/dt = noise + I + -L*V[i] + gs[i,t_spike1] + gi[i,t_spike2]

noise is Gaussian noise.

I is tonic input. (Only supplied to Cell 1 in the chain here)

L is leak factor. How quickly does cell return to rest (V=0)?

gs and gi are functions of the spike times in channel i / its neighbors.

gs defines how cell i integrates synaptic input from its neighbors.

gi defines how cell i varies based on its last spike. We can use it to make refractoriness and/or rebound. Resetting, too?

%}

CLOCK = tic();

% Set parameters
Nsec = 10;
Ncells = 30;
L = 1.0;
I0 = 5;
Vth = 1; % Voltage threshold for spiking
targets = [2:Ncells 0];
dt = 1e-4;
dt_ms = 1e-3/dt;
dt_sec = 1/dt;
Ntime = Nsec/dt;
time = (1:Ntime) * dt;
noise_std = 80;
noise = normrnd(0, noise_std, [Ncells,Ntime]);

% Arrays to hold simulated data
gs = zeros(Ncells, Ntime, 2); % all conductances
V = zeros(Ncells, Ntime); % all voltages
spikes = zeros(Ncells, Ntime); % raster plot
Vt = 0.2 * rand(Ncells,1);
I = zeros(Ncells,1); % input current
I(1) = I0; % Inject current into one end of chain
last_spike = zeros(Ncells, 1); % when was each cell's last spike?

% Spatial effects curve
s_knots = [0 5 10 20 50 100 200] * dt_ms;
s_knots(1) = 1;
max_s = s_knots(end);
% b_s = [2 2 1 0.2 0 0]';
b_s = 10*[20 20 10 2 0 0 0 0 0]';
% b_s = 10*[20 20 10 2 1 1 1 0 0]';
g_s = cubic_spline(s_knots, b_s);
gs_t = zeros(Ncells,1);

% Intrinsic effects curve
i_knots = [0 1 20 50 100 300 500 1000] * dt_ms;
i_knots(1) = 1;
max_i = i_knots(end);
% b_i = [-5 -5 -3 0 0 0 0 0 0 0]';
b_i = 1*[-10 -10 -5 -1 2 0.3 0 0 0 0]';
g_i = cubic_spline(i_knots, b_i);
gi_t = zeros(Ncells,1);

% main time loop----------------------------------------------
fprintf('Computing voltages\n');
for t = 1:Ntime
  
  if mod(t,dt_sec)==0, fprintf(['t = ' num2str(t/dt_sec) '\n']), end
  
  % update "conducatances" based on time since last spike
  idx1 = find(last_spike>0 & t-last_spike<max_i);  
  gi_t(idx1) = g_i(t-last_spike(idx1));
  
  % update "conducatance" based on time since last spike
  idx2 = find(last_spike>0 & t-last_spike<max_s);
  idx2(targets(idx2)==0)=[];
  gs_t(targets(idx2)) = g_s(t-last_spike(idx2));
  
  % save conductances
  gs(:,t,1) = gi_t;
  gs(:,t,2) = gs_t;
  
  % increment voltage
  Vt = Vt + dt*( -L*Vt + I + noise(:,t) + gs_t + gi_t );
  
  % update last_spike
  spiking = find(Vt >= Vth);
  spikes(spiking, t) = 1;
  last_spike(spiking) = t;
  
  % save voltage, reset spiking cells
  V(:,t) = Vt;
  Vt(spiking)=0;
    
end

toc(CLOCK);

end
