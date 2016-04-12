function [V,spikes,g,time] = sim_LIF_GLM(ext_knots,b_ext,int_knots,b_int)
%{

Enhanced LIF model

Model:
  dV[i]/dt = noise + I + -L*V[i] + g_ext[t - t_spike_j] + g_int[t - t_spike_i]

noise is Gaussian noise.

I is tonic input.

L is leak rate.

g_ext and g_int are functions of the spike times in channel i / its neighbors.
  (t_spike_i is the time of cell i's last spike. ditto for t_spike_j.)

g_ext defines how cell i integrates synaptic input from its neighbors.
  If there are multiple neighbors per cell (almost always), then this term
  will be a sum over multiple j

g_int defines how cell i varies based on its last spike. 
  We can use it to make refractoriness and/or rebound.

%}

CLOCK = tic();

% Set parameters
Nsec = 30;
Ncells = 20;
Vth = 1; % Voltage threshold for spiking

% Time parameters
% dt = 1e-4;
dt = 1e-3;
dt_ms = 1e-3/dt;
dt_sec = 1/dt;
Ntime = Nsec/dt;
time = (1:Ntime) * dt;

% Input drive + noise
I0 = 6;
I = I0*ones(Ncells,1); % tonic drive
noise_std = 60;
noise = normrnd(0, noise_std, [Ncells,Ntime]);

% Extrinsic effects curve
g_ext = cubic_spline(ext_knots, b_ext);
max_s = ext_knots(end); % duration of extrinsic effects

% Intrinsic effects curve
g_int = cubic_spline(int_knots, b_int);
max_i = int_knots(end); % duration of intrinsic effects

% Set network topology
N_feedforward = 3;
rho = unifrnd(0.5,1.5,Ncells);
sigma = 0.4;

% Connectivity matrix
W = zeros(Ncells);
for n = 1:Ncells
  
  % Which cells are connected?
  if n<=N_feedforward
    N = n-1;
  else
    N = N_feedforward;
  end
  neighbors = n-N:n-1;
  
  % Really simple uniform synapses:
  %W(n, neighbors) = 1;
  
  % Adding some distance-dependence, randomness to synapses
  %d = 1./(n-neighbors); % compute distance d
  d = abs(n-neighbors);
  W(n, neighbors) = rho(n, neighbors) .* exp(-(d ./ sigma).^2);
end

% Arrays to hold simulated data
g = zeros(Ncells, Ntime, 2); % all conductances
V = zeros(Ncells, Ntime); % all voltages
spikes = zeros(Ncells, Ntime); % raster plot
Vt = 1 * rand(Ncells,1); % initial voltages
last_spike = nan(Ncells, 1); % when was each cell's last spike?

% main time loop----------------------------------------------
fprintf('Computing voltages\n');
for t = 1:Ntime
  
%   if mod(t,dt_sec)==0, fprintf(['t = ' num2str(t/dt_sec) '\n']), end
  
  delays = t - last_spike;  

  % update "conducatances" based on time since last spike
  delays1 = delays; 
  delays1(delays1>max_i)=nan;
  g_int_t = nanzero(g_int, delays1);
  
  % synaptic conductances = weight matrix * conductance vector
  % conductance vector = g_s ( t - last_spike ) , if last_spike is recent
  % enough
  delays2 = delays;
  delays2(delays2>max_s) = nan;
  gs = nanzero(g_ext, delays2);
  g_ext_t = W * gs;

  % save conductances
  g(:,t,1) = g_int_t;
  g(:,t,2) = g_ext_t;
  
  % increment voltage
  Vt = Vt + dt*( noise(:,t) + I - Vt + g_ext_t + g_int_t );
  
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

function y = nanzero(x, i)
I = isnan(i);
i(I)=1;
y=x(i);
y(I)=0;
end
