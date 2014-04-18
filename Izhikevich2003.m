% time parameters
dt = 1; % time resolution [ms]
NT = 1e3;
time = (1:NT)*dt;

% network parameters
N_cells = 1;

% injected current
I = 5;

% cell parameters
a = 0.02; b = 0.2; c = -65; d = 8; % regular-spiking (RS) cell
% a = 0.02; b = 0.2; c = -50; d = 2; % fast-repetitive-bursting (FRB) cell
% a = 0.1; b = 0.2; c = -65; d = 2; % fast-spiking (FS) cell
% a = 0.02; b = 0.25; c = -65; d = 2; % low-threshold-spiking (LTS) cell
a = 0.02; b = 0.25; c = -65; d = 0.05; % thalamo-cortical (TC) cell
% a = 0.02; b = 0.2; c = -55; d = 4; % intrinsically-bursting (IB) cell

% declare arrays for state variables
u = zeros(N_cells,NT);
V = zeros(N_cells,NT);

% initial conditions
V(:,1) = -65;
u(:,1) = b.*V(:,1);

% time loop
for t = 1:NT-1
  % compute voltages (Eqn 1)
  % Fourth order Runge-Kutta method to integrate V
  % A's and B's come from rewriting dV/dt
  % and evaluating at t, t+dt/2, and t+dt, resp.
  % See Shelley & Tao, 2001 - Efficient and Accurate Time-Stepping
  % Schemes for Integrate-and-Fire NNs

  fired = find(V(:,t) >= 30);
  if ~isempty(fired)
    V(fired,t) = c(fired);
    u(fired,t) = u(fired,t) + d(fired);
  end

  du = a.*(b.*V(:,t) - u(:,t));
  u1 = u(:,t);
  u2 = u1 + dt/2*du;
  u3 = u2;
  u4 = u1 + dt*du;
  
  V1 = V(:,t);
  k1 = 0.04*V1.^2 + 5*V1 + 140 - u1 + I;
  V2 = V1 + 0.5*k1;
  k2 = 0.04*V2.^2 + 5*V2 + 140 - u2 + I;
  V3 = V1 + 0.5*k2;
  k3 = .04*V3.^2 + 5*V3 + 140 - u3 + I;
  V4 = V1 + k3;
  k4 = 0.04*V4.^2 + 5*V4 + 140 - u4 + I;
  
  u(:,t+1) = u(:,t) + dt*du;
  V(:,t+1) = V(:,t) + dt/6*(k1+2*k2+2*k3+k4);
  
end

subplot(2,1,1); plot(time,V);
subplot(2,1,2); plot(time,u);