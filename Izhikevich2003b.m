%----- declare global parameters
global dt
global NT
global params
global N_cells
global N_cell_types

params{1}.name = 'RS';
params{1}.dim = [20 20];
params{1}.N_cells = prod(params{1}.dim);
params{1}.a = [0.01 0.03];
params{1}.b = [0.2 0.2];
params{1}.c = [-65 -65];
params{1}.d = [8 8];

% params{1}.name = 'CH';
% params{1}.dim = [20 20];
% params{1}.N_cells = prod(params{1}.dim);
% params{1}.a = [0.01 0.03];
% params{1}.b = [0.2 0.2];
% params{1}.c = [-50 -50];
% params{1}.d = [2 2];

params{2}.name = 'FS';
params{2}.dim = [10 10];
params{2}.N_cells = prod(params{2}.dim);
params{2}.a = [0.05 0.15];
params{2}.b = [0.2 0.2];
params{2}.c = [-65 -65];
params{2}.d = [2 2];

% INTERESTING CONNECTIVITY PARAMETERS:
% -------------------------------------
%
% (1) results in outward waves w/out reverb
% params{1}.weights = [30 10]; params{1}.sigma = 1; params{2}.weights = [10 10]; params{2}.sigma = 1; 
% 
% (2) results in outward waves w/ some flickering
% params{1}.weights = [50 10]; params{1}.sigma = 1; params{2}.weights = [10 10]; params{2}.sigma = 1; 
% 
% (3) lots of flickering
% params{1}.weights = [60 10]; params{1}.sigma = 1; params{2}.weights = [10 10]; params{2}.sigma = 1; 
% 
% (4) multi-directional waves
% params{1}.weights = [60 10]; params{1}.sigma = 1; params{2}.weights = [30 10]; params{2}.sigma = 1; 
% 
% (5) multiple sources develop
% params{1}.weights = [40 10]; params{1}.sigma = 1; params{2}.weights = [10 30]; params{2}.sigma = 1;
% 
% (6) outward wave w/ reverb
params{1}.weights = [5 10]; params{1}.sigma = 4; params{2}.weights = [10 10]; params{2}.sigma = 1;
% params{1}.weights = [5 15]; params{1}.sigma = 4; params{2}.weights = [10 10]; params{2}.sigma = 1;
% params{1}.weights = [5 10]; params{1}.sigma = 4; params{2}.weights = [15 10]; params{2}.sigma = 1;
% params{1}.weights = [5 10]; params{1}.sigma = 4; params{2}.weights = [10 15]; params{2}.sigma = 1;

% params{1}.weights = [5 5]; params{1}.sigma = 4; params{2}.weights = [5 5]; params{2}.sigma = 1;

%----- time parameters
dt = 1; % time resolution [ms]
NT = 1e3; % number of time steps
time = (1:NT)*dt; % time axis

%------ cell parameters
N_cell_types = length(params);
N_cells = 0;
for i = 1:N_cell_types, N_cells = N_cells + params{i}.N_cells; end;

a = zeros(N_cells,1);
b = zeros(N_cells,1);
c = zeros(N_cells,1);
d = zeros(N_cells,1);

count = 0;
for i = 1:N_cell_types
  C = params{i}.N_cells;
  a(count+(1:C)) = unifrnd(params{i}.a(1),params{i}.a(2),[C,1]);
  b(count+(1:C)) = unifrnd(params{i}.b(1),params{i}.b(2),[C,1]);
  c(count+(1:C)) = unifrnd(params{i}.c(1),params{i}.c(2),[C,1]);
  d(count+(1:C)) = unifrnd(params{i}.d(1),params{i}.d(2),[C,1]);
  count = count + C;
end

%----- injected current to middle E cells
I0 = zeros(N_cells,1);
N_E_cells = prod(params{1}.dim);
temp = reshape(1:N_E_cells,params{1}.dim(1),params{1}.dim(2));
md = params{1}.dim(1)/2 + [-1:1];
input_cells = reshape(temp(md,md),1,[]);
I0(input_cells) = 4; % bifurcation from 3 to 4


%----- inter-cell geometry
S = zeros(N_cells, N_cells); % connectivity matrix
rho = unifrnd(0.5,1.5,[1,N_cells]); % random element to connectivity

% make lists mapping cell # to cell type #
cell_types = zeros(1,N_cells);
count = 0;
for i = 1:N_cell_types
  C = params{i}.N_cells;
  cell_types(count+(1:C)) = i;
  count = count + C;
end

%%------------------------------------------------------------------------
% % % NOTE1: Needs to be revised !!!!!
% % % NOTE2: Assumes for now that there are more E cells than I cells
% % % and remaps coord of I cells to approximate pos. relative to E cells
% % % if there are more I cells, this mapping should be reversed
E_cell_dim = params{1}.dim;
I_cell_dim = params{2}.dim;

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
%%------------------------------------------------------------------------


for i = 1:N_cells
  m = cell_types(i);
  r = rho(i);
  x1 = coord(i,1); y1 = coord(i,2);
  for j = 1:N_cells
    n = cell_types(j);
    sigma = params{m}.sigma;
    W = params{m}.weights(n);
    x2 = coord(j,1); y2 = coord(j,2);
    d0 = sqrt((x1-x2)^2+(y1-y2)^2);
    S(i,j) = W*r*exp(-d0^2/sigma^2);
  end
end

%----- declare arrays for state variables
u = zeros(N_cells,NT);
V = zeros(N_cells,NT);
dn = zeros(N_cells,NT);

%----- initial conditions
V(:,1) = -65;
u(:,1) = b.*V(:,1);

%----- time loop
% Fourth order Runge-Kutta integration
% See Shelley & Tao, 2001 - Efficient and Accurate Time-Stepping
% Schemes for Integrate-and-Fire NNs
for t = 1:NT-1

  % reset, save spiking cells 
  fired = find(V(:,t) >= 30);
  if ~isempty(fired)
    V(fired,t) = c(fired);
    u(fired,t) = u(fired,t) + d(fired);
  end
  dn(fired,t) = 1;
  
  % sum synaptic input from spiking cells
  if t<60
    I = I0;
  else
    I = 0;
  end
   I = I + sum(S(:,fired),2);
%   I = I0 + sum(S(:,fired),2);

  % increment recovery variable u
  du = a.*(b.*V(:,t) - u(:,t));
  u1 = u(:,t);
  u2 = u1 + dt/2*du;
  u3 = u2;
  u4 = u1 + dt*du;
  
  % increment voltage V (RK4)
  V1 = V(:,t);
  k1 = 0.04*V1.^2 + 5*V1 + 140 - u1 + I;
  V2 = V1 + 0.5*k1;
  k2 = 0.04*V2.^2 + 5*V2 + 140 - u2 + I;
  V3 = V1 + 0.5*k2;
  k3 = .04*V3.^2 + 5*V3 + 140 - u3 + I;
  V4 = V1 + k3;
  k4 = 0.04*V4.^2 + 5*V4 + 140 - u4 + I;
  
  % save state variables (u, V)
  u(:,t+1) = u(:,t) + dt*du;
  V(:,t+1) = V(:,t) + dt/6*(k1+2*k2+2*k3+k4);
  
end

% %%
%----- plot results
% subplot(2,1,1); plot(time(1:end-1),V(:,1:end-1));
% subplot(2,1,2); plot(time(1:end-1),u(:,1:end-1));
% %%
%----- results as movie
sim_movie_iz(V,[-80,-30],params);