%% Simulate data from LIF model

% Spatial effects curve
ext_knots = [0 5 10 20 50 100 125 150 175 200 300 400] * dt_ms;
ext_knots(1) = 1;
b_ext = 10*[5 5 2.5 1 0.5 0 0 0 0 0 0 0 0 0]';
% b_ext = 10*[5 5 2.5 1 0.5 0 0 1 2 1 0 0 0 0]';
% b_ext = 0*[5 5 3 0 0.3 1 0.1 0 0 0 0 0 0 0]';

% Intrinsic effects curve
% i_knots = [0 1 20 50 100 300 500 1000] * dt_ms;
int_knots = [0 5 20 50 100 125 150 175 200 300 400] * dt_ms;
% b_int = 1*[-5 -5 -1 0 0 0 0 0 0 0 0 0 0]';
b_int = 1*[-5 -5 -1 0 0 0 0.5 2 0.5 0 0 0 0]';
int_knots(1) = 1;

% Plot effects curves
g_int = cubic_spline(int_knots, b_int);
g_ext = cubic_spline(ext_knots, b_ext);
figure
subplot(211), plot(g_int);
subplot(212), plot(g_ext);

% Simulate LIF model
[V,dn,gs,time] = sim_LIF_GLM(ext_knots,b_ext,int_knots,b_int);
dt = time(2)-time(1);
Nchan = size(dn,1);
Ntime = length(time);

% % Simulate random (Poisson) data
% Nchan = 2; dt = 1e-4; Ntime = round(Nsec/dt);
% dn = poissrnd(10*dt, [Nchan, Ntime]);
% time = (1:Ntime)*dt;

% Create data object
d = pp_data(dn, time);

% Plot
figure, d.plot('raster');

%% Fit model
% Set parameters
% Choose knots and basis functions
dt_ms = round(.001 / dt);
T_knots = [0 1]; T_basis = 'indicator';

Q_knots = [0 50 100 200 500] * dt_ms;
% Q_knots = [0 100 200 500] * dt_ms;
Q_knots(1) = 1; % avoid using 0 lag
Q_basis = 'spline';

R_knots = [0 50 150 300 500]  * dt_ms;
% R_knots = [0 50 500] * dt_ms;
% R_knots = [0 50 150 250 500]  * dt_ms;
R_basis = 'spline';
R_knots(1) = 1; % avoid using 0 lag

% Create parameters object
p = pp_params();
response = Nchan;
neighbors = [Nchan-4 : Nchan-1];
p.response = response;
p = p.add_covar('rate', 0, T_knots, T_basis); 
p0 = p;
p = p.add_covar('intrinsic', response, Q_knots, Q_basis);
p = p.add_covar('spatial1',  neighbors(1), R_knots, R_basis);
p = p.add_covar('spatial2', neighbors(2), R_knots, R_basis);
p = p.add_covar('spatial3', neighbors(3), R_knots, R_basis);
% p = p.add_covar('spatial4', neighbors(4), R_knots, R_basis);

p1 = p0;
p1 = p1.add_covar('intrinsic', response, Q_knots, Q_basis);
p2 = p0;
p2 = p2.add_covar('spatial1',  neighbors(1), R_knots, R_basis);
p2 = p2.add_covar('spatial2',  neighbors(2), R_knots, R_basis);
p2 = p2.add_covar('spatial3',  neighbors(3), R_knots, R_basis);

% Estimate parameters
m = pp_model();
m0 = m.fit(d,p0);
m1 = m.fit(d,p1);
m2 = m.fit(d,p2);
m3 = m.fit(d,p);

% Deviance Table
dev_01 = m0.dev - m1.dev;
dev_02 = m0.dev - m2.dev;
dev_13 = m1.dev - m3.dev;
dev_23 = m2.dev - m3.dev;

delta_intrinsic = dev_23 ./ dev_01
delta_extrinsic = dev_13 ./ dev_02

% Np1 = length(p.covariate_ind{2});
% Np2 = length([p.covariate_ind{3:end}]);
% critical_f1 = finv(0.95,Np1,Np1);
% mn_f1 = Np1/(Np1-2);
% critical_f2 = finv(0.95,Np2,Np2);
% mn_f2 = Np2/(Np2-2);

% Hypothesis: delta ratios = 1
% Alternative: delta ratios != 1
% pval1 = 1-fcdf(delta_intrinsic, Np1, Np1);
% pval2 = 1-fcdf(delta_extrinsic, Np2, Np2);

% Plot resulting model
figure, m3.plot(d,p);