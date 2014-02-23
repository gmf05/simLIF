function get_sim_params()

% declare global variables
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

% time parameters
dt = 1e-4; % time step [s]
NT = 8e3; % length of sim [tauime steps]

% network size, topology
E_cell_dim = [50 50];
I_cell_dim = [25 25];
network_topology = 'sheet';

% PHYSICAL PARAMETERS:
% parameter | description [units] 

% input current
I0 = 6; % input current [nA]
ImVs = I0*1000; % adjusted current so dV/dt is in millivolts/sec

% capacitance & reversal potentials---------------------------
% See Table 1, Hall & Kuhlmann 2013
C_e = 0.5; % membrane capacitance for E cells [nF] <<<CHECK
C_i = 0.2; % membrane capacitance for I cells [nF] <<<CHECK
gL_e = 20; % leak conductance for E cells [nS] <<<CHECK
gL_i = 25; % leak conductance for I cells [nS] <<<CHECK
EL_e = -73.6; % leak reversal potential for E cells [mV]
EL_i = -81.6; % leak reversal potential for I cells [mV]
g_e_max = 80; % synaptic conductance with all channels open [nS] <<<CHECK
g_i_max = 120; % % synaptic conductance with all channels open [nS] <<<CHECK
E_e = 0; % synaptic reversal potential (E cells) [mV]
E_i = -70; % synaptic reversal potential (I cells) [mV]
% after-hyperpolarization parameters:
% g_tr_high_e = 60; % after-hyperpolarization conducatance during ref. [nS] <<<CHECK
% g_tr_high_i = 60; % after-hyperpolarization conductance during ref. [nS] <<<CHECK
g_tr_high = 60; % after-hyperpolarization conducatance during ref. [nS] <<<CHECK
tau_tr_e = 15e-3; % time constant for after-hyperpol. (spike rate adapt.) [s]
tau_tr_i = 5e-3; % time constant for after-hyperpol. (spike rate adapt.) [s]
E_tr = -90; % reversal potential for after-hypoerpol. [mV]

% membrane potential spiking/resetting params-------------------
Vmax = 0; % where membrane potential is set if > threshold [mV]
Vth = -55; % threshold to declare spike [mV]
V0 = -70; % where membrane potential is set during refractory period [mV]
T = 1e-3; % duration of a spike [s]
Trefract_e = 3e-3; % refractory peroid for E cells [s]
Trefract_i = 1.6e-3; % refractory period for I cells [s]

% synaptic kinetics parameters---------------------------------
beta_e = 0.667; % [ms^-1] <<<CHECK
beta_i = 0.19; % [ms^-1] <<<CHECK
alpha_max_e = 2.667; % [ms^-1] <<<CHECK
alpha_max_i = 3.143; % [ms^-1] <<<CHECK
k_e = 5; % saturation rate of synaptic conductance for E cells
k_i = 20; % saturation rate of synaptic conductance for I cells
P0 = 1; % resting prob. of synaptic release

% synaptic connectivity parameters ----------------------------
% sigma_e = 2; % width of E cell kernel
% sigma_i = 3.75; % width of I cell kernel (3.75 or 1)
% W_ee = 1.5; % synaptic strength E->E
% W_ei = 1.5; % synaptic strength E->I (1.5 or 4)
% W_ie = 4.5; % synaptic strength I->E (4.5 or 14)
% W_ii = 4.5; % synaptic strength I->I (4.5 or 1.5)
sigma_e = 4; % width of E cell kernel
sigma_i = 1; % width of I cell kernel (3.75 or 1)
W_ee = 1.5; % synaptic strength E->E
W_ei = 4; % synaptic strength E->I (1.5 or 4)
W_ie = 14; % synaptic strength I->E (4.5 or 14)
W_ii = 1.5; % synaptic strength I->I (4.5 or 1.5)

% W_ee = 2.5; W_ei = 0; W_ie = 0; W_ii = 0;

% presynaptic depression parameters---------------------------
tau_rel_e = 2e-1; % [s]
tau_rel_i = 4e-1; % [s]
fd_e = 0.9997;
fd_i = 0.9995;

end

