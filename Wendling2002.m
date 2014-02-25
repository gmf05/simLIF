function [ys,time] = Wendling2002()
  % physical parameters -----------------------------------------------
  % synaptic gains & time constants
  A = 3.25; % average E synaptic gain [mV]
  B = 22; % average slow I synaptic gain [mV]
  G = 10; % average fast I synaptic gain [mV]
  a = 100; % dendritic average time constant in the feedback E loop [Hz]
  b = 50; % dendritic average time constant in the slow feedback I loop [Hz]
  g = 500; % dendritic average time constant in the fast feedback I loop [Hz]
  % average numbers of synaptic contacts between the populations
  C = 135; % constant in terms of which other C's are defined
  C1 = C; % E feedback loop
  C2 = 0.8*C; % E feedback loop
  C3 = 0.25*C; % slow I feedback loop
  C4 = 0.25*C; % slow I feedback loop
  C5 = 0.3*C; % fast I feedback loop
  C6 = 0.1*C; % fast I feedback loop
  C7 = 0.8*C; % contacts between slow & fast I cells
  % parameters of the field -> spike transfer function:
  v0 = 6; % [mV]
  e0 = 2.5; % [Hz]
  r = 0.56; % [mV^-1]
  % transfer functions:
  S = @(v) 2*e0 ./ (1+exp(r*(v0-v)));
  
% %   % parameters for seizure (Figs 3 & 4)
% % %   A = 5;
% % %   B = 45; G = 5; % type 1
% % %   B = 35; G = 12.5; % type 2
% % %   B = 25; G = 5; % type 3
% % %   B = 10; G = 10; % type 4
% % %   B = 10; G = 25; % type 5
% % %   B = 15; G = 2; % type 6
  
  % time parameters---------------------------------------------------
  dt = 1e-4;
  NT = 1e4;
  time = (1:NT)*dt;
  p = normrnd(200,10,[1,NT]); % input current
  y = zeros(10, 1); % array for solving ODE's
  y0 = zeros(10, 1); % array for solving ODE's
  y1 = zeros(10, 1); % array for solving ODE's
  ys = zeros(10, NT); % array for solving ODE's
  
  % main time loop---------------------------------------------------
  for t = 1:NT
    
    y0(1:5) = y(1:5) + dt/2*(y(6:10));
    y0(6) = y(6) + dt/2*(A*a*S(y(2)-y(3)-y(4)) - 2*a*y(6) - a^2*y(1));
    y0(7) = y(7) + dt/2*(A*a*(p(t)-C2*S(C1*y(1))) - 2*a*y(7) - a^2*y(2));
    y0(8) = y(8) + dt/2*(B*b*C4*S(C3*y(1)) - 2*b*y(8) - b^2*y(3));
    y0(9) = y(9) + dt/2*(G*g*C7*S(C5*y(1) - C6*y(5)) - 2*g*y(9) - g^2*y(4));
    y0(10) = y(10) + dt/2*(B*b*S(C3*y(1)) - 2*b*y(10) - b^2*y(5));
    
    y1(1:5) = y(1:5) + dt*(y(6:10));
    y1(6) = y(6) + dt*(A*a*S(y(2)-y(3)-y(4)) - 2*a*y(6) - a^2*y(1));
    y1(7) = y(7) + dt*(A*a*(p(t)-C2*S(C1*y(1))) - 2*a*y(7) - a^2*y(2));
    y1(8) = y(8) + dt*(B*b*C4*S(C3*y(1)) - 2*b*y(8) - b^2*y(3));
    y1(9) = y(9) + dt*(G*g*C7*S(C5*y(1) - C6*y(5)) - 2*g*y(9) - g^2*y(4));
    y1(10) = y(10) + dt*(B*b*S(C3*y(1)) - 2*b*y(10) - b^2*y(5));            
    
    y = y1;
    ys(:,t) = y;
    
  end
  
end