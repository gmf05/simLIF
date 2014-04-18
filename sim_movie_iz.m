function sim_movie_iz(V,Cax,params)

count = 1;
filename = [date() '_' num2str(count) '.avi'];
while exist(filename,'file')
  count = count+1;
  filename = [date() '_' num2str(count) '.avi'];
end

% VW = VideoWriter(filename); 
% VW.FrameRate = 30;
% VW.open();

if nargin<2
  Cax = [min(min(V)), max(max(V))];
%   Cax = [-75 -55]; % good for voltage
%   Cax = [-80 -30]; % good for voltage
end

global dt
global NT
global N_cell_types

cell_dim = zeros(N_cell_types,2);
for i = 1:N_cell_types, cell_dim(i,:) = params{i}.dim; end
N_cells1 = prod(cell_dim(1,:));
time = (1:NT)*dt*1e-3;
t_delay = .001; % pause used in plotting [seconds]

% define handles sp1, sp2 which will be subplots?

for t = 1:NT
  
  subplot(N_cell_types,1,1)
  imagesc(reshape(V(1:N_cells1,t),cell_dim(1,:)));
  caxis(Cax);
  title(num2str(time(t)));

  subplot(N_cell_types,1,2)
  imagesc(reshape(V(N_cells1+1:end,t),cell_dim(2,:)));
  caxis(Cax);
  
  pause(t_delay);
%   VW.writeVideo(getframe());
end

% VW.close();

end