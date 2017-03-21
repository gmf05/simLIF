function sim_movie(V,params,Cax)
% Can attach Cax to params

if nargin<3
  Cax = [min(min(V)), max(max(V))];
%   Cax = [-75 -55]; % good for voltage
end

figure

dt = params.dt;
NT = params.NT; 
E_cell_dim = params.E_cell_dim;
I_cell_dim = params.I_cell_dim;

N_E_cells = prod(E_cell_dim);
time = (1:NT)*dt;
t_delay = .01; % pause used in plotting [seconds]



% % % Plot just E cells
% % for t = 1:NT
% %   imagesc(reshape(V(1:N_E_cells,t),E_cell_dim));
% %   caxis(Cax);
% %   xlabel([sprintf('%0.2g',time(t)) ' sec']);
% %   axis square;
% %   pause(t_delay);
% % end

% Plot E & I cells
for t = 1:NT
  subplot(2,1,1)
  imagesc(reshape(V(1:N_E_cells,t),E_cell_dim));
  caxis(Cax);
  xlabel([num2str(time(t)) ' sec']); 
  xlabel([sprintf('%0.2g',time(t)) ' sec']);
  ylabel('E cells'); axis normal;

  subplot(2,1,2)
  imagesc(reshape(V(N_E_cells+1:end,t),I_cell_dim));
  caxis(Cax); 
  ylabel('I cells'); axis normal;
  
  pause(t_delay);
end


end