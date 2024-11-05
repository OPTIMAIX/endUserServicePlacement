function [G, G_D_Cap, G_D_Dem, G_D_Cost, nLabels] = VSSPsimpleGRaph()

demands     = [1,2,3];
N           = 2;
G           = zeros(N,N);
G_D_Cost    = zeros(N,N);
G_D_Cap     = zeros(N,N,size(demands,1)*size(demands,2));
G_D_Dem     = zeros(N,N,size(demands,1)*size(demands,2));

nLabels = {'s','t'};

G_D_Dem(:,1) = demands;
G_D_Dem(:,end) = -demands;

G(1,2) = 1;
G_D_Cap(1,2,:) = [1,2,3];
G_D_Cost(1,2) = 0;


end


