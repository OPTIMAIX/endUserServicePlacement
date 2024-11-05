function [G, G_D_Cap, G_D_Dem, G_D_Cost, nLabels] = FlowNetConnGraphVectorial(servTypes, servTypeProfits, cloudTypes, servDemands, cloudElemPerType, flowDimensions)

% % serviceType is a matrix with a row per type and a column per requirement,
% % in this case, col1 is bandwidth, col2 storage, and col3 computing.
% serviceTypes = [5,5,5; 10,10,10; 15,15,15];
%
% %infrastructure elements are of a given type. this matrix has a a row per
% %type ans a column per resource disponibility.
% cloudTypes = [50,50,50;100,100,100; 150,150,150; 200,200,200; 250,250,250];
% %cloudTypes = [15,15,15];
%
% %this matrix contains the number of service demands per service type
% serviceDemands = [500,200,300];
%
% %this matrix contains the number of cloud elements per type at least 1
% %per type
% %cloudElementsPerType = [250,200,150,100,50];
% cloudElementsPerType = [1,1,1,1,1];

R = size(cloudTypes,1);
F = size(servTypes,1);
N = 1+F+R+1;
ceilProfit = max(servTypeProfits, [],'all'); %+1;

G        = zeros(N,N);
G_D_Cost = zeros(N,N);
G_D_Cap  = zeros(N,N, flowDimensions); %TODO this must be a generic 3D matrix of NxNxM, M the number of dimensions of the flow. 
G_D_Dem  = zeros(F,N);

% nodes are s = 1, F = [2..|F|+1], R = [|F|+2 ... |F|+|R|+1], t = N  
nLabels{1} = 's';
for i = 2:F+1
    nLabels{i} = char("F"+string(i)+"("+string(i)+")");
end
for i = F+2:F+R+1
    nLabels{i} = char("R"+string(i-(F+1))+"("+string(i)+")");
end
nLabels{N} = 't';

% demands at layer 1
G_D_Dem(:,1) = servDemands;
% demands at layer 3 
G_D_Dem(:,end) = -servDemands;

% edges layer 1 -> 2 
for j = 2:F+1
    G(1,j) = 1;
    G_D_Cap(1,j,:) = Inf;
    G_D_Cost(1,j) = 0;
end

for i = 1:R
    for j = 2:F+1
        % edges from layer 2 -> 3
        G(j,F+1+i) = 1;
        G_D_Cap(j,F+1+i,:) = Inf;
        G_D_Cost(j,F+1+i) = ceilProfit - servTypeProfits(j-1,i); %transforms costs to make it a minimization problem. 
    end
    % edges from layer 3 -> 4
    G(F+1+i, N) = 1;
    G_D_Cap(F+1+i,N,:) = cloudElemPerType(i)*cloudTypes(i,:);
    G_D_Cost(F+1+i,N) = 0;
end

end

