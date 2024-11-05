function [G, G_D_Cap, G_D_Dem, G_D_Cost, nLabels] = VSSPGRaph(F,R,demands)

% demands are a matrix of   ms1 c11 c12 ... c1n
%                           ms2 c21 c22 ... c2n
%                                   ...
%                           msm cm1 cm2 ... cmn

% flow is a matrix          m+1 x n+1
% capacity is a matrix      m+1 x n+1
% nodes are s = 1, F = [2..|F|+1], R = [|F|+2 ... |F|+|R|+1], t = N  


N           = 1+F+R+1;
G           = zeros(N,N);
G_D_Cost    = zeros(N,N);
G_D_Cap     = zeros(N,N,size(demands,1)*size(demands,2));
G_D_Dem     = zeros(N,N,size(demands,1)*size(demands,2));

nLabels = labelNodes(N, F, R); % node labels for ploting the graph

% demands at layer 1
G_D_Dem(:,1) = demands;
% demands at layer 3 
G_D_Dem(:,end) = -demands;

% edges layer 1 -> 2 
for j = 2:F+1
    G(1,j) = 1;
    G_D_Cap(1,j,:) = Inf;
    G_D_Cost(1,j) = 0;
end

for i = 1:R
    % edges from layer 2 -> 3
    for j = 2:F+1
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

%-----------------

function nLabels = labelNodes(N, F, R)
nLabels = cell(N); 
nLabels{1} = 's';
for i = 2:F+1
    nLabels{i} = char("F"+string(i)+"("+string(i)+")");
end
for i = F+2:F+R+1
    nLabels{i} = char("R"+string(i-(F+1))+"("+string(i)+")");
end
nLabels{N} = 't';
end