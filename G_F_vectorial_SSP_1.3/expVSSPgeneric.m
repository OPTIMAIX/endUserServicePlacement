%%
%
clear;

clear cloudTypes cloudElementsPerType serviceTypes serviceDemands serviceProfits;

layer1nodes = 3;
layer2nodes = 3;
flowDimensions = 3;
demands 

[G, G_D_Cap, G_D_Demands, G_D_Cost, nLabels] = VSSPGRaph(layer1nodes,layer2nodes,flowDimensions, demands);
elapsedTimeGraph = toc;

s = digraph(double(G));
p = plot(s,'Layout','subspace3','MarkerSize',15,'Marker','o','EdgeFontSize',12,'NodeFontSize',12, 'NodeFontWeight', 'bold','NodeLabel',nLabels);
labeledge(p,s.Edges.EndNodes(:,1),s.Edges.EndNodes(:,2),string(G_D_Cost(sub2ind(size(G_D_Cost),s.Edges.EndNodes(:,1),s.Edges.EndNodes(:,2)))));

tic
sol = succShortestPathBFVectorial(G, G_D_Cap, G_D_Demands, G_D_Cost, serviceTypes);
elapsedTimeSSP = toc;

totAlocServices = sum(sol,3);
totAlocServices = sum(totAlocServices(F+1:F+R,end));
ceilProfit = max(serviceProfits,[],'all');
computedProfit = 0;
negProfit = 0;

for i = 1:F
    negProfit = negProfit + sum(sol(:,:,i).*G_D_Cost,'All');
end

fval = ceilProfit*totAlocServices - negProfit;

%%
%
function p = profit(serviceTypes,cloudTypes,minProfit, maxProfit)
p = repmat(log(sum(serviceTypes,2))/log(max(sum(serviceTypes,2))),1,size(cloudTypes,1)).*(log(sum(cloudTypes,2))/log(max(sum(cloudTypes,2))))';
p = floor(minProfit + (maxProfit-minProfit)*(p - min(p(:)))/((max(p(:) - min(p(:))))));
end
%
%%%%%%%%%%%%%%
