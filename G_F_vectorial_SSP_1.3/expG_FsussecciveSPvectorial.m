%%
%
clear;
dateRun = datestr(now,'yyyymmddHHMM');


filename = sprintf('resultadosMinCostFlowSSPVect_%s.csv', dateRun);
fileID = fopen(filename, 'a');
        
fprintf(fileID,'%s, %s, %s, %s,%s, %s, %s, %s, %s, %s, %s, %s\n', ...
    'Exp', ...
    'demandPerServType', ...
    'offerDemandRatio', ...
    'minProfit', ...
    'maxProfit', ...
    'exitflag', ...
    'fval', ...
    'cloudTypes', ...
    'cloudElements', ...
    'serviceTypes', ...
    'serviceDemands', ...
    'elapsedTimeSSP');

demandPerServType = 100;
offerDemandRatio = 0.9;
minProfit = 1;
maxProfit = 100;
expectedProfit = 0; %FIXME
for Exp = 1:1
for R = 3:3 %cloud types
    for F = 3:3 % service types
        clear cloudTypes cloudElementsPerType serviceTypes serviceDemands serviceProfits;
        fprintf('%s%i%s%i\n','Starting vectorial ssp for F = ', F, ', R = ', R );
   
        serviceTypes(1:F,1:3) = [(3:3:3*F)' (3:3:3*F)' (3:3:3*F)'];
        cloudTypes(1:R,1:3) = [(25:25:25*R)' (25:25:25*R)' (25:25:25*R)'];
        serviceDemands = demandPerServType*ones(1,F);
        serviceProfits = profit(serviceTypes,cloudTypes, minProfit, maxProfit);
        demand = sum(serviceTypes.*serviceDemands','all');
        offer = demand*offerDemandRatio;
        cloudElementsPerType = repmat(round(offer/sum(cloudTypes, 'all')),1,R);
        flowDimensions = 3;
        
        tic
        [G, G_D_Cap, G_D_Demands, G_D_Cost, nLabels] = FlowNetConnGraphVectorial(serviceTypes,...
                                                                                    serviceProfits,...
                                                                                    cloudTypes,...
                                                                                    serviceDemands,...
                                                                                    cloudElementsPerType,...
                                                                                    flowDimensions);
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
        
        fprintf('%s%f%s%s%i%s%i\n','vectorial ssp computed in ', elapsedTimeSSP, 's', ', F = ', j, ', R = ', i );
        
        fprintf(fileID,'%i, %i, %f, %i, %i, %s, %f, %f, %f, %f, %f, %f\n', ...
                                                                            Exp, ...
                                                                            demandPerServType, ...
                                                                            offerDemandRatio, ...
                                                                            minProfit, ...
                                                                            maxProfit, ...
                                                                            " ", ...
                                                                            fval, ...
                                                                            size(cloudElementsPerType,2), ...
                                                                            sum(cloudElementsPerType), ...
                                                                            size(serviceTypes,1), ...
                                                                            sum(serviceDemands), ...
                                                                            elapsedTimeSSP) ;
        
    end
end
end


fclose(fileID);

%%
% 
function p = profit(serviceTypes,cloudTypes,minProfit, maxProfit)
p = repmat(log(sum(serviceTypes,2))/log(max(sum(serviceTypes,2))),1,size(cloudTypes,1)).*(log(sum(cloudTypes,2))/log(max(sum(cloudTypes,2))))';
p = floor(minProfit + (maxProfit-minProfit)*(p - min(p(:)))/((max(p(:) - min(p(:))))));
end
%
%%%%%%%%%%%%%%
