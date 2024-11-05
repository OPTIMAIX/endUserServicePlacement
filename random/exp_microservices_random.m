clear;
rng("default");

dateRun = datestr(now,'yyyymmddHHMM');

filename = sprintf('resultadosRandom_%s.csv', dateRun);
filenameLast = 'resultadosRandom_last.csv';
fileID = fopen(filename, 'a');

fprintf(fileID,'%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n', ...
    'Experiment', ...
    'ServiceTypes', ...
    'ElementTypes', ...
    'S',...
    'R',...
    'offerDemandRatio',...
    'minProfit', ...
    'maxProfit', ...
    'exitflag', ...
    'val', ...
    'elapsedTime');

numberOfMSperS = 3; % 3 identical microservices per service
demandPerStype = 10;
demandPerMStype = demandPerStype*numberOfMSperS;
offerDemandRatio = 0.9;
minProfit = 1;
maxProfit = 100;
infraDelays2pick = [10,20,30,40,50];

elapsedTime = 0;

for Exp = 1:10
    for ct = 3:3 %cloud types
        for st = 3:3 % service types,
            % clear everithing to start a new experiment
            clear cloudCapacityPerType cloudElementsPerType msReqsPerType msDemands microServiceTypeProfit F;
            fprintf('%s%i%s%i\n','Starting placement solver for st = ', st, ', ct = ', ct );

            % -------------------------------------------------------------
            % SERVICES ----------------------------------------------------
            % -------------------------------------------------------------
            
            % serviceTypes has the type of each service, each service has
            % numberOfMSperS microservices of the same type
            serviceTypes = zeros(st*demandPerStype,1);
            for i = 1:st
                for j = 1:demandPerStype
                    serviceTypes(i*demandPerStype-demandPerStype+j) = i; 
                end
            end

            % microServiceReqsPerType is a matrix with the requirements for ST, BW and
            % CPU for each one of the F microservice types
            msReqsPerType(1:st,1:3) = [(3:3:3*st)' (3:3:3*st)' (3:3:3*st)'];
            
            msDemands = demandPerMStype*ones(1,st);

            %The total demand of resources of the experiment
            totalDemand = sum(msReqsPerType.*msDemands','all');
            
            S = demandPerStype*st;         %number of services
            F = sum(msDemands);  %number of microservices
          
            b = zeros(S,numberOfMSperS);
            s = zeros(S,numberOfMSperS);
            c = zeros(S,numberOfMSperS);
            
            for i = 1:S
                for j = 1:numberOfMSperS
                    b(i,j) = msReqsPerType(serviceTypes(i),1);
                    s(i,j) = msReqsPerType(serviceTypes(i),2);
                    c(i,j) = msReqsPerType(serviceTypes(i),3);
                end
            end

            % delay restrictions between microservices of a service; they
            % are inf between microservices of different servces
            msDelayMatrix = matrixOfDelays(numberOfMSperS, infraDelays2pick);
            d = Inf(S,numberOfMSperS,S,numberOfMSperS); % FIXME: very sparce matrix
            
            for i = 1:S
                for j = 1:numberOfMSperS
                     for jp = 1:numberOfMSperS
                         d(i,j,i,jp) = msDelayMatrix(j,jp);
                     end
                end
            end
          
            % -------------------------------------------------------------
            % INFRASTRUCTURE ---------------------------------------------- 
            % -------------------------------------------------------------

            % cloudCapacityPerType is a matrix with the capacities of ST, BW and CPU for each type of
            % cloud infrastrucutre server
            cloudCapacityPerType(1:ct,1:3) = [(25:25:25*ct)' (25:25:25*ct)' (25:25:25*ct)'];

            % The number of available resouces are a ratio of the demands
            offer = totalDemand*offerDemandRatio;
            
            % For the previous offer, how many cloud servers of each type in
            % the system. FIXME: this only works because the resources are all
            % equivalent
            cloudElementsPerType = repmat(round(offer/sum(cloudCapacityPerType, 'all')),1,ct);

            % number of cloud elements R
            R = sum(cloudElementsPerType); 
 
            beta  = zeros(1,R);
            sigma = zeros(1,R);
            kapa  = zeros(1,R);
            cloudTypes = zeros(1,R);
            
            k = 1;
            for i = 1:ct
                for j = 1:cloudElementsPerType(i)
                    beta(k)  = cloudCapacityPerType(i,1);
                    sigma(k) = cloudCapacityPerType(i,2);
                    kapa(k)  = cloudCapacityPerType(i,3);
                    cloudTypes(k) = i; 
                    k = k+1;
                end
            end
            
            D = matrixOfDelays(R, infraDelays2pick);

            % -------------------------------------------------------------
            % PROFITS -----------------------------------------------------
            % -------------------------------------------------------------

            % microServiceTypeProfits represents the profit of allocating a service
            % of a given microServiceType on a server of a given cloudType.
            microServiceTypeProfit = profit(msReqsPerType, cloudCapacityPerType, minProfit, maxProfit);
            
            % p_ijk as in the paper, the profit of placing microservice j
            % of service i in the element k
            p = zeros(S,numberOfMSperS,R); 
            
            for i = 1:S
                for j = 1:numberOfMSperS
                    for k = 1:R
                        p(i,j,k) = microServiceTypeProfit(serviceTypes(i),cloudTypes(k));
                    end
                end
            end

            tic;

            % ACA YA ESTAN PREPARADAS TODAS LAS ENTRADAS 
            % SUSTITUIR LA SIGUIENTE LINEA POR EL SOLVER QUE CORRESPONDA
            [X,val] = randSolverMicroservices(b,s,c,beta,sigma,kapa,d,D,p); 

            elapsedTime = toc;

            fprintf(fileID,'%i, %i, %i, %i, %i, %f, %i, %i, %s, %f, %f\n', ...
                Exp, ...
                st, ...
                ct, ...
                S, ...
                R, ...
                offerDemandRatio, ...
                minProfit, ...
                maxProfit, ...
                " ", ...
                val, ...
                elapsedTime) ;
         
        end
    end
end

fclose(fileID);

copyfile(filename, filenameLast);

% -------------------------------------------------------------------------
% AUXILIAR FUNCTIONS ------------------------------------------------------
% -------------------------------------------------------------------------

function p = profit(serviceTypes,cloudTypes,minProfit, maxProfit)
p = repmat(log(sum(serviceTypes,2))/log(max(sum(serviceTypes,2))),1,size(cloudTypes,1)).*(log(sum(cloudTypes,2))/log(max(sum(cloudTypes,2))))';
p = floor(minProfit + (maxProfit-minProfit)*(p - min(p(:)))/((max(p(:) - min(p(:))))));
end
%

% asigns random delays between servers in a replicable manner
function m = matrixOfDelays(R, possible_values)
s = RandStream('mt19937ar', 'Seed', R);
num_possible_values = length(possible_values);
T = triu(possible_values(randi(s,num_possible_values, R, R)),1);
m = T + T';
end

function [X, val] = dummySolver(b,s,c,beta,sigma,kapa,d,D,p)
    X = zeros(size(p));
    val = 0;
end


%%%%%%%%%%%%%%
