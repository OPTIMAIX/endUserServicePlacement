clear;

dateRun = datestr(now,'yyyymmddHHMM');

numberOfMSperS = 3; % 3 identical microservices per service
demandPerStype = 10;
demandPerMStype = demandPerStype*numberOfMSperS;
offerDemandRatio = 0.9;
minProfit = 1;
maxProfit = 100;
infraDelays2pick = [10,20,30,40,50];

elapsedTime = 0;

%have to produce:
%

for Exp = 1:1
    for ct = 3:3 %cloud types
        for st = 3:3 % service types,
            % clear everithing to start a new experiment
            clear cloudCapacityPerType cloudElementsPerType msReqsPerType msDemands microServiceTypeProfit F;
            fprintf('%s%i%s%i\n','Starting data creation for MILP solver for st = ', st, ', ct = ', ct );

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

            F = repelem((1:1:numberOfMSperS),S,1);

            b = zeros(S,numberOfMSperS);
            s = zeros(S,numberOfMSperS);
            c = zeros(S,numberOfMSperS);

            for i = 1:S
                for j = 1:numberOfMSperS
                    s(i,j) = msReqsPerType(serviceTypes(i),2);
                    b(i,j) = msReqsPerType(serviceTypes(i),1);
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

            %--------------------------------------------------------------
            % genero archivo de datos
            %--------------------------------------------------------------

            instance = sprintf('%i_%i_%i', Exp, ct, st);
            filename = sprintf('input_%s.json', instance);
            fileID = fopen(filename, 'w');
            fprintf(fileID,'{\n');
            fprintf(fileID,'    "comment": "EDGE problem instance",\n');
            fprintf(fileID,'    "comment": "Version 1 %s Exp=%i, ct=%i, st=%i",\n', dateRun, Exp, ct, st);
            fprintf(fileID,'    "comment": "instname := d_%s",\n', instance);
            % fprintf(fileID,'\n');
            % fprintf(fileID,'param n := %i;\n', S);

            % services 
            fprintf(fileID,'    "services": [\n');
            for i = 1:S
                fprintf(fileID,'      {\n');
                fprintf(fileID,'        "name":"serv%i",\n', i);
                fprintf(fileID,'        "microservices": [\n');
                for j = 1:numberOfMSperS
                    fprintf(fileID,'          {\n');
                    fprintf(fileID,'            "name":"microServ%i",\n', F(i,j));
                    fprintf(fileID,'            "s": %i,\n', s(i,j));
                    fprintf(fileID,'            "b": %i,\n', b(i,j));
                    fprintf(fileID,'            "c": %i\n', c(i,j));
                    if j == numberOfMSperS
                        fprintf(fileID,'          }\n');
                    else
                        fprintf(fileID,'          },\n');
                    end
                end

                fprintf(fileID,'        ],\n');
                fprintf(fileID,'        "max_delay": [\n');
                for j = 1:numberOfMSperS-2
                    for k = j+1:numberOfMSperS
                        fprintf(fileID,'          ["microServ%i", "microServ%i", %i],\n', j, k, d(i,j,i,k));
                    end
                end
                fprintf(fileID,'          ["microServ%i", "microServ%i", %i]\n', numberOfMSperS-1, numberOfMSperS, d(i,numberOfMSperS-1,i,numberOfMSperS));
                
                fprintf(fileID,'        ]\n');
                if i == S
                    fprintf(fileID,'      }\n');
                else
                    fprintf(fileID,'      },\n');
                end
            end
            fprintf(fileID,'    ],\n');

            %network
            fprintf(fileID,'    "network": [\n');
            fprintf(fileID,'      {\n');
            fprintf(fileID,'        "elements": [\n');

            for i = 1:R
                fprintf(fileID,'          {\n');
                fprintf(fileID,'            "name":"inf%i",\n', i);
                fprintf(fileID,'            "s_max": %i,\n', sigma(i));
                fprintf(fileID,'            "b_max": %i,\n', beta(i));
                fprintf(fileID,'            "c_max": %i\n', kapa(i));
                if i == R
                    fprintf(fileID,'          }\n');
                else
                    fprintf(fileID,'          },\n');
                end
            end
            fprintf(fileID,'        ],\n');
            fprintf(fileID,'        "delays": [\n');
            for j = 1:R
                for k = 1:R
                    if (j == R && k == R)
                        fprintf(fileID,'          ["inf%i", "inf%i", %i]\n', j, k, D(j,k));
                    else
                        fprintf(fileID,'          ["inf%i", "inf%i", %i],\n', j, k, D(j,k));
                    end
                end
            end
            
            fprintf(fileID,'        ]\n');
            fprintf(fileID,'      }\n');
            fprintf(fileID,'    ]\n');
            fprintf(fileID,'  }\n');

            fprintf(fileID,'\n');
            
            % for i = 1:S
            %     for j = 1:numberOfMSperS
            %         for k = 1:R
            %             fprintf(fileID,' %i %i %i  %i\n',i,j,k,p(i,j,k));
            %         end
            %     end
            % end
            fprintf(fileID,'\n');

            fclose(fileID);
        end
    end
end



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
