function output = NFPlanningSolverStructs(json_data)

% Decode JSON string to a MATLAB struct
data = jsondecode(json_data);

services = dict_to_cell_array(data.services);
network = dict_to_cell_array(data.network);

% Services ----------------------------------------------------------------

S = length(services);
numberOfMSperS = 0;

for i = 1:S
    microservices = services{i}.microservices;
    numberOfMSperS = max(length(microservices),numberOfMSperS);
end

b = zeros(S,numberOfMSperS);
s = zeros(S,numberOfMSperS);
c = zeros(S,numberOfMSperS);

% delay restrictions between microservices of a service; they
% are inf between microservices of different servces
d = Inf(S,numberOfMSperS,S,numberOfMSperS);

for i = 1:S
    microservices = services{i}.microservices;
    numMicroservices = length(microservices);

    for j = 1:numMicroservices
        s(i,j) = microservices(j).s;
        b(i,j) = microservices(j).b;
        c(i,j) = microservices(j).c;
    end

    % Extract max delays
    max_delays = services{i}.max_delay;

    for j = 1:size(max_delays, 1)
        [~,from] = ismember(max_delays{j}{1},{microservices(:).name});
        [~,to] = ismember(max_delays{j}{2},{microservices(:).name});
        d(i,from,i,to) = max_delays{j}{3};
    end

end

% Network -----------------------------------------------------------------

R = length(network{1,1}.elements);

beta  = zeros(1,R);
sigma = zeros(1,R);
kapa  = zeros(1,R);

for i = 1:R
    beta(i) = network{1,1}.elements(i).b_max;
    sigma(i) = network{1,1}.elements(i).b_max;
    kapa(i) = network{1,1}.elements(i).b_max;
end

D = Inf(R,R);

for j = 1:size(network{1,1}.delays,1)
    [~,from]   = ismember(network{1,1}.delays{j}{1},{network{1,1}.elements(:).name});
    [~,to]     = ismember(network{1,1}.delays{j}{2},{network{1,1}.elements(:).name});
    D(from,to) = network{1,1}.delays{j}{3};
    D(to,from) = network{1,1}.delays{j}{3};
end

%FIXME: should add profit to json

% microServiceTypeProfits represents the profit of allocating a service
% of a given microServiceType on a server of a given cloudType.

% p_ijk as in the paper, the profit of placing microservice j
% of service i in the element k
p = ones(S,numberOfMSperS,R);

%calling the solver
[X,val] = randSolverMicroservices(b,s,c,beta,sigma,kapa,d,D,p);

structOutput = XtoStruct(X);

%output = jsonencode(structOutput);

output = '{   "comment": "EDGE problem planning, Version 1 202403311716 Exp=1, ct=5, st=5, instname := d_1_3_3",   "network":     {       "elements": [         {           "name": "inf1",           "microservices": [             {               "service_name":"serv1",               "microservice_name":"micrServ1"             },             {               "service_name":"serv1",               "microservice_name":"micrServ2"             },             {               "service_name":"serv2",               "microservice_name":"micrServ1"             }           ]         },         {           "name": "inf1",           "microservices": [             {               "service_name":"serv3",               "microservice_name":"micrServ1"             },             {               "service_name":"serv3",               "microservice_name":"micrServ2"             }           ]         }        ]     } }';

end


%--------------------------------------------------------------------------
% Aux Functions
%--------------------------------------------------------------------------

function cell_array = dict_to_cell_array(dict_list)
% Converts a list of structures to a cell array
cell_array = {};
for i = 1:length(dict_list)
    item = dict_list(i);
    fields = fieldnames(item);
    cell_item = struct();
    for j = 1:length(fields)
        field = fields{j};
        value = item.(field);
        if isnumeric(value) || ischar(value)
            cell_item.(field) = {value};
        else
            cell_item.(field) = value;
        end
    end
    cell_array{end+1} = cell_item;
end
end

function p = profit(serviceTypes,cloudTypes,minProfit, maxProfit)
p = repmat(log(sum(serviceTypes,2))/log(max(sum(serviceTypes,2))),1,size(cloudTypes,1)).*(log(sum(cloudTypes,2))/log(max(sum(cloudTypes,2))))';
p = floor(minProfit + (maxProfit-minProfit)*(p - min(p(:)))/((max(p(:) - min(p(:))))));
end

function data = XtoStruct(X)

data = struct;
data.comment = 'EDGE problem planning';
data.comment_1 = 'Version 1 202403311716 Exp=1, ct=5, st=5';
data.comment_2 = 'instname := d_1_3_3';

data.network = struct;
data.network.elements = struct;

% X(service, microservice, device)

for k = 1:size(X,3)
    data.network.elements(k, 1).name = 'inf1';
    data.network.elements(k, 1).microservices = struct;
    for i = 1:size(X,1)
        for j = 1:size(X,2)
            if X(i,j,k)
                data.network.elements(k, 1).microservices(i, 1).service_name = 'servI';
                data.network.elements(k, 1).microservices(i, 1).microservice_name = 'micrServJ';
            end
        end
    end
end

end
