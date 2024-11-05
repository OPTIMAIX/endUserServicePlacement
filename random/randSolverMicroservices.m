function [x, val] =  randSolverMicroservices(b,s,c,beta,sigma,kapa,d,D,p)

n = size(b,1);
m = size(b,2);
r = size(beta,2);
x = zeros(n,m,r);
val = 0;

for i = 1:n
    for j = 1:m
        valpre = val; 
        %find all suitable nodes in terms of resources
        suitable = suitableNodes(1:length(beta),b(i,j),s(i,j),c(i,j),beta,sigma,kapa);
        if ~isempty(suitable)
            k = suitable(randi(length(suitable)));
            x(i,j,k) = 1;
            val = val + p(i,j,k); 
            % extract the assigned resources
            beta(k) = beta(k)-b(i,j);
            sigma(k) = sigma(k)-s(i,j);
            kapa(k) = kapa(k)-c(i,j);
            for l = j+1:m
                %find all possible nodes in terms of delay
                possible = findCloseNodes(D, k, d(i,j,i,l));
                %from possible nodes chose those with resources
                suitable = suitableNodes(possible,b(i,l),s(i,l),c(i,l),beta,sigma,kapa);
                if ~isempty(suitable)
                    newk = suitable(randi(length(suitable)));
                    x(i,l,newk) = 1;
                    val = val + p(i,j,newk); 
                    % extract the assigned resources
                    beta(k) = beta(k)-b(i,j);
                    sigma(k) = sigma(k)-s(i,j);
                    kapa(k) = kapa(k)-c(i,j);
                else
                    x(i,:,:) = 0;
                    val = valpre;
                    return
                end
            end
        else 
            x(i,:,:) = 0;
            val = valpre;
            return
        end
    end
end

end

function closeNodes = findCloseNodes(distanceMatrix, targetNode, maxDelay)
    % Get the distances from the target node to all other nodes
    distancesToTarget = distanceMatrix(targetNode, :);  
    % Find nodes that are closer than the given maximum distance
    closeNodes = find(distancesToTarget < maxDelay);
end

function suitable = suitableNodes(nodes,b,s,c,beta,sigma,kapa)
        suitable = find(beta(nodes) >= b);
        suitable = find(sigma(suitable) >= s);
        suitable = find(kapa(suitable) >= c);
end
