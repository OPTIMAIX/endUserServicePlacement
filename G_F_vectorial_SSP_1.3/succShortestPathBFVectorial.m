function x_fs = succShortestPathBFVectorial(G,r, b, c, d)
% costs are all positive

% G is the conectivity graph
% r(i) is the capacity of the graph in dimension i
% flow is composed of multi-dimensional quantums and is described by the matrix [q_i,d1_i,...,dn_i]; the flow has q units of quantum i, which is composed of dj_i units in dimension j
% b is a matrix where $b_ij$ is the number of services of type $i$
% "produced" by node $j$
% c is the matrix of costs $c_ij$ of each edge $ij$ in G
% d is a matrix of service types, $d_ij$ is the amount of resources $j$ used by
% one service of type $i$

n = length(b);

x_fs = zeros(n,n,size(d,1));  %flow is in terms of services,
x_fr = zeros(n,n,size(r,3));  %flow is in terms of resources,

%
G_x = G;
r   = r;
b   = b;
d   = d;
G_x_f = r;
G_x_c = c;

pie = zeros(1,n);

e = b;
% E = 1;
% D = n;

[~, E] = find(e>0);
[~, D] = find(e<0);

while ~isempty(E)
    
    k = E(end);
    l = D(end); 
    skipPath = false;
    
    %determine shortest path distances d(j) from node s to all other nodes in G(x) with respect to the reduced costs $c_{ij}^\pi$;
    [p,g,~] = BFMSpathOT(G_x_c,k,G_x); % FIXME dijkstra would work faster

    if g(l) > intmax('uint32') - 10 
        break; 
    end
    
    pie = pie - g'; 
    
    % one min per dimension
    minP = Inf(1,size(r,3));
    
    v = l; %it's L the letter, not one the number :-)
    
    while p(v) ~= 0
        u = p(v);
        minP = min(minP,squeeze(G_x_f(u,v,:))');
        if find(minP < d(find(e(:,k)),:))
            G_x(u,v) = 0;
            skipPath = true;
        end
        v = u;
    end
    
    if skipPath; continue; end
    
    gamma = min(min(sum(e(:,k).*d), -sum(e(:,l).*d)), minP);
    
    % computes the optimal combination of services for this minimums 
    [x,fval,exitflag,output] = ksp(e(:,k),gamma,ones(size(e,1),1),d'); 
    
    if isempty(find(x))
        break;
    end
    
    v = l;
    while p(v) ~= 0
        u = p(v);
        for i = 1:size(x_fs,3)
            if(G(u,v)) % it's a forwards edge
                x_fs(u,v,i) = x_fs(u,v,i) + x(i); 
                x_fr(u,v,:) = x_fr(u,v,:) + x(i).*d(i);
            else       %it's a backwards edge 
                x_fs(v,u,i) = x_fs(v,u,i) - x(i);   %FIXME lo que tengo que sacar en los links backwards es servicios que equivalgan a la capacidad de los que hag pasar por ahi. 
                x_fr(u,v,:) = x_fr(u,v,:) - x(i).*d(i);
            end
        end
        v = u;
    end

    % residual graph 
    G_x_f = (r - x_fr);
    G_x_f = G_x_f + permute(x_fr, [2 1 3]);
    
    G_x = zeros(n,n);  
    for i = 1: size(G_x_f,3) 
        G_x = or(G_x,G_x_f(:,:,i));
    end
        
    [i,j,v] = find(G_x); %capaz que s epuede poner G_x_c aca
    
    for k = 1:length(i)
        if G_x(i(k),j(k)) 
            if G(i(k),j(k))
                G_x_c(i(k),j(k)) = G_x_c(i(k),j(k)) - pie(i(k)) + pie(j(k));
            else
                G_x_c(i(k),j(k)) = 0;
            end
        end
    end
    
    for i = 1:size(x_fs,3)
        e(i,:) = b(i,:) + (sum(x_fs(:,:,i),1) - sum(x_fs(:,:,i),2)'); 
    end 
    
    [row, E] = find(e>0);
    [row, D] = find(e<0);
   
end

end
