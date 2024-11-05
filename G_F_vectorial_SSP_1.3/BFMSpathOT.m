function [p,D,iter] = BFMSpathOT(G,r,C)

[tail,head,W] = find(C);           % Get arc list {u,v,duv, 1:m} from conectivity graph.
[~,n] = size(C); 
m = nnz(C);
p(1:n,1) = 0;                      % Shortest path tree of parent pointers
%p = uint32(zeros(n,1));
D(1:n,1) = Inf;                    % Sh. path distances from node i=1:n to the root of the sp tree
%D = uint32(Inf(n,1));
p(r)=0; D(r)=0;                    % Set the root of the SP tree (p,D)
for iter = 1:n-1                   % Converges in <= n-1 iters if no 
    optimal = true;                % negative cycles exist in G
    for arc = 1:m                  % O(m) optimality test. 
        u = tail(arc); 
        v = head(arc);
        duv = G(u,v);
        if D(v) > D(u) + duv;      %
           D(v) = D(u) + duv;      % Sp Tree not optimal: Update (p,D) 
           p(v) = u;               %
           optimal = false;
        end
    end 
    if optimal
        return                     % SP Tree p is optimal;
    end
end 

