function [x_opt, f_opt, exitflag, output] = milpsolver(bw,s,c,beta,sigma,kapa,d,D,p)


% Decision variables
x = binvar(n, max(cellfun(@numel, F)), m, 'full'); % If microservice j of service i is assigned to infrastructure element k
y = binvar(n, 1, 'full'); % If all microservices of service i are assigned to any infrastructure element
z = binvar(n, max(cellfun(@numel, F)), n-1, m, m, 'full'); % If the pair of microservices (j1,j2) of service i is assigned to the pair of infrastructure elements (k1,k2)

% Objective function
f = reshape(p, [n * max(cellfun(@numel, F)) * m, 1]);
intcon = reshape(x, [n * max(cellfun(@numel, F)) * m, 1]);
A = [];
b = [];
Aeq = [];
beq = [];
lb = zeros(size(intcon));
ub = ones(size(intcon));

% Storage, Bandwidth, and Computation Constraints
for k = R
    A = [A; reshape(s, [n * max(cellfun(@numel, F)), 1])' * reshape(x(:, :, k), [n * max(cellfun(@numel, F)), 1]) <= sigma(k)];
    A = [A; reshape(bw, [n * max(cellfun(@numel, F)), 1])' * reshape(x(:, :, k), [n * max(cellfun(@numel, F)), 1]) <= beta(k)];
    A = [A; reshape(c, [n * max(cellfun(@numel, F)), 1])' * reshape(x(:, :, k), [n * max(cellfun(@numel, F)), 1]) <= kapa(k)];
end

% Exclusion Constraints
for i = S
    for j = F{i}
        A = [A; sum(x(i, j, :)) <= 1];
    end
end

% Inclusion Constraints
for i = S
    for j = F{i}
        A = [A; x(i, j, :) <= y(i)];
    end
end

for i = S
    A = [A; sum(reshape(x(i, :, :), [1, n * m])) >= numel(F{i}) * y(i)];
end

% Delay Constraints
for i = S
    for pIdx = 1:length(P{i})
        (j1, j2) = P{i}(pIdx);
        for qIdx = 1:length(Q)
            (k1, k2) = Q{qIdx};
            A = [A; D(k1, k2) * z(i, j1, j2, k1, k2) <= d(i, pIdx)];
            A = [A; D(k1, k2) * z(i, j1, j2, k2, k1) <= d(i, pIdx)];
            A = [A; z(i, j1, j2, k1, k2) <= x(i, j1, k1)];
            A = [A; z(i, j1, j2, k1, k2) <= x(i, j2, k2)];
            A = [A; z(i, j1, j2, k1, k2) >= x(i, j1, k1) + x(i, j2, k2) - 1];
            A = [A; z(i, j1, j2, k2, k1) <= x(i, j1, k2)];
            A = [A; z(i, j1, j2, k2, k1) <= x(i, j2, k1)];
            A = [A; z(i, j1, j2, k2, k1) >= x(i, j1, k2) + x(i, j2, k1) - 1];
        end
    end
end

% Solve the MILP problem
options = optimoptions('intlinprog', 'Display', 'iter');
[x_opt, f_opt, exitflag, output] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub, options);

% Extract the results
% You can access the assignment matrix x_opt to get the optimal assignment of microservices to infrastructure elements.


end

