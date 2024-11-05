function [x,fval,exitflag,output] = ksp(d,c,p,w)
%KSP Summary of this function goes here
%   Detailed explanation goes here

%%
% maximize $\sum_{j=1}^{n} p_{j} x_{j}$
%%
% subject to $\sum_{j=1}^{n} w_{i j} x_{j} \leq c_{i}, \quad i=1, \ldots,
% d$
%%
% $0 \leq x_{j} \leq d_{j}$
%%
% $\quad x_{j} \in\{0,1\}, \quad j=1, \ldots, n$
%%

% d is the vector with a number of services per type 
% c is the vector with the max capacity at each dimension
% p is a vector with the profit of each service type 
% w_ij is a matrix with the weight in dimension i of service j 

% linprog solves 
%% 
% $\min_{x} f^{T} x$  such that
%%
% $A \cdot x \leq b$
%%
% $Aeq \cdot x = beq$
%%
% $lb \leq x \leq ub$
% 
% 

f = -p;
intcon = 1:length(f);
A = w;
b = double(c);
lb = zeros(size(f));
ub = d;

options = optimoptions('intlinprog','Display','off');

[x,fval,exitflag,output] = intlinprog(f,intcon,A,b,[],[],lb,ub, [],options);
        
end

