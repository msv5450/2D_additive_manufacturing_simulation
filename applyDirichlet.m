function [LHSbc,RHSbc] = applyDirichlet(xnodes, surfconn, f, LHS, RHS, value)

LHSbc = LHS;
RHSbc = RHS;

% Get unique nodes on surface
nodeList = unique(surfconn);

% Get x and y values
x = xnodes(nodeList,1);
y = xnodes(nodeList,2);

% Compute dirichlet values from function f
ubc = f(x,y,value);

% Modify LHS rows
LHSbc(nodeList,:) = 0;
LHSbc(nodeList,nodeList) = speye(length(nodeList));

% Modify RHS
RHSbc(nodeList) = ubc;
