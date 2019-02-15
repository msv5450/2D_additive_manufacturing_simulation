function BlockLHSbc = applyZeroDirichletBC_LHS(surfconn, BlockLHS)

BlockLHSbc = BlockLHS;

% Get unique nodes on surface
nodeList = unique(surfconn);

% Modify LHS rows only
% RHS of micro equation will be modified in the iteration loop
BlockLHSbc(nodeList,:) = 0;
BlockLHSbc(nodeList,nodeList) = speye(length(nodeList));
