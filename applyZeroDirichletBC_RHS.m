function BlockRHSbc = applyZeroDirichletBC_RHS(surfconn, BlockRHS, nnEnrichment, nMacro)

BlockRHSbc = BlockRHS;

% Get unique nodes on surface
nodeList = unique(surfconn);

% loop over Macro elements (blocks)
for MacroEID = 1:nMacro
    
    % find the offSet index
    offSet = (MacroEID - 1) * nnEnrichment;
    nodeListOffSet = nodeList + offSet;
    
    % Modify RHS
    BlockRHSbc(nodeListOffSet) = 0.0;
end
    