function newMicroPhi = rearrangeMicroPhi(phi_micro, nEnrichment, nnEnrichment, slaveMicroElemList,...
                       micro_nconn, enrich_nconn, nMacro, nnMicro)

newMicroPhi = zeros(nnMicro,1);
micro_nen = size(micro_nconn,2);  % number of micro nodes per micro element

% Loop over Macro elements
for MacroEID = 1:nMacro
    
    % Get element ID of micro elements inside this Macro element
    microElemIDs = slaveMicroElemList(MacroEID,:);
    offSet = (MacroEID-1) * nnEnrichment;
    
    % Loop over micro elements inside this Macro element
    for iMicroElem = 1:nEnrichment
        
        % Get global micro node IDs of this micro element
        micro_nodes = micro_nconn(microElemIDs(iMicroElem),:);
        % Get local micro node IDs; ranging from 1 to nnEnrich
        enrich_nodes = enrich_nconn(iMicroElem,:);
        
        % Loop over nodes of micro element
        for iMicroNode = 1:micro_nen
            
            indx_local = enrich_nodes(iMicroNode);
            indx_global = micro_nodes(iMicroNode);
            phi = phi_micro(offSet + indx_local);
            newMicroPhi(indx_global) = phi;
        end
    end
end