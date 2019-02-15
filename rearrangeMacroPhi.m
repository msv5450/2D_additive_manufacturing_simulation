function [newMacroPhi,newLf] = rearrangeMacroPhi(phi_Macro, Lf, Macro_nconn)

Macro_ne = size(Macro_nconn,1);   % number of Macro elements
Macro_nen = size(Macro_nconn,2);  % number of Macro nodes per Macro element

newMacroPhi = zeros(Macro_nen*Macro_ne,1);
newLf = zeros(Macro_nen*Macro_ne,1);

% Loop over Macro elements
for MacroEID = 1:Macro_ne
    
    MacroNodeIDs = Macro_nconn(MacroEID,:);
    offSet = (MacroEID - 1)*Macro_nen;
    
    % Loop over Macro elements
    for i = 1:Macro_nen
        
        % get value of phi at each element node and assemble to global
        phi = phi_Macro(MacroNodeIDs(i));
        liquidFraction = Lf(MacroNodeIDs(i));
        newMacroPhi(offSet + i) = phi;
        newLf(offSet + i) = liquidFraction;
    end
end