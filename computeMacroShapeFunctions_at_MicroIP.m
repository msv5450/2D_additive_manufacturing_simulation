function [N, dNdx] = computeMacroShapeFunctions_at_MicroIP(hxMacro, hyMacro, nEnrichment,...
                     Macro_nconn, micro_nconn, Macro_xnodes, micro_xnodes, slaveMicroElemList)

% This function takes number of Macro elements as input and returns Macro 
% shape functions and their derivatives at the integration
% points of the micro elements that make up each Macro element.
% N [MacroEID, iq, inode, Enrichment ID]
% dN/dx [MacroEID, iq, inode, Enrichment ID, direction]

aa(1,1) = -sqrt(3)^-1;
aa(1,2) = -sqrt(3)^-1;
aa(2,1) = sqrt(3)^-1;
aa(2,2) = -sqrt(3)^-1;
aa(3,1) = sqrt(3)^-1;
aa(3,2) = sqrt(3)^-1;
aa(4,1) = -sqrt(3)^-1;
aa(4,2) = sqrt(3)^-1;

nq = 4;
NN = 0.0;
ne = size(Macro_nconn,1);   % number of Macro elements
p = zeros(2,1);
d = zeros(2,1);
J = zeros(2,2);
N = zeros(ne,4,4,nEnrichment);
dNdx = zeros(ne,4,4,nEnrichment,2);

% loop over Macro elements
for MacroEID=1:ne
    % Coordinates for Macro element nodes (size: 4X2)
    coordsMacro = Macro_xnodes(Macro_nconn(MacroEID,:)',:);
    % coordinates of the first node
    x0 = coordsMacro(1,1);
    y0 = coordsMacro(1,2);
    % get eId of micro elemets inside this Macro element (1 X nEnrichment)
    microEIDs = slaveMicroElemList(MacroEID,:);
    
    % loop over micro elements
    for iEnrich=1:nEnrichment
        microEID = microEIDs(iEnrich);
        % Coordinates for micro element nodes (size: 4X2)
        coords_micro = micro_xnodes(micro_nconn(microEID,:)',:);
        
        % loop over micro integration points
        for iq = 1:nq 
            a = aa(iq,1);
            b = aa(iq,2);
    
            % x coordinate of the micro quadrature point
            xq = coords_micro(1,1)*0.25*(1-a)*(1-b) + coords_micro(2,1)*0.25*(1+a)*(1-b) +...
                 coords_micro(3,1)*0.25*(1+a)*(1+b) + coords_micro(4,1)*0.25*(1-a)*(1+b);
    
            % y coordinate of the micro quadrature point
            yq = coords_micro(1,2)*0.25*(1-a)*(1-b) + coords_micro(2,2)*0.25*(1+a)*(1-b) +...
                 coords_micro(3,2)*0.25*(1+a)*(1+b) + coords_micro(4,2)*0.25*(1-a)*(1+b);
             
            % get xq and yq in (xi, eta) parent domain of Macro element
            % (a,b) are coordinates of ip with respect to the Micro parent
            % domain
            dx = xq - x0;
            dy = yq - y0;
            xi_ip  = -1.0 + 2.0*dx/hxMacro;
            eta_ip = -1.0 + 2.0*dy/hyMacro;
            
            % compute Jacobian of Macro element at the micro ip
            J(1,1) = (coordsMacro(3,1)+coordsMacro(2,1)-coordsMacro(4,1)-coordsMacro(1,1)) +...
                     eta_ip*(coordsMacro(1,1)-coordsMacro(2,1)+coordsMacro(3,1)-coordsMacro(4,1));
                 
            J(2,2) = (coordsMacro(3,2)+coordsMacro(4,2)-coordsMacro(2,2)-coordsMacro(1,2)) +...
                     xi_ip*(coordsMacro(1,2)-coordsMacro(2,2)+coordsMacro(3,2)-coordsMacro(4,2));
                 
            J(2,1) = (coordsMacro(3,1)+coordsMacro(4,1)-coordsMacro(2,1)-coordsMacro(1,1)) +...
                     xi_ip*(coordsMacro(1,1)-coordsMacro(2,1)+coordsMacro(3,1)-coordsMacro(4,1));
                 
            J(1,2) = (coordsMacro(3,2)+coordsMacro(2,2)-coordsMacro(4,2)-coordsMacro(1,2)) +...
                     eta_ip*(coordsMacro(1,2)-coordsMacro(2,2)+coordsMacro(3,2)-coordsMacro(4,2));
            J = J/4;
            
            % loop over Macro nodes
            for MacroNode=1:4
                
                if (MacroNode == 1)
                    NN = 0.25*(1-xi_ip)*(1-eta_ip);
                    p(1) = 0.25*(-1+eta_ip);
                    p(2) = 0.25*(-1+xi_ip);
                elseif(MacroNode == 2)
                    NN = 0.25*(1+xi_ip)*(1-eta_ip);
                    p(1) = 0.25*(1-eta_ip);
                    p(2) = 0.25*(-1-xi_ip);   
                elseif(MacroNode == 3)
                    NN = 0.25*(1+xi_ip)*(1+eta_ip);
                    p(1) = 0.25*(1+eta_ip);
                    p(2) = 0.25*(1+xi_ip);
                else
                   NN = 0.25*(1-xi_ip)*(1+eta_ip);
                   p(1) = 0.25*(-1-eta_ip);
                   p(2) = 0.25*(1-xi_ip);
                end
                
                d = inv(J)*p;
                N(MacroEID, iq, MacroNode, iEnrich) = NN;
                dNdx(MacroEID, iq, MacroNode, iEnrich, 1) = d(1);
                dNdx(MacroEID, iq, MacroNode, iEnrich, 2) = d(2);
            end
        end    
    end
end









