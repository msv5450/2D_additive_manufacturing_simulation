function [LHSbc,RHSbc] = applyLaserFlux(xnodes, surfconn, LHS, RHS)

LHSbc = LHS;
RHSbc = RHS;

ne = size(surfconn,1);   % number of surface elements
nen = size(surfconn,2);  % number of nodes per surface element
nq  = 2;                 % number of element integration points
beamEff = 0.75;
f = 2.0;
rb = 0.05;
Q = 2000.0;

% Loop over elements
for ielt = 1:ne
    %[1x2] X coordinates of two nodes
    coords = xnodes(surfconn(ielt,:)',:);  
    nodeID = surfconn(ielt,:);
    [xq, wq, N] = compute1dFemShapeFunctions(coords);
    
    for iq = 1:nq 
        x_ip = xq(iq);
        % Laser is fixed at 1.1 
        r2 = (1.1 - x_ip) * (1.1 - x_ip);
        qLaser = -f * beamEff * Q / (rb*sqrt(2.0*pi)) * exp(-f*r2 / (2.0*rb*rb));
        
        for I = 1:nen
            RHSbc(nodeID(I)) = RHSbc(nodeID(I)) + N(iq,I) * qLaser * wq(iq);
        end
    end
    
end