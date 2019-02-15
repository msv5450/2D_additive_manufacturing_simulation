function [N] =...
    computeQuad2dCFDshapeFunctions(ipCoords, nodesCoordList, elementLength)

% transform ipCoords from physical (x,y) coordinates to parent (xi,eta)
% coordinates. Lineraly scale because elelemnts are structured.
dx = ipCoords(1) - nodesCoordList(1,1);
dy = ipCoords(2) - nodesCoordList(1,2);

xi_ip = -1.0 + 2.0*dx/elementLength;
eta_ip = -1.0 + 2.0*dy/elementLength;

N = zeros(4,1);
% shape functions at the given point inside element
N(1) = 0.25*(1-xi_ip)*(1-eta_ip);
N(2) = 0.25*(1+xi_ip)*(1-eta_ip);
N(3) = 0.25*(1+xi_ip)*(1+eta_ip);
N(4) = 0.25*(1-xi_ip)*(1+eta_ip);

