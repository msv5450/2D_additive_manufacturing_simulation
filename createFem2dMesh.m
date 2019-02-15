function [xnodes, nconn, surfconnD, surfconnR, surfconnU, surfconnL] =...
    createFem2dMesh(Lx, Ly, nx, ny)

nxNodes = nx + 1;
nyNodes = ny + 1;

dx = Lx/nx;
dy = Ly/ny;

nn = nxNodes * nyNodes;
ne = nx * ny;

xnodes = zeros(nn, 2);  % node coordinates
nconn  = zeros(ne, 4);  % Element connectivity
surfconnD = zeros(nx, 2);
surfconnR = zeros(ny, 2);
surfconnU = zeros(nx, 2);
surfconnL = zeros(ny, 2);

% origin [0.0 0.0] is at the center of the domain
% xnodes gives xy coordinates given the nodeID
for j = 1:nyNodes
    for i = 1:nxNodes
        xnodes((j-1)*nxNodes + i, :) = [((i-1)*dx - Lx*0.5) ((j-1)*dy - Ly*0.5)];
    end
end

% nconn gives 4 nodes that belong to an element, given the eID  
for je = 1:ny
    for ie = 1:nx
        nconn((je-1)*nx + ie, 1) = (je-1)*nxNodes + ie;
        nconn((je-1)*nx + ie, 2) = (je-1)*nxNodes + ie + 1;
        nconn((je-1)*nx + ie, 3) = (je  )*nxNodes + ie + 1;
        nconn((je-1)*nx + ie, 4) = (je  )*nxNodes + ie;
    end
end

% Surface D (down)
for ie = 1:nx
    surfconnD(ie, :) = [ie ie+1];
end

% Surface R (right)
for ie = 1:ny
    surfconnR(ie, :) = [ie*nxNodes (ie+1)*nxNodes];
end

% Surface U (up)
for ie = 1:nx
    surfconnU(ie, :) = [nn-ie+1 nn-ie];
end

% Surface L (left)
for ie = 1:ny
    surfconnL(ie, :) = [(nyNodes-ie)*nxNodes+1 ...
                        (nyNodes-ie-1)*nxNodes+1];
end


    