function [xnodes, nconn] = createCFD2dMesh(Lx, Ly, nx, ny)

    nxNodes = nx + 1;
    nyNodes = ny + 1;

    dx = Lx/nx;
    dy = Ly/ny;

    nn = nxNodes * nyNodes;
    ne = nx * ny;

    xnodes = zeros(nn, 2);              % node coordinates
    nconn  = zeros(ne, 4);              % Element connectivity
    centroids = zeros(ne, 3);           % eID + Element centroid (x,y)

    % origin [0.0 0.0] is at the center of the domain
    % xnodes gives xy coordinates of NALU mesh given the nodeID
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
    
end
