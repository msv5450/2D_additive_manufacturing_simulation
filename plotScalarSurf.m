function h = plotScalarSurf(xnodes,phi)

x = xnodes(:,1);
y = xnodes(:,2);
xu = unique(x);
yu = unique(y);
nx = length(xu);
ny = length(yu);
[x,y] = meshgrid(xu,yu);
phiPlot = reshape(phi,nx,ny)';
h = surf(x,y,phiPlot);
colorbar;
