function h = plotScalar(xnodes,phi)

x = xnodes(:,1);
y = xnodes(:,2);
xu = unique(x);
yu = unique(y);
nx = length(xu);
ny = length(yu);
[x,y] = meshgrid(xu,yu);
phiPlot = reshape(phi,nx,ny)';
h = contourf(x,y,phiPlot,30,'edgecolor','none');
colorbar;
colormap jet;
