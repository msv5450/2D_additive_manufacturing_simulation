function [xq, wq, N] = compute1dFemShapeFunctions(coords)

% using 2 quadrature (not 3) points for applying laser flux BC
nn = 2;
nip = 2;
xin = [-1 1];

a = 1./sqrt(3);
xiq = [-a a];

N = zeros(nip,nn);
N(:,1) = 1/2 * (1 - xiq(:));
N(:,2) = 1/2 * (1 + xiq(:));

x1 = coords(1);
x2 = coords(2);
h = x2 - x1;

wq = [h/2 h/2];

xq = [(N(1,1)*x1 + N(1,2)*x2) (N(2,1)*x1 + N(2,2)*x2)];