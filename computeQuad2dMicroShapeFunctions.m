function [xq, wq, N, dNdx] = computeQuad2dMicroShapeFunctions(coords)

% This function takes the element nodal coordinates (4x2) as input,
% and returns the integration point coordinates xq (4x1) and
% weights wq (4x1), shape functions for each node at each
% quadrature point, and shape function derivatives. On return,
% N(iq,inode) gives the shape function for node inode at quadrature
% point iq, and dNdx(iq,inode,idim) is the derivative in the idim
% direction of the shape function for node inode at quadrature
% point iq.

aa(1,1) = -sqrt(3)^-1;
aa(1,2) = -sqrt(3)^-1;
aa(2,1) = sqrt(3)^-1;
aa(2,2) = -sqrt(3)^-1;
aa(3,1) = sqrt(3)^-1;
aa(3,2) = sqrt(3)^-1;
aa(4,1) = -sqrt(3)^-1;
aa(4,2) = sqrt(3)^-1;

J = zeros(2,2);
N = zeros(4,4);
xq = zeros(4,2);
wq = zeros(4,1);
dNdx = zeros(4,4,2);

for i=1:4
    a = aa(i,1);
    b = aa(i,2);
    
    % x coordinate of 4 quadrature points xq
    xq(i,1) = coords(1,1)*0.25*(1-a)*(1-b) + coords(2,1)*0.25*(1+a)*(1-b) + coords(3,1)*0.25*(1+a)*(1+b)+...
              coords(4,1)*0.25*(1-a)*(1+b);
    
    % y coordinate of 4 quadrature points xq
    xq(i,2) = coords(1,2)*0.25*(1-a)*(1-b) + coords(2,2)*0.25*(1+a)*(1-b) + coords(3,2)*0.25*(1+a)*(1+b)+...
              coords(4,2)*0.25*(1-a)*(1+b);
    
    % Defining Jacobian Matrix
    J(1,1) = (coords(3,1)+coords(2,1)-coords(4,1)-coords(1,1)) + b*(coords(1,1)-coords(2,1)+coords(3,1)-coords(4,1));
    J(2,2) = (coords(3,2)+coords(4,2)-coords(2,2)-coords(1,2)) + a*(coords(1,2)-coords(2,2)+coords(3,2)-coords(4,2));
    J(2,1) = (coords(3,1)+coords(4,1)-coords(2,1)-coords(1,1)) + a*(coords(1,1)-coords(2,1)+coords(3,1)-coords(4,1));
    J(1,2) = (coords(3,2)+coords(2,2)-coords(4,2)-coords(1,2)) + b*(coords(1,2)-coords(2,2)+coords(3,2)-coords(4,2));
    J = J/4;
    
    % Weight function equals det(J)
    wq(i) = det(J);
    
end

% constructing shape functions on each quadrature point
for i=1:4
    for j=1:4
        a = aa(j,1);
        b = aa(j,2);
        
        if (i==1)
            N(j,i) = 0.25*(1-a)*(1-b);
        elseif(i==2)
            N(j,i) = 0.25*(1+a)*(1-b); 
        elseif(i==3)
            N(j,i) = 0.25*(1+a)*(1+b);
        else
            N(j,i) = 0.25*(1-a)*(1+b);
        end
    end
end

d = zeros(2,1);
p = zeros(2,1);
% constructing shape functions derivatives on each quadrature point
for i=1:4
    for j=1:4
        a = aa(j,1);
        b = aa(j,2);
        
        J(1,1) = (coords(3,1)+coords(2,1)-coords(4,1)-coords(1,1)) + b*(coords(1,1)-coords(2,1)+coords(3,1)-coords(4,1));
        J(2,2) = (coords(3,2)+coords(4,2)-coords(2,2)-coords(1,2)) + a*(coords(1,2)-coords(2,2)+coords(3,2)-coords(4,2));
        J(2,1) = (coords(3,1)+coords(4,1)-coords(2,1)-coords(1,1)) + a*(coords(1,1)-coords(2,1)+coords(3,1)-coords(4,1));
        J(1,2) = (coords(3,2)+coords(2,2)-coords(4,2)-coords(1,2)) + b*(coords(1,2)-coords(2,2)+coords(3,2)-coords(4,2));
        J = J/4;
        
        if (i==1)
            p(1) = 0.25*(-1+b);
            p(2) = 0.25*(-1+a);
        elseif (i==2)
            p(1) = 0.25*(1-b);
            p(2) = 0.25*(-1-a);
        elseif (i==3)
            p(1) = 0.25*(1+b);
            p(2) = 0.25*(1+a);
        else
            p(1) = 0.25*(-1-b);
            p(2) = 0.25*(1-a);
        end

        d = inv(J)*p;
        dNdx(j,i,1) = d(1);
        dNdx(j,i,2) = d(2);
    end
end


