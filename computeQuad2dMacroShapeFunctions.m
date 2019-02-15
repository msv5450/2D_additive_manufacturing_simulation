function [xq, wq, N, dNdx] = computeQuad2dMacroShapeFunctions(coords)

% This function takes the element nodal coordinates (4x2) as input,
% and returns the integration point coordinates xq (16x2) and
% weights wq (16x1), shape functions for each node at each
% quadrature point, and shape function derivatives. On return,
% N(iq,inode) gives the shape function for node inode at quadrature
% point iq, and dNdx(iq,inode,idim) is the derivative in the idim
% direction of the shape function for node inode at quadrature
% point iq.

coord_1 = sqrt(3.0/7.0 - 2.0/7.0*sqrt(6.0/5.0));
coord_2 = sqrt(3.0/7.0 + 2.0/7.0*sqrt(6.0/5.0));
coords16 = [-coord_2, -coord_1, coord_1, coord_2];

w1 = (18.0 + sqrt(30.0)) / 36.0;
w2 = (18.0 - sqrt(30.0)) / 36.0;
w = [w2, w1, w1, w2];

aa = zeros(16,2);
indx = 1;

for i=1:4
    for j=1:4
        aa(indx,1) = coords16(j);
        aa(indx,2) = coords16(i);
        indx = indx + 1;
    end
end

J = zeros(2,2);
N = zeros(16,4);
xq = zeros(16,2);
wq = zeros(16,1);
dNdx = zeros(16,4,2);

for i=1:16
    a = aa(i,1);
    b = aa(i,2);
    ii = rem((i-1),4) + 1;
    jj = floor((i-1)/4) + 1;
    
    % x coordinate of 9 quadrature points xq
    xq(i,1) = coords(1,1)*0.25*(1-a)*(1-b) + coords(2,1)*0.25*(1+a)*(1-b) + coords(3,1)*0.25*(1+a)*(1+b)+...
              coords(4,1)*0.25*(1-a)*(1+b);
    
    % y coordinate of 9 quadrature points xq
    xq(i,2) = coords(1,2)*0.25*(1-a)*(1-b) + coords(2,2)*0.25*(1+a)*(1-b) + coords(3,2)*0.25*(1+a)*(1+b)+...
              coords(4,2)*0.25*(1-a)*(1+b);
    
    % Defining Jacobian Matrix
    J(1,1) = (coords(3,1)+coords(2,1)-coords(4,1)-coords(1,1)) + b*(coords(1,1)-coords(2,1)+coords(3,1)-coords(4,1));
    J(2,2) = (coords(3,2)+coords(4,2)-coords(2,2)-coords(1,2)) + a*(coords(1,2)-coords(2,2)+coords(3,2)-coords(4,2));
    J(2,1) = (coords(3,1)+coords(4,1)-coords(2,1)-coords(1,1)) + a*(coords(1,1)-coords(2,1)+coords(3,1)-coords(4,1));
    J(1,2) = (coords(3,2)+coords(2,2)-coords(4,2)-coords(1,2)) + b*(coords(1,2)-coords(2,2)+coords(3,2)-coords(4,2));
    J = J/4;
    
    % Weight function equals det(J)
    wq(i) = det(J) * w(ii) * w(jj);
    
end

% constructing shape functions on each quadrature point
for i=1:4      %loop over nodes
    for j=1:16  %loop over integration points
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
for i=1:4       %loop over nodes
    for j=1:16  %loop over integration points
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


      
