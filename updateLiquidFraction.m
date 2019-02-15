function Lf = updateLiquidFraction(phi, solidus, liquidus, nnMacro)

Lf = zeros (nnMacro,1);
% loop over Macro nodes
for inode = 1:nnMacro
    if (phi(inode) >= liquidus)
        Lf(inode) = 1.0;
    elseif (phi(inode) > solidus && phi(inode) < liquidus)
        Lf(inode) = (phi(inode) - solidus) / (liquidus - solidus);
    else
        Lf(inode) = 0.0;
    end
end