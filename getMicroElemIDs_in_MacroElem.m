function slaveMicroElemList = getMicroElemIDs_in_MacroElem(nxMacro, nyMacro, nxMicro, nyMicro)

% given the Macro element ID you can get the element IDs of the micro
% elements that are inside a Macro element
 
rx = nxMicro / nxMacro;
ry = nyMicro / nyMacro;
Macro_ne = nxMacro * nyMacro;    % total number of Macro elements

% each Macro element has rx*ry micro elements
slaveMicroElemList = zeros(Macro_ne, rx*ry);

% Loop over Macro elements
for ielt = 1:Macro_ne
    y = floor((ielt-1) / nxMacro);
    x = floor(mod(ielt-1, nxMacro));
    starting_microID = y*nxMacro*rx*ry + x*rx;
    
    c = 1;
    for i = 1:ry
        k = starting_microID + (i-1)*nxMacro*rx;
        for j = 1:rx
            slaveMicroElemList(ielt,c) = k+j;
            c = c+1;
        end
    end
end
