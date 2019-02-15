clear all;
close all;

% This code solves a homogeneous advcetion diffusion for a melt pool problem
% in 2D using variational multi-scale FEM. The strong form of the equation solved is:

%   d/dx_i * (Kappa T/dx_i) = rho*Cp*u.Div(T) on 0<=x<=L, 0<=y<=L
 
% The velocity field is solved by NALU for the same mesh (240x40)
% and is mapped to quadrature points of the coarser Matlab (60x10) Macro
% grid. Each Matlab Macro element has 16 integration points and several 
% micro elements. Each micro element has 4 integration points. The domain
% of integration points will be limited to the NALU melt pool region in
% another code. For now, the micro and Macro domains are the same.

% Latent heat is advected instead of using the equivalent heat capacitance
% model. Both micro and Macro scale are solved iteratively.

%% Get the velocity field from NALU output
velocity_fileName = '/home/MJ/MATLAB_codes/variational_multiscale/VMS_meltpool/velocity_field/meltpool_80.csv';
NALU_data = csvread(velocity_fileName, 1);
NALU_data = round(NALU_data, 6);
sort_NALU_data = sortrows(NALU_data,[7 5]);   
vnodes = zeros(size(NALU_data,1), 2);
tempNodes = zeros(size(NALU_data,1), 1);
tempNodes(:,1) = sort_NALU_data(:,1);
vnodes(:,1) = sort_NALU_data(:,2);   % x velocity component
vnodes(:,2) = sort_NALU_data(:,4);   % y velocity component

%% Problem Parameters
Lx = 3.0;
Ly = 0.5;
kappa  = 0.3;           % thermal conductivity
cp = 0.7;               % specific heat
rho = 1.0;              % density
latent = 297.0;         % latent heat of fusion
solidus = 1563.0;       % solidus temperature
liquidus = 1623.0;      % liquidus temperature
v_beam = -2.5;          % scan speed

%% Mesh Parameters
nxMacro = 48;                 % number of Macro elements in x direction
nyMacro = 8;                  % number of Macro elements in y direction
nxMicro = 240;                 % number of micro elements in x direction
nyMicro = 40;                  % number of micro elements in y direction
hxMacro = Lx / nxMacro;        % Macro element length in x direction
hyMacro = Ly / nyMacro;        % Macro element length in y direction
rx = nxMicro / nxMacro;
ry = nyMicro / nyMacro;
nEnrichment = rx * ry;          % number of micro elements in each Macro element
nnEnrichment = (rx+1) * (ry+1); % number of micro nodes in each Macro element

nnxMacro = nxMacro + 1;         % number of Macro nodes in x direction
nnyMacro = nyMacro + 1;         % number of Macro nodes in y direction
nnMacro  = nnxMacro * nnyMacro; % total number of Macro nodes (unknowns)
nnxMicro = nxMicro + 1;         % number of micro nodes in x direction
nnyMicro = nyMicro + 1;         % number of micro nodes in y direction
nnMicro  = nnxMicro * nnyMicro; % total number of micro nodes (unknowns)
nMacro = nxMacro * nyMacro;     % total number of Macro elements

% Assume that the micro mesh is not the same as fine NALU mesh
nx_CFD = 240;                   % number of elements in x direction; NALU mesh
ny_CFD = 40;                    % number of elements in y direction; NALU mesh
h = Lx / nxMacro;               % Macro element length
h_CFD = Lx / nx_CFD;            % Nalu element length

%% Create mesh
% Macro grid
[Macro_xnodes, Macro_nconn, Macro_surfconnD, Macro_surfconnR, Macro_surfconnU, Macro_surfconnL] =...
    createFem2dMesh(Lx, Ly, nxMacro, nyMacro);

% micro grid
[micro_xnodes, micro_nconn, micro_surfconnD, micro_surfconnR, micro_surfconnU, micro_surfconnL] =...
    createFem2dMesh(Lx, Ly, nxMicro, nyMicro);

% Enrichment grid: one Macro element as the domain with its micro elements
[enrich_xnodes, enrich_nconn, enrich_surfconnD, enrich_surfconnR, enrich_surfconnU, enrich_surfconnL] =...
    createFem2dMesh(hxMacro, hyMacro, rx, ry);

% Reconstruct NALU mesh
[xnodes_CFD, nconn_CFD] = createCFD2dMesh(Lx, Ly, nx_CFD, ny_CFD);

% List of all micro element IDs that are in a given Macro element ID
slaveMicroElemList = getMicroElemIDs_in_MacroElem(nxMacro, nyMacro, nxMicro, nyMicro);

% get shape functions and their derivatives of Macro elements at the
% integration of points of their micro elements
[N_Macro_at_microIP, dNdx_Macro_at_microIP] = computeMacroShapeFunctions_at_MicroIP(hxMacro, hyMacro,...
                                              nEnrichment, Macro_nconn, micro_nconn, Macro_xnodes,...
                                              micro_xnodes, slaveMicroElemList);
                                          
%% Initialize sparse matrix-vector problem
% Macro equation matrices
LHS_Macro = sparse(nnMacro,nnMacro);              % gets multiplied by phi_Macro
RHS_Macro = sparse(nnMacro,nnEnrichment*nMacro);  % gets multiplied by phi_micro
RHS_Macro_latent = sparse(nnMacro,nnMacro);       % gets multiplied by Lf
mainRHS_Macro = zeros(nnMacro,1);
Lf = zeros(nnMacro,1);                            % liquid fraction
phi_Macro = 295.0 * ones(nnMacro,1);              % Macro initial condition (unknowns)
phi_Macro_old = phi_Macro;
rearranged_phi_Macro = zeros(4*nMacro,1);         % rearranged temp for micro RHS
rearranged_Lf = zeros(4*nMacro,1);                % rearranged Lf for micro RHS


% micro equation matrices
LHS_micro = sparse(nnEnrichment*nMacro,nnEnrichment*nMacro); % gets multiplied by phi_micro
RHS_micro = sparse(nnEnrichment*nMacro,4*nMacro);            % gets multiplied by rearranged_phi_Macro
RHS_micro_latent = sparse(nnEnrichment*nMacro,4*nMacro);     % gets multiplied by rearranged_Lf
mainRHS_micro = zeros(nnEnrichment*nMacro,1);
phi_micro = zeros(nnEnrichment*nMacro,1);                    % micro initial condition (unknowns)

%% Populate global Diffusion LHS/RHS matrices of coefficients 
% Macro matrices
[LHS_Macro,RHS_Macro,RHS_Macro_latent] = computeFemMacroMatrices(Macro_nconn, micro_nconn, Macro_xnodes, micro_xnodes, nconn_CFD,...
                     enrich_nconn, xnodes_CFD, slaveMicroElemList, N_Macro_at_microIP, dNdx_Macro_at_microIP, nnEnrichment,...
                     nEnrichment, h_CFD, nx_CFD, ny_CFD, vnodes, v_beam, latent, cp, rho, kappa, LHS_Macro, RHS_Macro, RHS_Macro_latent);
                    
% micro matrices                   
[LHS_micro,RHS_micro,RHS_micro_latent] = computeFemMicroMatrices(Macro_nconn, micro_nconn, micro_xnodes, nconn_CFD,enrich_nconn, xnodes_CFD,...
                     slaveMicroElemList, h_CFD, nx_CFD, ny_CFD, dNdx_Macro_at_microIP, enrich_surfconnD, enrich_surfconnR, enrich_surfconnU,... 
                     enrich_surfconnL, vnodes, v_beam, cp, rho, kappa, latent, nnEnrichment, nEnrichment, LHS_micro, RHS_micro, RHS_micro_latent);
                  
%% Solve Macro and micro equations iteratively until convergence
tolerance = 1.0e-3;
maxIter = 10;
iter = 1;
delta = 1.0;        % norm_2 of of consecutive Macro solutions

while (delta > tolerance && iter <= maxIter)
    
% Micro:      
                     
    % Update the main LHS/RHS of micro equation
    [rearranged_phi_Macro,rearranged_Lf] = rearrangeMacroPhi(phi_Macro, Lf, Macro_nconn);
    mainRHS_micro = -(RHS_micro * rearranged_phi_Macro) - (RHS_micro_latent * rearranged_Lf);

    % Apply micro boundary conditions
    mainRHS_micro = applyZeroDirichletBC_RHS(enrich_surfconnD, mainRHS_micro, nnEnrichment, nMacro);
    mainRHS_micro = applyZeroDirichletBC_RHS(enrich_surfconnU, mainRHS_micro, nnEnrichment, nMacro);
    mainRHS_micro = applyZeroDirichletBC_RHS(enrich_surfconnR, mainRHS_micro, nnEnrichment, nMacro);
    mainRHS_micro = applyZeroDirichletBC_RHS(enrich_surfconnL, mainRHS_micro, nnEnrichment, nMacro);

    % Solve micro equation 
    phi_micro = LHS_micro \ mainRHS_micro;
    
% Macro:                     
    % Update the main LHS/RHS of Macro equation
    mainRHS_Macro = -(RHS_Macro * phi_micro) - (RHS_Macro_latent * Lf);
    
    % Apply Macro boundary conditions
    [LHS_Macro,mainRHS_Macro] = applyLaserFlux(Macro_xnodes, Macro_surfconnU, LHS_Macro, mainRHS_Macro);
    [LHS_Macro,mainRHS_Macro] = applyDirichlet(Macro_xnodes, Macro_surfconnD, @fixedBcFunction, LHS_Macro, mainRHS_Macro, 295.0);
    [LHS_Macro,mainRHS_Macro] = applyDirichlet(Macro_xnodes, Macro_surfconnR, @fixedBcFunction, LHS_Macro, mainRHS_Macro, 295.0);
    [LHS_Macro,mainRHS_Macro] = applyDirichlet(Macro_xnodes, Macro_surfconnL, @fixedBcFunction, LHS_Macro, mainRHS_Macro, 295.0);
    
    % Solve Macro equation and update liquid fraction 
    phi_Macro = LHS_Macro \ mainRHS_Macro;
    Lf = updateLiquidFraction(phi_Macro, solidus, liquidus, nnMacro);
    
    delPhi = phi_Macro - phi_Macro_old;
    delta = norm(delPhi, 2)/ ((nxMacro-1)*(nyMacro-1));
    X = ['Iteration = ',num2str(iter),'   norm_2 = ', num2str(delta)];
    disp(X)
    phi_Macro_old = phi_Macro;    
    iter = iter + 1;
   
end

%% Plot data

display_microPhi = rearrangeMicroPhi(phi_micro, nEnrichment, nnEnrichment, slaveMicroElemList,...
                                     micro_nconn, enrich_nconn, nMacro, nnMicro);
figure(1);
plotScalar(micro_xnodes,display_microPhi);
pbaspect([6 1 1]);
title('micro \phi''');

figure(2);
plotScalar(Macro_xnodes,phi_Macro);
pbaspect([6 1 1]);
title('VMS Macro \phi');

figure(3);
plotScalar(Macro_xnodes,Lf);
pbaspect([6 1 1]);
title('VMS liquid fraction');



