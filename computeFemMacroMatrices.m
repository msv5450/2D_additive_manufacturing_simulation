function [LHS,RHS,advect_latent] = computeFemMacroMatrices(Macro_nconn, micro_nconn, Macro_xnodes, micro_xnodes, nconn_CFD,...
                     enrich_nconn, xnodes_CFD, slaveMicroElemList, N_Macro_at_microIP, dNdx_Macro_at_microIP, nnEnrichment,...
                     nEnrichment, h_CFD, nx_CFD, ny_CFD, vnodes, v_beam, latent, cp, rho, kappa, LHS, RHS, advect_latent)

% compute global LHS and RHS matrices of Macro equation 

Macro_ne = size(Macro_nconn,1);   % number of Macro elements
Macro_nen = size(Macro_nconn,2);  % number of Macro nodes per Macro element
Macro_nq  = 16;                   % number of Macro element integration points
ndim = 2;                         % number of spatial dimensions (2 for 2D)
micro_nen = size(micro_nconn,2);  % number of micro nodes per micro element
micro_nq  = 4;                    % number of micro element integration points
CFD_nen = size(nconn_CFD,2);      % number of nodes per NALU element

% Loop over Macro elements
for MacroEID = 1:Macro_ne
% first form LHS based on only Macro elements
% then form RHS combining micro and Macro
    
%% FORM LHSe
    % LHS Element matrices
    LHSeDif = zeros(Macro_nen,Macro_nen);   % Difffusion part
    LHSeAdv = zeros(Macro_nen,Macro_nen);   % Advection part
    RHSe_Macro = zeros(Macro_nen,nnEnrichment);
    RHSe_advect_LfMacro = zeros(Macro_nen,Macro_nen);  
    
    % Coordinates for element nodes (size: 4X2)
    MacroCoords = Macro_xnodes(Macro_nconn(MacroEID,:)',:);
    
    % get Macro shape functions at the Macro element integration points
    [xq, wq, N, dNdx] = computeQuad2dMacroShapeFunctions(MacroCoords);
    
    % get eIDs of all micro elemets inside this Macro element (1 X nEnrichment)
    microEIDs = slaveMicroElemList(MacroEID,:);
    
    % loop over Macro integration points and form element LHS matrix
    for iq = 1:Macro_nq
        
        % coordinates of integration point
        xy_ip = xq(iq,:);
        
        % Find which NALU element this Macro integration point belongs to
        NALU_eID = locate_ip_inCFDmesh(xy_ip, h_CFD, nx_CFD, ny_CFD);

        % Get the shape functions of the located CFD element at this ip
        CFD_nodeCoords = xnodes_CFD(nconn_CFD(NALU_eID,:)',:);
        N_CFD_ip = computeQuad2dCFDshapeFunctions(xy_ip, CFD_nodeCoords, h_CFD);
        
        % Velocity components for NALU element nodes (size: 4X2)
        u_CFD_nodes = vnodes(nconn_CFD(NALU_eID,:)',:);
        
        % Compute velocity components at this ip Sigma (N_CFD_ip*u_CFD_nodes)
        u_MacroIp = zeros(ndim,1);
        for I = 1:CFD_nen
            for idim = 1:ndim
                u_MacroIp(idim) = u_MacroIp(idim) + N_CFD_ip(I) * u_CFD_nodes(I,idim);
            end
        end
        
        % change the frame of reference by adding scan velocity to V_x
        u_MacroIp(1) = u_MacroIp(1) + v_beam;
        
        % Loop over Macro node pairs
        for i = 1:Macro_nen               
            for j = 1:Macro_nen           
               
                % Diffusion
                for idim = 1:ndim   
                    LHSeDif(i,j) = ...
                        LHSeDif(i,j) +...
                        dNdx(iq,i,idim) * kappa * dNdx(iq,j,idim) * wq(iq); 
                end
       
                %Advection of temp and latent heat
                for idim = 1:ndim  
                    LHSeAdv(i,j) = ...
                        LHSeAdv(i,j) +...
                        N(iq,i) * rho * cp * u_MacroIp(idim) * dNdx(iq,j,idim) * wq(iq);
                    
                    RHSe_advect_LfMacro(i,j) = ...
                        RHSe_advect_LfMacro(i,j) +...
                        N(iq,i) * rho * latent * u_MacroIp(idim) * dNdx(iq,j,idim) * wq(iq);
                end
            end
        end
    end
    
%% Form RHSe
    % loop over micro elements that are inside this Macro element
    for iEnrich=1:nEnrichment
        
        % RHS Element matrices
        RHSeDif_micro = zeros(Macro_nen,micro_nen);   % Difffusion part
        RHSeAdv_micro = zeros(Macro_nen,micro_nen);   % Advection part
        
        microEID = microEIDs(iEnrich);
        % Coordinates for micro element nodes (size: 4X2)
        coords_micro = micro_xnodes(micro_nconn(microEID,:)',:);
        
        % get micro shape functions at micro integration points
        [xq_micro, wq_micro, n_micro, dndx_micro] = computeQuad2dMicroShapeFunctions(coords_micro);
        
        % loop over micro integration points and form element RHS matrix
        for iq = 1:micro_nq
            
            % coordinates of integration point
            xy_ip_micro = xq_micro(iq,:);
        
            % Find which NALU element this micro integration point belongs to
            NALU_eID = locate_ip_inCFDmesh(xy_ip_micro, h_CFD, nx_CFD, ny_CFD);

            % Get the shape functions of the located CFD element at this ip
            CFD_nodeCoords = xnodes_CFD(nconn_CFD(NALU_eID,:)',:);
            N_CFD_ip = computeQuad2dCFDshapeFunctions(xy_ip_micro, CFD_nodeCoords, h_CFD);

            % Velocity components for NALU element nodes (size: 4X2)
            u_CFD_nodes = vnodes(nconn_CFD(NALU_eID,:)',:);

            % Compute velocity components at this ip Sigma (N_CFD_ip*u_CFD_nodes)
            u_microIp = zeros(ndim,1);
            
            for I = 1:CFD_nen
                for idim = 1:ndim
                    u_microIp(idim) = u_microIp(idim) + N_CFD_ip(I) * u_CFD_nodes(I,idim);
                end
            end
            
            % change the frame of reference by adding scan velocity to V_x
            u_microIp(1) = u_microIp(1) + v_beam;
            
            % Loop over Macro nodes
            for I = 1:Macro_nen 
                % Loop over micro nodes
                for A = 1:micro_nen 
                    
                    % Diffusion
                    for idim = 1:ndim   
                        RHSeDif_micro(I,A) = ...
                            RHSeDif_micro(I,A) +...
                            dNdx_Macro_at_microIP(MacroEID,iq,I,iEnrich,idim) * kappa * dndx_micro(iq,A,idim) * wq_micro(iq); 
                    end

                    %Advection
                    for idim = 1:ndim  
                        RHSeAdv_micro(I,A) = ...
                            RHSeAdv_micro(I,A) +...
                            N_Macro_at_microIP(MacroEID,iq,I,iEnrich) * rho * cp * u_microIp(idim) * dndx_micro(iq,A,idim) * wq_micro(iq);
                    end
                end
            end
        end
        
        % Assemble RHSe_micro to RHSe_Macro  
        for i = 1:Macro_nen                    
            I = i;
            for j = 1:micro_nen
                J = enrich_nconn(iEnrich,j);
                RHSe_Macro(I,J) = RHSe_Macro(I,J) + RHSeDif_micro(i,j) + RHSeAdv_micro(i,j); 
            end
        end
    end

%% Assemble Macro LHSe/RHSe to global micro matrices

    % Assemble Macro LHSe to global Macro LHS matrix  
    % Assemble global latent heat advection matrix
    for i = 1:Macro_nen                    
        I = Macro_nconn(MacroEID,i);
        for j = 1:Macro_nen
            J = Macro_nconn(MacroEID,j);
            LHS(I,J) = LHS(I,J) + LHSeDif(i,j) + LHSeAdv(i,j); 
            advect_latent(I,J) = advect_latent(I,J) + RHSe_advect_LfMacro(i,j);
        end
    end
    
    % Assemble RHSe_Macro to global Macro RHS matrix
    for i = 1:Macro_nen                    
        I = Macro_nconn(MacroEID,i);
        for j = 1:nnEnrichment
            J = (MacroEID - 1) * nnEnrichment + j;
            RHS(I,J) = RHS(I,J) + RHSe_Macro(i,j); 
        end
   end
end