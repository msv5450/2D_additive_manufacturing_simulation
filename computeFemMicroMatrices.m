function [LHS,RHS,advect_latent] = computeFemMicroMatrices(Macro_nconn, micro_nconn, micro_xnodes, nconn_CFD,enrich_nconn, xnodes_CFD,...
                     slaveMicroElemList, h_CFD, nx_CFD, ny_CFD, dNdx_Macro_at_microIP, enrich_surfconnD, enrich_surfconnR, enrich_surfconnU,... 
                     enrich_surfconnL, vnodes, v_beam, cp, rho, kappa, latent, nnEnrichment, nEnrichment, LHS, RHS, advect_latent)

% compute global LHS and RHS matrices of micro equation 

Macro_ne = size(Macro_nconn,1);   % number of Macro elements
Macro_nen = size(Macro_nconn,2);  % number of Macro nodes per Macro element
ndim = 2;                         % number of spatial dimensions (2 for 2D)
micro_nen = size(micro_nconn,2);  % number of micro nodes per Macro element
micro_nq  = 4;                    % number of micro element integration points
CFD_nen = size(nconn_CFD,2);      % number of nodes per NALU element

% Loop over Macro elements
for MacroEID = 1:Macro_ne

    % LHS Macro element matrix
    LHSe_Macro = zeros(nnEnrichment,nnEnrichment);      
    % RHS Macro element matrix
    RHSe_Macro = zeros(nnEnrichment,Macro_nen);
    RHSe_advect_latent = zeros(nnEnrichment,Macro_nen); 
    offSet = (MacroEID-1) * nnEnrichment;
    
    % get eIDs of all micro elemets inside this enriched Macro element
    microElemIDs = slaveMicroElemList(MacroEID,:);
    
    % loop over micro elements that are inside this Macro element
    for local_microEID = 1:nEnrichment
        
        % LHS micro element matrices
        LHSeDif_micro = zeros(micro_nen,micro_nen);     % Difffusion part
        LHSeAdv_micro = zeros(micro_nen,micro_nen);     % Advection part
        % RHS micro element matrices
        RHSeDif_micro = zeros(micro_nen,Macro_nen);     % Difffusion part
        RHSeAdv_micro = zeros(micro_nen,Macro_nen);     % Advection part
        RHSe_advect_latent_micro = zeros(micro_nen,Macro_nen);
        
        global_microEID = microElemIDs(local_microEID);
        % Coordinates for micro element nodes (size: 4X2)
        coords_micro = micro_xnodes(micro_nconn(global_microEID,:)',:);
        
        % get micro shape functions at micro integration points
        [xq_micro, wq_micro, n_micro, dndx_micro] = computeQuad2dMicroShapeFunctions(coords_micro);

        % loop over micro integration points 
        for iq = 1:micro_nq
            
            % coordinates of integration point
            xy_ip_micro = xq_micro(iq,:);
        
            % Find which NALU element this Macro integration point belongs to
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

%% Form LHSe
            % Loop over micro node pairs
            for A = 1:micro_nen 
                for B = 1:micro_nen 

                    % Diffusion
                    for idim = 1:ndim   
                        LHSeDif_micro(A,B) = ...
                            LHSeDif_micro(A,B) +...
                            dndx_micro(iq,A,idim) * kappa * dndx_micro(iq,B,idim) * wq_micro(iq); 
                    end

                    %Advection
                    for idim = 1:ndim  
                        LHSeAdv_micro(A,B) = ...
                            LHSeAdv_micro(A,B) +...
                            n_micro(iq,A) * rho * cp * u_microIp(idim) * dndx_micro(iq,B,idim) * wq_micro(iq);
                    end
                end
            end
%% Form RHSe
            % Loop over micro nodes
            for A = 1:micro_nen 
                % Loop over Macro nodes from master Macro element
                for I = 1:Macro_nen 

                    % Diffusion
                    for idim = 1:ndim   
                        RHSeDif_micro(A,I) = ...
                            RHSeDif_micro(A,I) +...
                            dndx_micro(iq,A,idim) * kappa * dNdx_Macro_at_microIP(MacroEID,iq,I,local_microEID,idim) * wq_micro(iq); 
                    end

                    %Advection of Macro temp and latent heat
                    for idim = 1:ndim  
                        RHSeAdv_micro(A,I) = ...
                            RHSeAdv_micro(A,I) +...
                            n_micro(iq,A) * rho * cp * u_microIp(idim) * dNdx_Macro_at_microIP(MacroEID,iq,I,local_microEID,idim) * wq_micro(iq);
                        
                        RHSe_advect_latent_micro(A,I) = ...
                            RHSe_advect_latent_micro(A,I) +...
                            n_micro(iq,A) * rho * latent * u_microIp(idim) * dNdx_Macro_at_microIP(MacroEID,iq,I,local_microEID,idim) * wq_micro(iq);
                    end
                end
            end
        end
%% Assemble local micro LHSe/RHSe to Macro matrices

        % Assemble LHSe_micro to LHSe_Macro 
        for i = 1:micro_nen                    
            I = enrich_nconn(local_microEID,i);
            for j = 1:micro_nen
                J = enrich_nconn(local_microEID,j);
                LHSe_Macro(I,J) = LHSe_Macro(I,J) + LHSeDif_micro(i,j) + LHSeAdv_micro(i,j); 
            end
        end

        % Assemble RHSe_micro to RHSe_Macro  
        for i = 1:micro_nen                    
            I = enrich_nconn(local_microEID,i);
            for j = 1:Macro_nen
                J = j;
                RHSe_Macro(I,J) = RHSe_Macro(I,J) + RHSeDif_micro(i,j) + RHSeAdv_micro(i,j); 
                RHSe_advect_latent(I,J) = RHSe_advect_latent(I,J) + RHSe_advect_latent_micro(i,j);
            end
        end
        
        % Apply zero Dirichlet BC on LHSe of micro equation 
        LHSe_Macro = applyZeroDirichletBC_LHS(enrich_surfconnD, LHSe_Macro);
        LHSe_Macro = applyZeroDirichletBC_LHS(enrich_surfconnL, LHSe_Macro);
        LHSe_Macro = applyZeroDirichletBC_LHS(enrich_surfconnR, LHSe_Macro);
        LHSe_Macro = applyZeroDirichletBC_LHS(enrich_surfconnU, LHSe_Macro);
    end
%% Assemble Macro LHSe/RHSe to global micro matrices

    % Assemble LHSe_Macro to global micro LHS matrix 
    for i = 1:nnEnrichment                    
        I = offSet + i;
        for j = 1:nnEnrichment
            J = offSet + j;
            LHS(I,J) = LHS(I,J) + LHSe_Macro(i,j); 
        end
    end
    
    % Assemble RHSe_Macro to global micro RHS matrix
    for i = 1:nnEnrichment                    
        I = offSet + i;
        for j = 1:Macro_nen
            J = (MacroEID - 1) * Macro_nen + j;
            RHS(I,J) = RHS(I,J) + RHSe_Macro(i,j); 
            advect_latent(I,J) = advect_latent(I,J) + RHSe_advect_latent(i,j);
        end
   end
end