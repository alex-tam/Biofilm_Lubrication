function lubrication
%--------------------------------------------------------------------------
% LUBRICATION Solve the lubrication model.
%   Solves the lubrication model for biofilm expansion.
%   Alex Tam, 13/03/2019.
%--------------------------------------------------------------------------
%-------------------------- Model assumptions -----------------------------
epsilon = 0.1; % [-] biofilm aspect ratio
H0 = 0.1; % [-] maximum initial biofilm height
Psi_m = 1/9; % [-] ECM production rate
Psi_d = 0; % [-] cell death rate
precursor = 1e-4; % [-] precursor film height
threshold = 1.1*precursor; % [-] source term threshold

%------------------------ Dimensional parameters --------------------------
R_b = 2.875; % [mm] Initial biofilm size
R_p = 41.5; % [mm] Petri dish radius
D_s = (1 - 0.023*0.3)*4.04e-2; % [mm^2/min] Diffusivity of glucose in substratum
D_b = 0.24*4.04e-2; % [mm^2/min] Diffusivity of glucose in biofilm
Q = 2.92e-3; % [mm/min] Mass transfer coefficient
G = 0.5/(pi*R_p^2); % [g/mm^2] Initial nutrient concentration (0.5g)

%-------------------------- Fitted parameters -----------------------------
psi_n = 8; % [mm^2/g/min] cell production rate
eta = 1e-2; % [/min] nutrient consumption rate

%----------------------- Dimensionless parameters -------------------------
R_dim = R_p/R_b; % [-] dimensionless Petri dish radius
Upsilon = eta*R_b^2/D_b; % [-] nutrient consumption rate
Qs = Q*R_b/(epsilon*D_s); % [-] substratum permeability
Qb = Q*R_b/(epsilon*D_b); % [-] biofilm permeability
D = D_s/(psi_n*G*R_b^2); % [-] nutrient diffusion coefficient
Pe = (psi_n*G*R_b^2)/D_b; % [-] Peclet number
gamma = 1; % [-] surface tension coefficient
T = 14208*psi_n*G; % [-] dimensionless time in the experiment

%------------------------- Numerical parameters ---------------------------
nR = 201; j = (3:nR-2)'; jd = (2:nR-1)';
nZ = 101;
nTimes = 10001;
r = linspace(0, R_dim, nR)'; dr = r(2) - r(1);
R = repmat(r, 1, nZ); 
Xi = R; jXi = (2:nR-1)';
jZeta = (2:nZ-1)';
Zeta = repmat(linspace(0, 1, nZ), nR, 1);
dxi = dr;
dzeta = Zeta(1,2) - Zeta(1,1);

t = linspace(0, T, nTimes); dt = t(2) - t(1);
output_files = 10; % number of files to output
dlmwrite('r.csv', r, 'precision', 20);
dlmwrite('nR.csv', nR); dlmwrite('nZ.csv', nZ);
dlmwrite('dish_size.csv', R_dim, 'precision', 20);
dlmwrite('t.csv', t, 'precision', 20);
dlmwrite('R_mat.csv', reshape(R, nR*nZ,1));
dlmwrite('Zeta_mat.csv', reshape(Zeta, nR*nZ,1));

%-------------------------- Initial conditons -----------------------------
H = precursor + (H0-precursor)*(ones(size(r)) - r.^2).*(r < 1);  % [-] initial biofilm height
Vol_Frac = 1*ones(nR, nZ); % [-] initial cell volume fraction
volume_fraction = reshape(Vol_Frac, nR*nZ, 1);
Gs = ones(nR, 1); % [-] initial substratum nutrient concentration within the biofilm
Gb = zeros(nR, 1); % [-] initial biofilm nutrient concentration

% Solve for depth-averaged cell volume fraction
Z = Zeta.*repmat(H, 1, nZ);
Total_Phi = nan(nR,1);
for row = 1:nR
    Total_Phi(row) = trapz(Z(row,:), Vol_Frac(row,:));
end
Bar_Phi = Total_Phi./H;

% Solve for velocity
H_r = [ (-3*H(1) + 4*H(2) - H(3))/(2*dr) ; (H(3:nR) - H(1:nR-2))/(2*dr) ; (3*H(nR) - 4*H(nR-1) + H(nR-2))/(2*dr) ];
H_rr = [ (-3*H_r(1) + 4*H_r(2) - H_r(3))/(2*dr) ; (H_r(3:nR) - H_r(1:nR-2))/(2*dr) ; (3*H_r(nR) - 4*H_r(nR-1) + H_r(nR-2))/(2*dr) ];
H_rrr = [ (-3*H_rr(1) + 4*H_rr(2) - H_rr(3))/(2*dr) ; (H_rr(3:nR) - H_rr(1:nR-2))/(2*dr) ; (3*H_rr(nR) - 4*H_rr(nR-1) + H_rr(nR-2))/(2*dr) ];
Surface_Tension = H_rrr + H_rr./r - H_r./(r.^2);
Surface_Tension(1) = 0;

Integral_uz = nan(nR,nZ);
for row = 1:nR
    Integral_uz(row, :) = cumtrapz(Z(row,:), Vol_Frac(row,:));
end
theta_uz = Z.^2/2.*(Z/3 - repmat(H, 1, nZ)).*repmat(Surface_Tension, 1, nZ);
theta_uz(1,:) = 0;
derivative_uz = nan(nR, nZ);
derivative_uz(1,:) = (-3*theta_uz(1,:) + 4*theta_uz(2,:) - theta_uz(3,:))/dr;
derivative_uz(jd,:) = (R(jd+1,:).*theta_uz(jd+1,:) - R(jd-1,:).*theta_uz(jd-1,:))./(2*R(jd,:)*dr);
derivative_uz(nR,:) = (3*R(nR,:).*theta_uz(nR,:) - 4*R(nR-1,:).*theta_uz(nR-1,:) + R(nR-2,:).*theta_uz(nR-2,:))./(2*R(nR,:)*dr);
Uz = (1 + Psi_m)*repmat(Gb, 1, nZ).*Integral_uz + gamma*derivative_uz;

%------------------------------ Solve PDEs --------------------------------
for i = 1:nTimes-1
    % Extract data
    if sum(i == 1:(nTimes-1)/output_files:nTimes) == 1
        dlmwrite(['biofilm_height-',num2str(i),'.csv'], H, 'precision', '%.20f');
        dlmwrite(['depth-averaged_vol_frac-',num2str(i),'.csv'], Bar_Phi, 'precision', '%.20f');
        dlmwrite(['biofilm_nutrient-',num2str(i),'.csv'], Gb, 'precision', '%.20f');
        dlmwrite(['substratum_nutrient-',num2str(i),'.csv'], Gs, 'precision', '%.20f');
        dlmwrite(['u_z-',num2str(i),'.csv'], Uz, 'precision', '%.5f');
        dlmwrite(['vol_frac-',num2str(i),'.csv'], Vol_Frac, 'precision', '%.5f');
    end
    % Store variables
    h = H; gs = Gs; gb = Gb; bar_phi = Bar_Phi; surface_tension = Surface_Tension; vol_frac = Vol_Frac; uz = Uz;
    J = source(Psi_m, bar_phi, gb, h, threshold);
    % 1. Biofilm height
    Ah = matrix_h(j, r, nR, dt, dr, gamma, h); bh = rhs_h(j, r, nR, precursor, J, dt, dr, gamma, h);
    H = Ah\bh;
    % 2. Cell volume fraction
    % Construct matrices
    H_mat = repmat(h, 1, nZ);
    h_mat = repmat(h, 1, nZ);
%     bar_phi_mat = repmat(bar_phi, 1, nZ);
    theta_mat = repmat(Surface_Tension, 1, nZ);
    % Compute dHdt
    dHdt = (H_mat - h_mat)/dt;
    % Compute dUdZ
    dUdZ = nan(nR, nZ);
    dUdZ(:,1) = (-3*uz(:,1) + 4*uz(:,2) - uz(:,3))/(2*dzeta);
    dUdZ(:,jZeta) = (uz(:,jZeta+1) - uz(:,jZeta-1))/(2*dzeta);
    dUdZ(:,nZ) = (3*uz(:,nZ) - 4*uz(:,nZ-1) + uz(:,nZ-2))/(2*dzeta);
    % Compute dHdR
    dHdX = nan(nR, nZ);
    dHdX(1,:) = 0;
    dHdX(jXi,:) = (h_mat(jXi+1,:)- h_mat(jXi-1,:))./(2*dxi);
    dHdX(nR,:) = (3*h_mat(nR,:) - 4*h_mat(nR-1,:) + h_mat(nR-2,:))./(2*dxi);
%     % Compute d/dR(RH^3T)/R
%     ddR = nan(nR, nZ);
%     ddR(1,:) = 2*h_mat(1,:).^3.*(-3*theta_mat(1,:) + 4*theta_mat(2,:) - theta_mat(3,:))/(2*dxi); % Replace with L'Hop
%     ddR(jXi,:) = (Xi(jXi+1,:).*h_mat(jXi+1,:).^3.*theta_mat(jXi+1,:) - Xi(jXi-1,:).*h_mat(jXi-1,:).^3.*theta_mat(jXi-1,:))./(2*dxi*Xi(jXi,:));
%     ddR(nR,:) = (3*h_mat(nR,:) - 4*h_mat(nR-1,:) + h_mat(nR-2,:))./(2*dxi*Xi(nR,:));
    % Compute advection coefficients
    a1 = -gamma*Zeta.*(Zeta/2 - ones(nR, nZ));
    a2 = uz./h_mat - Zeta./h_mat.*dHdt + gamma*Zeta.^2.*(Zeta/2 - ones(nR, nZ)).*h_mat.*theta_mat.*dHdX;
    a3 = Gb - Psi_d - dUdZ./h_mat - gamma*Zeta.*(Zeta - ones(nR, nZ)).*h_mat.*theta_mat.*dHdX;
    % Solve for volume fraction
    b = rhs(nR, nZ, jXi, jZeta, dxi, dzeta, dt, Xi, theta_mat, h_mat, volume_fraction, a1, a2, a3);
    volume_fraction = gmres(@(volume_fraction) lhs(nR, nZ, jXi, jZeta, dxi, dzeta, dt, Xi, theta_mat, h_mat, volume_fraction, a1, a2, a3), b, 10, 1e-8, 10, [], [], volume_fraction);
    Vol_Frac = reshape(volume_fraction, nR, nZ);
    % 3. Nutrient concentration (substratum)
    Ags = matrix_gs(jd, r, dr, dt, D, Qs, nR); bgs = rhs_gs(jd, r, dr, dt, D, Qs, nR, gs, gb);
    Gs = Ags\bgs;
    % Compute surface tension term: expanded form, sequential differentiation (2nd order)
    Integral_gb = nan(nR,1);
    for row = 1:nR
        Integral_gb(row) = trapz(Z(row,:), Z(row,:).*(Z(row,:)/2 - H(row)).*(1 - vol_frac(row,:)));
    end
    theta_gb = Integral_gb.*(surface_tension);
    theta_gb(1) = 0; % enforced by velocity condition
    % 4. Nutrient concentration (biofilm)
    Agb = matrix_gb(h, r, dr, nR, gamma, Pe, Qb, Upsilon, bar_phi, jd, dt, theta_gb); bgb = rhs_gb(h, r, dr, nR, gamma, Pe, Qb, Upsilon, bar_phi, jd, dt, theta_gb, gb, gs);
    Gb = Agb\bgb;
    % 5. Depth-averaged cell volume fraction
    Z = Zeta.*repmat(H, 1, nZ);
    Total_Phi = nan(nR, 1);
    for row = 1:nR
        Total_Phi(row) = trapz(Z(row,:), Vol_Frac(row,:));
    end
    Bar_Phi = Total_Phi./H;
    % 6. Vertical velocity
    H_r = [ (-3*H(1) + 4*H(2) - H(3))/(2*dr) ; (H(3:nR) - H(1:nR-2))/(2*dr) ; (3*H(nR) - 4*H(nR-1) + H(nR-2))/(2*dr) ];
    H_rr = [ (-3*H_r(1) + 4*H_r(2) - H_r(3))/(2*dr) ; (H_r(3:nR) - H_r(1:nR-2))/(2*dr); (3*H_r(nR) - 4*H_r(nR-1) + H_r(nR-2))/(2*dr) ];
    H_rrr = [ (-3*H_rr(1) + 4*H_rr(2) - H_rr(3))/(2*dr) ; (H_rr(3:nR) - H_rr(1:nR-2))/(2*dr); (3*H_rr(nR) - 4*H_rr(nR-1) + H_rr(nR-2))/(2*dr) ];
    Surface_Tension = H_rrr + H_rr./r - H_r./(r.^2);
    Surface_Tension(1) = 0;
    Integral_uz = nan(nR,nZ);
    for row = 1:nR
        Integral_uz(row, :) = cumtrapz(Z(row,:), Vol_Frac(row,:));
    end
    theta_uz = Z.^2/2.*(Z/3 - repmat(H, 1, nZ)).*repmat(Surface_Tension, 1, nZ);
    derivative_uz = nan(nR, nZ);
    derivative_uz(1,:) = (-3*theta_uz(1,:) + 4*theta_uz(2,:) - theta_uz(3,:))/dr;
    derivative_uz(jd,:) = (R(jd+1,:).*theta_uz(jd+1,:) - R(jd-1,:).*theta_uz(jd-1,:))./(2*R(jd,:)*dr);
    derivative_uz(nR,:) = (3*theta_uz(nR,:) - 4*theta_uz(nR-1,:) + theta_uz(nR-2,:))./(2*R(nR,:)*dr);
    Uz = (1 + Psi_m)*repmat(Gb, 1, nZ).*Integral_uz + gamma*derivative_uz;
end

dlmwrite(['biofilm_height-',num2str(i+1),'.csv'], H, 'precision', '%.20f');
dlmwrite(['biofilm_nutrient-',num2str(i+1),'.csv'], Gb, 'precision', '%.20f');
dlmwrite(['substratum_nutrient-',num2str(i+1),'.csv'], Gs, 'precision', '%.20f');
dlmwrite(['depth-averaged_vol_frac-',num2str(i+1),'.csv'], Bar_Phi, 'precision', '%.20f');
dlmwrite(['u_z-',num2str(i+1),'.csv'], reshape(Uz, nR*nZ,1), 'precision', '%.5f');
dlmwrite(['vol_frac-',num2str(i+1),'.csv'], reshape(Vol_Frac, nR*nZ,1), 'precision', '%.5f');
        
%% Functions
%----------------------------- Source term --------------------------------
function J = source(Psi_m, Bar_Phi, Gb, h, threshold)
    J = (1 + Psi_m)*Bar_Phi.*Gb.*(h >= threshold); % modified source term
end

%------------------- Linear system for biofilm height ---------------------
function A = matrix_h(j, r, nPoints, dt, dr, gamma, h)
    % common coefficients for interior points
    cp = ((r(j+1) + r(j)).*(h(j+1) + h(j)).^3)./(16*r(j)*dr^2); cn = ((r(j) + r(j-1)).*(h(j) + h(j-1)).^3)./(16*r(j)*dr^2);
    rpp = r(j+2) + r(j+1); rp = r(j+1) + r(j); rn = r(j) + r(j-1); rnn = r(j-1) + r(j-2);
    dp = 2*r(j+1)*dr^2; d = 2*r(j)*dr^2; dn = 2*r(j-1)*dr^2;
    % coefficients for N-1 boundary
    Cp = ((r(nPoints) + r(nPoints-1))*(h(nPoints) + h(nPoints-1))^3)/(16*r(nPoints-1)*dr^2); Cn = ((r(nPoints-1) + r(nPoints-2))*(h(nPoints-1) + h(nPoints-2))^3)/(16*r(nPoints-1)*dr^2);
    Rpp = 2*r(nPoints) + dr; Rp = r(nPoints) + r(nPoints-1); Rn = r(nPoints-1) + r(nPoints-2); Rnn = r(nPoints-2) + r(nPoints-3);
    Dp = 2*r(nPoints)*dr^2; Dm = 2*r(nPoints-1)*dr^2; Dn = 2*r(nPoints-2)*dr^2;
    % construct matrix for linear system
    index1 = [ 1 ; 1 ; 1 ; 2 ; 2 ; 2 ; 2 ; 2 ; j ; j ; j ; j ; j ; nPoints-1 ; nPoints-1 ; nPoints-1 ; nPoints-1 ; nPoints ];
    index2 = [ 1 ; 2 ; 3 ; 1 ; 2 ; 3 ; 4 ; 5 ; j+2 ; j+1 ; j ; j-1 ; j-2 ; nPoints-3 ; nPoints-2 ; nPoints-1 ; nPoints ; nPoints ];
    entries = [ -3/(2*dr) ; 4/(2*dr) ; -1/(2*dr) ; ... % no-flux condition
       -5/(2*dr^3) ; 18/(2*dr^3) ; -24/(2*dr^3) ; 14/(2*dr^3) ; -3/(2*dr^3) ; ... % third derivative condition
       gamma*dt/6*cp.*rpp./dp ; gamma*dt/6*(-cp.*rpp./dp - cp.*rp./dp - cp.*rp./d - cn.*rp./d) ; 1 + gamma*dt/6*(cp.*rp./dp + cp.*rp./d + cp.*rn./d + cn.*rp./d + cn.*rn./d + cn.*rn./dn) ; gamma*dt/6*(-cp.*rn./d - cn.*rn./d - cn.*rn./dn -cn.*rnn./dn) ; gamma*dt/6*cn.*rnn./dn ; ... % interior grid points
       gamma*dt/6*(Cn*Rnn/Dn) ; gamma*dt/6*(-Cp*Rn/Dm -Cn*Rn/Dm - Cn*Rn/Dn - Cn*Rnn/Dn) ; 1 + gamma*dt/6*(Cp*Rpp/Dp + Cp*Rp/Dp + Cp*Rp/Dm + Cp*Rn/Dm + Cn*Rp/Dm + Cn*Rn/Dm + Cn*Rn/Dn) ; gamma*dt/6*(-Cp*Rpp/Dp - Cp*Rp/Dp - Cp*Rp/Dm - Cn*Rp/Dm) ; ... % ghost points for no-flux condition
       1 ]; % Dirichlet condition
    A = sparse(index1, index2, entries, nPoints, nPoints);
end

function b = rhs_h(j, r, nPoints, precursor, J, dt, dr, gamma, h)
    % pre-allocate
    b = nan(nPoints,1);
    % common coefficients for interior points
    cp = ((r(j+1) + r(j)).*(h(j+1) + h(j)).^3)./(16*r(j)*dr^2); cn = ((r(j) + r(j-1)).*(h(j) + h(j-1)).^3)./(16*r(j)*dr^2);
    rpp = r(j+2) + r(j+1); rp = r(j+1) + r(j); rn = r(j) + r(j-1); rnn = r(j-1) + r(j-2);
    dp = 2*r(j+1)*dr^2; d = 2*r(j)*dr^2; dn = 2*r(j-1)*dr^2;
    % coefficients for N-1 boundary
    Cp = ((r(nPoints) + r(nPoints-1))*(h(nPoints) + h(nPoints-1))^3)/(16*r(nPoints-1)*dr^2); Cn = ((r(nPoints-1) + r(nPoints-2))*(h(nPoints-1) + h(nPoints-2))^3)/(16*r(nPoints-1)*dr^2);
    Rpp = 2*r(nPoints) + dr; Rp = r(nPoints) + r(nPoints-1); Rn = r(nPoints-1) + r(nPoints-2); Rnn = r(nPoints-2) + r(nPoints-3);
    Dp = 2*r(nPoints)*dr^2; Dm = 2*r(nPoints-1)*dr^2; Dn = 2*r(nPoints-2)*dr^2;
    % construct RHS
    b(1) = 0; % no-flux condition
    b(2) = 0; % third derivative condition
    b(j) = h(j) + dt*J(j) ...
         + (-gamma*dt/6*cp.*rpp./dp).*h(j+2)...
         + (-gamma*dt/6*(-cp.*rpp./dp - cp.*rp./dp - cp.*rp./d - cn.*rp./d)).*h(j+1)...
         + (-gamma*dt/6*(cp.*rp./dp + cp.*rp./d + cp.*rn./d + cn.*rp./d + cn.*rn./d + cn.*rn./dn)).*h(j)...
         + (-gamma*dt/6*(-cp.*rn./d - cn.*rn./d - cn.*rn./dn -cn.*rnn./dn)).*h(j-1)...
         + (-gamma*dt/6*cn.*rnn./dn).*h(j-2);
    b(nPoints-1) = h(nPoints-1) + dt*J(nPoints-1) ...
                 + (-gamma*dt/6*(-Cp*Rpp/Dp - Cp*Rp/Dp - Cp*Rp/Dm - Cn*Rp/Dm))*h(nPoints)...
                 + (-gamma*dt/6*(Cp*Rpp/Dp + Cp*Rp/Dp + Cp*Rp/Dm + Cp*Rn/Dm + Cn*Rp/Dm + Cn*Rn/Dm + Cn*Rn/Dn))*h(nPoints-1)...
                 + (-gamma*dt/6*(-Cp*Rn/Dm - Cn*Rn/Dm - Cn*Rn/Dn - Cn*Rnn/Dn))*h(nPoints-2)...
                 + (-gamma*dt/6*(Cn*Rnn/Dn))*h(nPoints-3); % ghost points for no-flux condition
    b(nPoints) = precursor; % Dirichlet condition
end

%--------------- Linear system for cell volume fraction -------------------
function Ax = lhs(nR, nZ, jXi, jZeta, dxi, dzeta, dt, Xi, theta_mat, h_mat, volume_fraction, a1, a2, a3)
    Phi = reshape(volume_fraction, nR, nZ);
    dPdR = sparse(nR, nZ); dPdZ = sparse(nR, nZ);
    % Compute d/dR(XTH^2P)/R
    % dPdR(1,:) = 0;
    dPdR(1,:) = 2*h_mat(1,:).^2.*Phi(1,:).*(-3*theta_mat(1,:) + 4*theta_mat(2,:) - theta_mat(3,:))/(2*dxi); % Replace with L'Hop
    dPdR(jXi,:) = (Xi(jXi+1,:).*theta_mat(jXi+1,:).*h_mat(jXi+1,:).^2.*Phi(jXi+1,:) - Xi(jXi-1,:).*theta_mat(jXi-1,:).*h_mat(jXi-1,:).^2.*Phi(jXi-1,:))./(2*dxi*Xi(jXi,:));
    dPdR(nR,:) = (3*Xi(nR,:).*theta_mat(nR,:).*h_mat(nR,:).^2.*Phi(nR,:) - 4*Xi(nR-1,:).*theta_mat(nR-1,:).*h_mat(nR-1,:).^2.*Phi(nR-1,:) + Xi(nR-2,:).*theta_mat(nR-2,:).*h_mat(nR-2,:).^2.*Phi(nR-2,:))./(2*dxi*Xi(nR,:));
    % Compute dP/dZ
    dPdZ(:,1) = zeros(nR,1);
    dPdZ(:,jZeta) = (Phi(:,jZeta+1) - Phi(:,jZeta-1))/(2*dzeta);
    dPdZ(:,nZ) = (3*Phi(:,nZ) - 4*Phi(:,nZ-1) + Phi(:,nZ-2))/(2*dzeta);
    % Construct Ax matrix vector product
    AX = Phi + dt/2*a1.*dPdR + dt/2*a2.*dPdZ - dt/2*a3.*Phi; % Crank--Nicolson
%     AX = Phi + dt*a1.*dPdR + dt*a2.*dPdZ - dt*a3.*Phi; % Fully implicit
    AX(1,2:nZ) = (-3*Phi(1,2:nZ) + 4*Phi(2,2:nZ) - Phi(3,2:nZ))/(2*dxi); % replace term with explicit BC
    AX(:,1) = (-3*Phi(:,1) + 4*Phi(:,2) - Phi(:,3))/(2*dzeta); % replace term with explicit BC
    Ax = reshape(AX, nR*nZ, 1);
end

function b = rhs(nR, nZ, jXi, jZeta, dxi, dzeta, dt, Xi, theta_mat, h_mat, volume_fraction, a1, a2, a3)
    Phi = reshape(volume_fraction, nR, nZ);
    dPdR = sparse(nR, nZ); dPdZ = sparse(nR, nZ);
    % Compute d/dR(XTH^2P)/R
    % dPdR(1,:) = 0;
    dPdR(1,:) = 2*h_mat(1,:).^2.*Phi(1,:).*(-3*theta_mat(1,:) + 4*theta_mat(2,:) - theta_mat(3,:))/(2*dxi); % Replace with L'Hop
    dPdR(jXi,:) = (Xi(jXi+1,:).*theta_mat(jXi+1,:).*h_mat(jXi+1,:).^2.*Phi(jXi+1,:) - Xi(jXi-1,:).*theta_mat(jXi-1,:).*h_mat(jXi-1,:).^2.*Phi(jXi-1,:))./(2*dxi*Xi(jXi,:));
    dPdR(nR,:) = (3*Xi(nR,:).*theta_mat(nR,:).*h_mat(nR,:).^2.*Phi(nR,:) - 4*Xi(nR-1,:).*theta_mat(nR-1,:).*h_mat(nR-1,:).^2.*Phi(nR-1,:) + Xi(nR-2,:).*theta_mat(nR-2,:).*h_mat(nR-2,:).^2.*Phi(nR-2,:))./(2*dxi*Xi(nR,:));
    % Compute dP/dZ
    dPdZ(:,1) = zeros(nR,1);
    dPdZ(:,jZeta) = (Phi(:,jZeta+1) - Phi(:,jZeta-1))/(2*dzeta);
    dPdZ(:,nZ) = (3*Phi(:,nZ) - 4*Phi(:,nZ-1) + Phi(:,nZ-2))/(2*dzeta);
    % Construct Ax matrix vector product
    B = Phi - dt/2*a1.*dPdR - dt/2*a2.*dPdZ + dt/2*a3.*Phi; % Crank--Nicolson
%     B = Phi; % fully implicit
    B(1,2:nZ) = 0; % replace term with explicit BC
    B(:,1) = 0; % replace term with explicit BC
    b = reshape(B, nR*nZ, 1);
end

%---------------- Linear system for substratum nutrient -------------------
function A = matrix_gs(jd, r, dr, dt, D, Qs, nPoints)
    index1 = [ 1 ; 1 ; jd ; jd ; jd ; nPoints ; nPoints ];
    index2 = [ 1 ; 2 ; jd+1 ; jd ; jd-1 ; nPoints-1 ; nPoints ];
    entries = [ 1 + 2*dt*D/dr^2 + dt*D*Qs/2 ; -2*dt*D/dr^2 ; ... % ghost points for no-flux condition
       -dt*D/(4*dr^2)*(r(jd+1) + r(jd))./r(jd) ; ones(nPoints-2, 1) + dt*D/(dr^2) + dt*D*Qs/2 ; -dt*D/(4*dr^2)*(r(jd) + r(jd-1))./r(jd) ; ... % interior grid points
       -dt*D/(dr^2) ; 1 + dt*D/(dr^2) + dt*D*Qs/2 ]; % ghost points for no-flux condition
    A = sparse(index1, index2, entries, nPoints, nPoints);
end

function b = rhs_gs(jd, r, dr, dt, D, Qs, nPoints, gs, gb)
    b = nan(nPoints,1);
    b(1) = gs(1)*(1 - 2*dt*D/dr^2 - dt*D*Qs/2) + gs(2)*(2*dt*D/dr^2) + dt*D*Qs*gb(1); % ghost points for no-flux condition
    b(jd) = gs(jd).*(1 - dt*D/(dr^2) - dt*D*Qs/2)...
        + gs(jd+1).*(dt*D/(4*dr^2)*(r(jd+1) + r(jd))./r(jd))...
        + gs(jd-1).*(dt*D/(4*dr^2)*(r(jd) + r(jd-1))./r(jd)) + dt*D*Qs*gb(jd); % interior grid points
    b(nPoints) = gs(nPoints)*(1 - dt*D/(dr^2) - dt*D*Qs/2) + gs(nPoints-1)*(dt*D/(dr^2)) + dt*D*Qs*gb(nPoints); % ghost points for no-flux condition
end

%------------------ Linear system for biofilm nutrient --------------------
function A = matrix_gb(h, r, dr, nPoints, gamma, Pe, Qb, Upsilon, bar_phi, jd, dt, theta)
    % Construct linear system
    index1 = [ 1 ; 1 ; jd ; jd ; jd ; nPoints ; nPoints ];
    index2 = [ 1 ; 2 ; jd+1 ; jd ; jd-1 ; nPoints-1 ; nPoints ];
    entries = [ Pe*h(1) + dt*gamma*Pe*(-3*theta(1) + 4*theta(2) - theta(3))/(6*dr) + 2*dt*h(1)/(dr^2) + dt*Qb/2 + dt*Upsilon*bar_phi(1)*h(1)/2 ; -2*dt*h(1)/(dr^2) ; ... % ghost points for no-flux condition
        dt*gamma*Pe*(r(jd+1).*theta(jd+1))./(12*dr*r(jd)) - dt*(r(jd) + r(jd+1)).*(h(jd) + h(jd+1))./(8*r(jd)*dr^2) ; Pe*h(jd) + dt*(r(jd) + r(jd+1)).*(h(jd) + h(jd+1))./(8*r(jd)*dr^2) + dt*(r(jd) + r(jd-1)).*(h(jd) + h(jd-1))./(8*r(jd)*dr^2) + dt*Qb/2 + dt*Upsilon*bar_phi(jd).*h(jd)/2 ; -dt*gamma*Pe*(r(jd-1).*theta(jd-1))./(12*dr*r(jd)) - dt*(r(jd) + r(jd-1)).*(h(jd) + h(jd-1))./(8*r(jd)*dr^2) ; ... % interior grid points
        -dt*h(nPoints)/(dr^2) ; Pe*h(nPoints) + dt*gamma*Pe*(3*theta(nPoints) - 4*theta(nPoints-1) + theta(nPoints-2))/(12*r(nPoints)*dr) + dt*h(nPoints)/dr^2 + dt*Qb/2 + dt*Upsilon*bar_phi(nPoints)*h(nPoints)/2 ]; % ghost points for no-flux condition
    A = sparse(index1, index2, entries, nPoints, nPoints);
end

function b = rhs_gb(h, r, dr, nPoints, gamma, Pe, Qb, Upsilon, bar_phi, jd, dt, theta, gb, gs)
    b = nan(nPoints,1);
    b(1) = gb(1)*(Pe*h(1) - dt*gamma*Pe*(-3*theta(1) + 4*theta(2) - theta(3))/(6*dr) - 2*dt*h(1)/(dr^2) - dt*Qb/2 - dt*Upsilon*bar_phi(1)*h(1)/2)...
        + gb(2)*(2*dt*h(1)/(dr^2)) + dt*Qb*gs(1); % ghost points for no-flux condition
    b(jd) = gb(jd).*(Pe*h(jd) - dt*(r(jd) + r(jd+1)).*(h(jd) + h(jd+1))./(8*r(jd)*dr^2) - dt*(r(jd) + r(jd-1)).*(h(jd) + h(jd-1))./(8*r(jd)*dr^2) - dt*Qb/2 - dt*Upsilon*bar_phi(jd).*h(jd)/2) ...
        + gb(jd+1).*(-dt*gamma*Pe*(r(jd+1).*theta(jd+1))./(12*dr*r(jd)) + dt*(r(jd) + r(jd+1)).*(h(jd) + h(jd+1))./(8*r(jd)*dr^2)) ...
        + gb(jd-1).*( dt*gamma*Pe*(r(jd-1).*theta(jd-1))./(12*dr*r(jd)) + dt*(r(jd) + r(jd-1)).*(h(jd) + h(jd-1))./(8*r(jd)*dr^2))...
        + dt*Qb*gs(jd); % interior grid points
    b(nPoints) = gb(nPoints)*(Pe*h(nPoints) - dt*gamma*Pe*(3*theta(nPoints) - 4*theta(nPoints-1) + theta(nPoints-2))/(12*r(nPoints)*dr) - dt*h(nPoints)/dr^2 - dt*Qb/2 - dt*Upsilon*bar_phi(nPoints)*h(nPoints)/2)...
        + gb(nPoints-1)*(dt*h(nPoints)/(dr^2)) + dt*Qb*gs(nPoints); % ghost points for no-flux condition
end
%--------------------------------------------------------------------------
end