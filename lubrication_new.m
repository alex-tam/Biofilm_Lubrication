function lubrication_new
%--------------------------------------------------------------------------
% LUBRICATION_NEW Solve the lubrication model.
%   Solves the lubrication model for biofilm expansion.
%   Alex Tam, 13/03/2019.
%--------------------------------------------------------------------------
%-------------------------- Model assumptions -----------------------------
epsilon = 0.1; % [-] biofilm aspect ratio
H0 = 0.1; % [-] maximum initial biofilm height
Psi_m = 1/9; % [-] ECM production rate
Psi_d = 0; % [-] cell death rate
precursor = 1e-4; % [-] precursor film height
threshold = 0.002; %1.1*precursor; % [-] source term threshold

%------------------------ Dimensional parameters --------------------------
R_b = 2.875; % [mm] Initial biofilm size
R_p = 41.5; % [mm] Petri dish radius
D_s = (1 - 0.023*0.3)*4.04e-2; % [mm^2/min] Diffusivity of glucose in substratum
D_b = 0.24*4.04e-2; % [mm^2/min] Diffusivity of glucose in biofilm
Q = 2.92e-3; % [mm/min] Mass transfer coefficient
G = 0.5/(pi*R_p^2); % [g/mm^2] Initial nutrient concentration (0.5g)

%-------------------------- Fitted parameters -----------------------------
psi_n = 50; % [mm^2/g/min] cell production rate
eta = 3.7e-3; % [/min] nutrient consumption rate

%----------------------- Dimensionless parameters -------------------------
R_dim = R_p/R_b; % [-] dimensionless Petri dish radius
Upsilon = eta*R_b^2/D_b; % [-] nutrient consumption rate
Qs = Q*R_b/(epsilon*D_s); % [-] substratum permeability
Qb = Q*R_b/(epsilon*D_b); % [-] biofilm permeability
D = D_s/(psi_n*G*R_b^2); % [-] nutrient diffusion coefficient
Pe = (psi_n*G*R_b^2)/D_b; % [-] Peclet number
gamma = 10; % [-] surface tension coefficient
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

% Indices for the 2D advection equation
Index = reshape(1:nR*nZ, nR, nZ); % Matrix of each index
% Construct the seven regions
rl_z = reshape(Index(1, 1:nZ), nZ, 1); % r = 0, z
rl_z_rp = reshape(Index(2, 1:nZ), nZ, 1);
rl_z_rpp = reshape(Index(3, 1:nZ), nZ, 1);
ri_zi = reshape(Index(jXi, jZeta), (nR-2)*(nZ-2), 1); % interior r, interior z
ri_zi_rp = reshape(Index(jXi+1, jZeta), (nR-2)*(nZ-2), 1);
ri_zi_rn = reshape(Index(jXi-1, jZeta), (nR-2)*(nZ-2), 1);
ri_zi_zp = reshape(Index(jXi, jZeta+1), (nR-2)*(nZ-2), 1);
ri_zi_zn = reshape(Index(jXi, jZeta-1), (nR-2)*(nZ-2), 1);
rr_zi = reshape(Index(nR, jZeta), nZ-2, 1); % r = R, interior z
rr_zi_rn = reshape(Index(nR-1, jZeta), nZ-2, 1);
rr_zi_rnn = reshape(Index(nR-2, jZeta), nZ-2, 1);
rr_zi_zp = reshape(Index(nR, jZeta+1), nZ-2, 1);
rr_zi_zn = reshape(Index(nR, jZeta-1), nZ-2, 1);
ri_zb = reshape(Index(jXi, 1), nR-2, 1); % interior r, z = 0
ri_zb_rp = reshape(Index(jXi-1, 1), nR-2, 1);
ri_zb_rn = reshape(Index(jXi-1, 1), nR-2, 1);
ri_zb_zp = reshape(Index(jXi, 2), nR-2, 1);
ri_zb_zpp = reshape(Index(jXi, 3), nR-2, 1);
ri_zt = reshape(Index(jXi, nZ), nR-2, 1); % interior r, z = h
ri_zt_rp = reshape(Index(jXi+1, nZ), nR-2, 1);
ri_zt_rn = reshape(Index(jXi-1, nZ), nR-2, 1);
ri_zt_zn = reshape(Index(jXi, nZ-1), nR-2, 1);
ri_zt_znn = reshape(Index(jXi, nZ-2), nR-2, 1);
rr_zb = Index(nR, 1); % r = R, z = 0
rr_zb_rn = Index(nR-1, 1);
rr_zb_rnn = Index(nR-2, 1);
rr_zb_zp = Index(nR, 2);
rr_zb_zpp = Index(nR, 3);
rr_zt = Index(nR, nZ); % r = R, z = h
rr_zt_rn = Index(nR-1, nZ);
rr_zt_rnn = Index(nR-2, nZ);
rr_zt_zn = Index(nR, nZ-1);
rr_zt_znn = Index(nR, nZ-2);

% Indices of entries in the sparse matrix    
Index2D_1 = [ rl_z ; rl_z ; rl_z ; ...
    ri_zb ; ri_zb ; ri_zb ; ri_zb ; ri_zb ; ...
    ri_zi ; ri_zi ; ri_zi ; ri_zi ; ri_zi ; ...
    rr_zi ; rr_zi ; rr_zi ; rr_zi ; rr_zi ; ...
    ri_zt ; ri_zt ; ri_zt ; ri_zt ; ri_zt ; ...
    rr_zt ; rr_zt ; rr_zt ; rr_zt ; rr_zt ; ...
    rr_zb ; rr_zb ; rr_zb ; rr_zb ; rr_zb ];
Index2D_2 = [ rl_z ; rl_z_rp ; rl_z_rpp ; ...
    ri_zb ; ri_zb_rp ; ri_zb_rn ; ri_zb_zp ; ri_zb_zpp ; ...
    ri_zi ; ri_zi_rp ; ri_zi_rn ; ri_zi_zp ; ri_zi_zn ; ...
    rr_zi ; rr_zi_rn ; rr_zi_rnn ; rr_zi_zp ; rr_zi_zn ; ...
    ri_zt ; ri_zt_rp ; ri_zt_rn ; ri_zt_zn ; ri_zt_znn ; ...
    rr_zt ; rr_zt_rn ; rr_zt_rnn ; rr_zt_zn ; rr_zt_znn ; ...
    rr_zb ; rr_zb_rn ; rr_zb_rnn ; rr_zb_zp ; rr_zb_zpp ];

%-------------------------- Initial conditons -----------------------------
H = precursor + (H0-precursor)*(ones(size(r)) - (r/1).^2).*(r <= 1);  % [-] initial biofilm height (parabolic)
% H = (precursor + (H0-precursor)*(1 - (r/1).^2).^4).*(r <= 1) + precursor*(r > 1); % [-] initial biofilm height (Ward/King)
contact_line = r(find(H <= threshold, 1))*ones(nTimes,1); % [-] initial contact line position
thickness_index = max(H)/contact_line(1)*ones(nTimes, 1); % [-] initial thickness index
% c = 50; r0 = 0.9; Vol_Frac = ((exp(-c*r0) + 2*exp(-c*R))./(exp(-c*r0) + exp(-c*R)) - 1).*(R <= 1); % [-] Sigmoidal initial condition
% Vol_Frac = (1*ones(nR, nZ) - 3*R.^2 + 2*R.^3).*(R <=1); % [-] smooth step function
Vol_Frac = (3*Zeta.^2 - 2*Zeta.^3).*(1*ones(nR, nZ) - 3*R.^2 + 2*R.^3).*(R <=1); % [-] smooth step function
% Z = Zeta.*repmat(H, 1, nZ);
% dist = sqrt(R.^2 + (H0*ones(nR,nZ) - Z).^2);
% Vol_Frac = (1*ones(nR, nZ) - 3*dist.^2 + 2*dist.^3).*(R <=1); % [-] smooth step function
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
Surface_Tension = surf_tens(H, r, dr, nR, threshold);

H_mat = repmat(H, 1, nZ);
Integral_uz = nan(nR,nZ);
for row = 1:nR
    Integral_uz(row, :) = cumtrapz(Z(row,:), Vol_Frac(row,:));
end
theta_uz = Z.^2/2.*(Z/3 - H_mat).*repmat(Surface_Tension, 1, nZ);
derivative_uz = nan(nR, nZ);
derivative_uz(1,:) = (-3*theta_uz(1,:) + 4*theta_uz(2,:) - theta_uz(3,:))/dr;
derivative_uz(jd,:) = (R(jd+1,:).*theta_uz(jd+1,:) - R(jd-1,:).*theta_uz(jd-1,:))./(2*R(jd,:)*dr);
derivative_uz(nR,:) = (3*R(nR,:).*theta_uz(nR,:) - 4*R(nR-1,:).*theta_uz(nR-1,:) + R(nR-2,:).*theta_uz(nR-2,:))./(2*R(nR,:)*dr);
Uz = (1 + Psi_m)*repmat(Gb, 1, nZ).*Integral_uz.*(H_mat >= threshold) + gamma*derivative_uz.*(H_mat >= threshold);
Uz(1,:) = (4*Uz(2,:) - Uz(3,:))/3; % artificially set derivative to zero

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
    h_mat = repmat(h, 1, nZ);
    gb_mat = repmat(gb, 1, nZ);
    bar_phi_mat = repmat(bar_phi, 1, nZ);
    theta_mat = repmat(surface_tension, 1, nZ);
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
    % Compute d/dR(RH^3T)/R
    ddR = nan(nR, nZ);
    ddR(1,:) = 2*h_mat(1,:).^3.*(-3*theta_mat(1,:) + 4*theta_mat(2,:) - theta_mat(3,:))/(2*dxi); % Replace with L'Hop
    ddR(jXi,:) = (Xi(jXi+1,:).*h_mat(jXi+1,:).^3.*theta_mat(jXi+1,:) - Xi(jXi-1,:).*h_mat(jXi-1,:).^3.*theta_mat(jXi-1,:))./(2*dxi*Xi(jXi,:));
    ddR(nR,:) = (3*Xi(nR,:).*h_mat(nR,:) - 4*Xi(nR-1,:).*h_mat(nR-1,:) + Xi(nR-2,:).*h_mat(nR-2,:))./(2*dxi*Xi(nR,:));
    % Compute dHdt
    dHdt = (1+Psi_m)*bar_phi_mat.*gb_mat.*h_mat - gamma/3*ddR;
    % Compute advection coefficients
    a1 = -gamma*Zeta.*(Zeta/2 - ones(nR, nZ));
    a2 = uz./h_mat - Zeta./h_mat.*dHdt + gamma*Zeta.^2.*(Zeta/2 - ones(nR, nZ)).*h_mat.*theta_mat.*dHdX;
    a3 = gb_mat - Psi_d - dUdZ./h_mat - gamma*Zeta.*(Zeta - ones(nR, nZ)).*h_mat.*theta_mat.*dHdX;
%     a2 = 0.1*Zeta.^2; % artificial
    a1 = a1.*(h_mat >= threshold); a2 = a2.*(h_mat >= threshold); a3 = a3.*(h_mat >= threshold);
    % Solve for volume fraction (Crank--Nicolson)  
    A = lhs_phi(nR, nZ, Index2D_1, Index2D_2, dxi, dzeta, dt, Xi, theta_mat, h_mat, a1, a2, a3, rl_z, ri_zb, ri_zb_rp, ri_zb_rn, ri_zi, ri_zi_rp, ri_zi_rn, rr_zi, rr_zi_rn, rr_zi_rnn, ri_zt, ri_zt_rp, ri_zt_rn, rr_zt, rr_zt_rn, rr_zt_rnn, rr_zb, rr_zb_rn, rr_zb_rnn);
    b = rhs_phi(nR, nZ, jXi, jZeta, dxi, dzeta, dt, Xi, theta_mat, h_mat, volume_fraction, a1, a2, a3);
    [precL, precU] = ilu(A);
    volume_fraction = gmres(A, b, 10, 1e-6, 10, precL, precU, volume_fraction);
    Vol_Frac = reshape(volume_fraction, nR, nZ);
    % 3. Nutrient concentration (substratum)
    Ags = matrix_gs(jd, r, dr, dt, D, Qs, nR, h, threshold); bgs = rhs_gs(jd, r, dr, dt, D, Qs, nR, gs, gb, h, threshold);
    Gs = Ags\bgs;
    % Compute surface tension term: expanded form, sequential differentiation (2nd order)
    Integral_gb = nan(nR,1);
    for row = 1:nR
        Integral_gb(row) = trapz(Z(row,:), Z(row,:).*(Z(row,:)/2 - H(row)).*(1 - vol_frac(row,:)));
    end
    theta_gb = Integral_gb.*surface_tension;
    theta_gb(1) = 0; % enforced by velocity condition
    % 4. Nutrient concentration (biofilm)
    Agb = matrix_gb(h, r, dr, nR, gamma, Pe, Qb, Upsilon, bar_phi, jd, dt, theta_gb, threshold); bgb = rhs_gb(h, r, dr, nR, gamma, Pe, Qb, Upsilon, bar_phi, jd, dt, theta_gb, gb, gs, threshold);
    Gb = Agb\bgb;
    % 5. Depth-averaged cell volume fraction
    Z = Zeta.*repmat(H, 1, nZ);
    Total_Phi = nan(nR, 1);
    for row = 1:nR
        Total_Phi(row) = trapz(Z(row,:), Vol_Frac(row,:));
    end
    Bar_Phi = Total_Phi./H;
    % 6. Vertical velocity
    H_mat = repmat(H, 1, nZ);
    Surface_Tension = surf_tens(H, r, dr, nR, threshold);
%     Integral_uzeta = nan(nR,nZ);
%     for row = 1:nR
%         Integral_uzeta(row, :) = cumtrapz(Zeta(row,:), H(row)*Vol_Frac(row,:));
%     end
%     theta_uzeta = H_mat.^3.*repmat(Surface_Tension, 1, nZ);
%     derivative_uzeta = nan(nR, nZ); % pre-allocate (rz^2/2(z/3-h)Theta)_r/r
%     derivative_uzeta(1,:) = (-3*theta_uzeta(1,:) + 4*theta_uzeta(2,:) - theta_uzeta(3,:))/dxi;
%     derivative_uzeta(jd,:) = (Xi(jd+1,:).*theta_uzeta(jd+1,:) - Xi(jd-1,:).*theta_uzeta(jd-1,:))./(2*Xi(jd,:)*dxi);
%     derivative_uzeta(nR,:) = (3*theta_uzeta(nR,:) - 4*theta_uzeta(nR-1,:) + theta_uzeta(nR-2,:))./(2*Xi(nR,:)*dxi);
%     % Compute dHdX
%     dHdX = nan(nR, nZ);
%     dHdX(1,:) = 0;
%     dHdX(jXi,:) = (H_mat(jXi+1,:)- H_mat(jXi-1,:))./(2*dxi);
%     dHdX(nR,:) = (3*H_mat(nR,:) - 4*H_mat(nR-1,:) + H_mat(nR-2,:))./(2*dxi);
%     % Compute Uz
%     Uz = (1 + Psi_m)*repmat(Gb, 1, nZ).*Integral_uz.*(H_mat >= threshold) ...
%         + gamma/2.*Zeta.^2.*(Zeta/3 - ones(nR, nZ)).*derivative_uzeta.*(H_mat >= threshold)...
%         - Xi.*Zeta.*(Zeta/2 - ones(nR, nZ)).*H_mat.^2.*dHdX.*repmat(Surface_Tension, 1, nZ); 
    Integral_uz = nan(nR,nZ);
    for row = 1:nR
        Integral_uz(row, :) = cumtrapz(Z(row,:), Vol_Frac(row,:));
    end
    theta_uz = Z.^2/2.*(Z/3 - H_mat).*repmat(Surface_Tension, 1, nZ);
    derivative_uz = nan(nR, nZ); % pre-allocate (rz^2/2(z/3-h)Theta)_r/r
    derivative_uz(1,:) = (-3*theta_uz(1,:) + 4*theta_uz(2,:) - theta_uz(3,:))/dr;
    derivative_uz(jd,:) = (R(jd+1,:).*theta_uz(jd+1,:) - R(jd-1,:).*theta_uz(jd-1,:))./(2*R(jd,:)*dr);
    derivative_uz(nR,:) = (3*R(nR,:).*theta_uz(nR,:) - 4*R(nR-1,:).*theta_uz(nR-1,:) + R(nR-2,:).*theta_uz(nR-2,:))./(2*R(nR,:)*dr);
    Uz = (1 + Psi_m)*repmat(Gb, 1, nZ).*Integral_uz.*(H_mat >= threshold) + gamma*derivative_uz.*(H_mat >= threshold);
    Uz(1,:) = (4*Uz(2,:) - Uz(3,:))/3; % artificially set derivative to zero
    % Store contact line position
    contact_line(i+1) = r(find(H <= threshold, 1));
    thickness_index(i+1) = max(H)/contact_line(i+1);
end

dlmwrite(['biofilm_height-',num2str(i+1),'.csv'], H, 'precision', '%.20f');
dlmwrite(['biofilm_nutrient-',num2str(i+1),'.csv'], Gb, 'precision', '%.20f');
dlmwrite(['substratum_nutrient-',num2str(i+1),'.csv'], Gs, 'precision', '%.20f');
dlmwrite(['depth-averaged_vol_frac-',num2str(i+1),'.csv'], Bar_Phi, 'precision', '%.20f');
dlmwrite(['u_z-',num2str(i+1),'.csv'], reshape(Uz, nR*nZ,1), 'precision', '%.5f');
dlmwrite(['vol_frac-',num2str(i+1),'.csv'], reshape(Vol_Frac, nR*nZ,1), 'precision', '%.5f');
dlmwrite('contact_line.csv', contact_line, 'precision', '%.5f');
dlmwrite('thickness_index.csv', thickness_index, 'precision', '%.5f');
       
%% Functions

%--------------------------- Surface tension ------------------------------
function theta = surf_tens(H, r, dr, nR, threshold)
%     % Second-order sequential
%     H_r = [ (-3*H(1) + 4*H(2) - H(3))/(2*dr) ; (H(3:nR) - H(1:nR-2))/(2*dr) ; (3*H(nR) - 4*H(nR-1) + H(nR-2))/(2*dr) ];
%     H_rr = [ (-3*H_r(1) + 4*H_r(2) - H_r(3))/(2*dr) ; (H_r(3:nR) - H_r(1:nR-2))/(2*dr) ; (3*H_r(nR) - 4*H_r(nR-1) + H_r(nR-2))/(2*dr) ];
%     H_rrr = [ (-3*H_rr(1) + 4*H_rr(2) - H_rr(3))/(2*dr) ; (H_rr(3:nR) - H_rr(1:nR-2))/(2*dr) ; (3*H_rr(nR) - 4*H_rr(nR-1) + H_rr(nR-2))/(2*dr) ];
%     theta = H_rrr + H_rr./r - H_r./(r.^2);
%     theta(1) = 0;
    % Conservative form
    theta = nan(nR,1);
    theta(1) = 0;
    theta(2) = (r(3)*(H(4) - H(2)) - r(2)*(H(3) - H(1)))/(dr^3*(r(3)+r(2))) - r(2)*(H(3) - H(1))/(r(2)*dr^3);
    theta(3:nR-2) = ((r(5:nR) + r(4:nR-1)).*(H(5:nR) - H(4:nR-1)) - (r(4:nR-1) + r(3:nR-2)).*(H(4:nR-1) - H(3:nR-2)))./(4*r(4:nR-1)*dr^3) ...
        - ((r(3:nR-2) + r(2:nR-3)).*(H(3:nR-2) - H(2:nR-3)) - (r(2:nR-3) + r(1:nR-4)).*(H(2:nR-3) - H(1:nR-4)))./(4*r(2:nR-3)*dr^3);
    theta(nR-1) = (r(nR)*(3*H(nR) - 4*H(nR-1) + H(nR-2)) - r(nR-1)*(H(nR) - H(nR-2)))/(dr^3*(r(nR) + r(nR-1)))...
        + (r(nR-2)*(H(nR-1) - H(nR-3)) - r(nR-1)*(H(nR) - H(nR-2)))/(dr^3*(r(nR-1) + r(nR-2)));
    theta(nR) = 0;
    theta = smooth(r, theta);
    theta = theta.*(H >= threshold);
%     % Conservative form, avoid discontinuity
%     nB = find(H < threshold, 1) - 1;
%     theta = zeros(nR,1);
%     theta(1) = 0;
%     theta(2) = (r(3)*(H(4) - H(2)) - r(2)*(H(3) - H(1)))/(dr^3*(r(3)+r(2))) - r(2)*(H(3) - H(1))/(r(2)*dr^3);
%     theta(3:nB-2) = ((r(5:nB) + r(4:nB-1)).*(H(5:nB) - H(4:nB-1)) - (r(4:nB-1) + r(3:nB-2)).*(H(4:nB-1) - H(3:nB-2)))./(4*r(4:nB-1)*dr^3) ...
%        - ((r(3:nB-2) + r(2:nB-3)).*(H(3:nB-2) - H(2:nB-3)) - (r(2:nB-3) + r(1:nB-4)).*(H(2:nB-3) - H(1:nB-4)))./(4*r(2:nB-3)*dr^3);
%     theta(nB-1) = (r(nB)*(3*H(nB) - 4*H(nB-1) + H(nB-2)) - r(nB-1)*(H(nB) - H(nB-2)))/(dr^3*(r(nB) + r(nB-1)))...
%        + (r(nB-2)*(H(nB-1) - H(nB-3)) - r(nB-1)*(H(nB) - H(nB-2)))/(dr^3*(r(nB-1) + r(nB-2)));
%     theta(nB) = 3/r(nB)*(9*r(nB)*H(nB) - 12*r(nB)*H(nB-1) + 3*r(nB)*H(nB-2) - 4*(r(nB) + r(nB-1))*(H(nB) - H(nB-1)) + r(nB-1)*(H(nB) - H(nB-2)))/(2*dr^3) ...
%         - 2/(r(nB) + r(nB-1))*(3*r(nB)*H(nB) - 4*r(nB)*H(nB) + r(nB)*H(nB-2) - r(nB-1)*(H(nB) - H(nB-2)))/(2*dr^3) ...
%         + 1/(r(nB-1))*((r(nB) + r(nB-1))*(H(nB) - H(nB-1)) + (r(nB-1) + r(nB-2))*(H(nB-1) - H(nB-2)))/(2*dr^3);
%     theta = theta.*(H >= threshold);
%     % Sixth-order sequential
%     H_r = [(-147*H(1) + 360*H(2) - 450*H(3) + 400*H(4) - 225*H(5) + 72*H(6) - 10*H(7))/(60*dr) ; ...
%           (-10*H(1) - 77*H(2) + 150*H(3) - 100*H(4) + 50*H(5) - 15*H(6) + 2*H(7))/(60*dr) ; ...
%           (2*H(1) - 24*H(2) - 35*H(3) + 80*H(4) - 30*H(5) + 8*H(6) - H(7))/(60*dr) ; ...
%           (-H(1:nR-6) + 9*H(2:nR-5) - 45*H(3:nR-4) + 45*H(5:nR-2) - 9*H(6:nR-1) + H(7:nR))/(60*dr) ; ...
%           (-2*H(nR) + 24*H(nR-1) + 35*H(nR-2) - 80*H(nR-3) + 30*H(nR-4) - 8*H(nR-5) + H(nR-6))/(60*dr) ; ...
%           (10*H(nR) + 77*H(nR-1) - 150*H(nR-2) + 100*H(nR-3) - 50*H(nR-4) + 15*H(nR-5) - 2*H(nR-6))/(60*dr) ; ...
%           (147*H(nR) - 360*H(nR-1) + 450*H(nR-2) - 400*H(nR-3) + 225*H(nR-4) - 72*H(nR-5) + 10*H(nR-6))/(60*dr)];
%     H_rr = [(-147*H_r(1) + 360*H_r(2) - 450*H_r(3) + 400*H_r(4) - 225*H_r(5) + 72*H_r(6) - 10*H_r(7))/(60*dr) ; ...
%           (-10*H_r(1) - 77*H_r(2) + 150*H_r(3) - 100*H_r(4) + 50*H_r(5) - 15*H_r(6) + 2*H_r(7))/(60*dr) ; ...
%           (2*H_r(1) - 24*H_r(2) - 35*H_r(3) + 80*H_r(4) - 30*H_r(5) + 8*H_r(6) - H_r(7))/(60*dr) ; ...
%           (-H_r(1:nR-6) + 9*H_r(2:nR-5) - 45*H_r(3:nR-4) + 45*H_r(5:nR-2) - 9*H_r(6:nR-1) + H_r(7:nR))/(60*dr) ; ...
%           (-2*H_r(nR) + 24*H_r(nR-1) + 35*H_r(nR-2) - 80*H_r(nR-3) + 30*H_r(nR-4) - 8*H_r(nR-5) + H_r(nR-6))/(60*dr) ; ...
%           (10*H_r(nR) + 77*H_r(nR-1) - 150*H_r(nR-2) + 100*H_r(nR-3) - 50*H_r(nR-4) + 15*H_r(nR-5) - 2*H_r(nR-6))/(60*dr) ; ...
%           (147*H_r(nR) - 360*H_r(nR-1) + 450*H_r(nR-2) - 400*H_r(nR-3) + 225*H_r(nR-4) - 72*H_r(nR-5) + 10*H_r(nR-6))/(60*dr)];
%     H_rrr = [(-147*H_rr(1) + 360*H_rr(2) - 450*H_rr(3) + 400*H_rr(4) - 225*H_rr(5) + 72*H_rr(6) - 10*H_rr(7))/(60*dr) ; ...
%           (-10*H_rr(1) - 77*H_rr(2) + 150*H_rr(3) - 100*H_rr(4) + 50*H_rr(5) - 15*H_rr(6) + 2*H_rr(7))/(60*dr) ; ...
%           (2*H_rr(1) - 24*H_rr(2) - 35*H_rr(3) + 80*H_rr(4) - 30*H_rr(5) + 8*H_rr(6) - H_rr(7))/(60*dr) ; ...
%           (-H_rr(1:nR-6) + 9*H_rr(2:nR-5) - 45*H_rr(3:nR-4) + 45*H_rr(5:nR-2) - 9*H_rr(6:nR-1) + H_rr(7:nR))/(60*dr) ; ...
%           (-2*H_rr(nR) + 24*H_rr(nR-1) + 35*H_rr(nR-2) - 80*H_rr(nR-3) + 30*H_rr(nR-4) - 8*H_rr(nR-5) + H_rr(nR-6))/(60*dr) ; ...
%           (10*H_rr(nR) + 77*H_rr(nR-1) - 150*H_rr(nR-2) + 100*H_rr(nR-3) - 50*H_rr(nR-4) + 15*H_rr(nR-5) - 2*H_rr(nR-6))/(60*dr) ; ...
%           (147*H_rr(nR) - 360*H_rr(nR-1) + 450*H_rr(nR-2) - 400*H_rr(nR-3) + 225*H_rr(nR-4) - 72*H_rr(nR-5) + 10*H_rr(nR-6))/(60*dr)];
%     theta = H_rrr + H_rr./r - H_r./(r.^2);
%     theta(1) = 0;
%     theta = theta.*(H >= threshold);
%     % Sixth-order sequential, avoid discontinuity
%     theta = zeros(nR,1);
%     nB = find(H < threshold, 1) - 1;
%     H_r = [(-147*H(1) + 360*H(2) - 450*H(3) + 400*H(4) - 225*H(5) + 72*H(6) - 10*H(7))/(60*dr) ; ...
%           (-10*H(1) - 77*H(2) + 150*H(3) - 100*H(4) + 50*H(5) - 15*H(6) + 2*H(7))/(60*dr) ; ...
%           (2*H(1) - 24*H(2) - 35*H(3) + 80*H(4) - 30*H(5) + 8*H(6) - H(7))/(60*dr) ; ...
%           (-H(1:nB-6) + 9*H(2:nB-5) - 45*H(3:nB-4) + 45*H(5:nB-2) - 9*H(6:nB-1) + H(7:nB))/(60*dr) ; ...
%           (-2*H(nB) + 24*H(nB-1) + 35*H(nB-2) - 80*H(nB-3) + 30*H(nB-4) - 8*H(nB-5) + H(nB-6))/(60*dr) ; ...
%           (10*H(nB) + 77*H(nB-1) - 150*H(nB-2) + 100*H(nB-3) - 50*H(nB-4) + 15*H(nB-5) - 2*H(nB-6))/(60*dr) ; ...
%           (147*H(nB) - 360*H(nB-1) + 450*H(nB-2) - 400*H(nB-3) + 225*H(nB-4) - 72*H(nB-5) + 10*H(nB-6))/(60*dr)];
%     H_rr = [(-147*H_r(1) + 360*H_r(2) - 450*H_r(3) + 400*H_r(4) - 225*H_r(5) + 72*H_r(6) - 10*H_r(7))/(60*dr) ; ...
%           (-10*H_r(1) - 77*H_r(2) + 150*H_r(3) - 100*H_r(4) + 50*H_r(5) - 15*H_r(6) + 2*H_r(7))/(60*dr) ; ...
%           (2*H_r(1) - 24*H_r(2) - 35*H_r(3) + 80*H_r(4) - 30*H_r(5) + 8*H_r(6) - H_r(7))/(60*dr) ; ...
%           (-H_r(1:nB-6) + 9*H_r(2:nB-5) - 45*H_r(3:nB-4) + 45*H_r(5:nB-2) - 9*H_r(6:nB-1) + H_r(7:nB))/(60*dr) ; ...
%           (-2*H_r(nB) + 24*H_r(nB-1) + 35*H_r(nB-2) - 80*H_r(nB-3) + 30*H_r(nB-4) - 8*H_r(nB-5) + H_r(nB-6))/(60*dr) ; ...
%           (10*H_r(nB) + 77*H_r(nB-1) - 150*H_r(nB-2) + 100*H_r(nB-3) - 50*H_r(nB-4) + 15*H_r(nB-5) - 2*H_r(nB-6))/(60*dr) ; ...
%           (147*H_r(nB) - 360*H_r(nB-1) + 450*H_r(nB-2) - 400*H_r(nB-3) + 225*H_r(nB-4) - 72*H_r(nB-5) + 10*H_r(nB-6))/(60*dr)];
%     H_rrr = [(-147*H_rr(1) + 360*H_rr(2) - 450*H_rr(3) + 400*H_rr(4) - 225*H_rr(5) + 72*H_rr(6) - 10*H_rr(7))/(60*dr) ; ...
%           (-10*H_rr(1) - 77*H_rr(2) + 150*H_rr(3) - 100*H_rr(4) + 50*H_rr(5) - 15*H_rr(6) + 2*H_rr(7))/(60*dr) ; ...
%           (2*H_rr(1) - 24*H_rr(2) - 35*H_rr(3) + 80*H_rr(4) - 30*H_rr(5) + 8*H_rr(6) - H_rr(7))/(60*dr) ; ...
%           (-H_rr(1:nB-6) + 9*H_rr(2:nB-5) - 45*H_rr(3:nB-4) + 45*H_rr(5:nB-2) - 9*H_rr(6:nB-1) + H_rr(7:nB))/(60*dr) ; ...
%           (-2*H_rr(nB) + 24*H_rr(nB-1) + 35*H_rr(nB-2) - 80*H_rr(nB-3) + 30*H_rr(nB-4) - 8*H_rr(nB-5) + H_rr(nB-6))/(60*dr) ; ...
%           (10*H_rr(nB) + 77*H_rr(nB-1) - 150*H_rr(nB-2) + 100*H_rr(nB-3) - 50*H_rr(nB-4) + 15*H_rr(nB-5) - 2*H_rr(nB-6))/(60*dr) ; ...
%           (147*H_rr(nB) - 360*H_rr(nB-1) + 450*H_rr(nB-2) - 400*H_rr(nB-3) + 225*H_rr(nB-4) - 72*H_rr(nB-5) + 10*H_rr(nB-6))/(60*dr)];
%     theta(1:nB) = H_rrr + H_rr./r(1:nB) - H_r./(r(1:nB).^2);
%     theta(1) = 0;
end

%----------------------------- Source term --------------------------------
function J = source(Psi_m, Bar_Phi, Gb, h, threshold)
    J = (1 + Psi_m)*Bar_Phi.*Gb.*(h > threshold); % modified source term
end

%------------------- Linear system for biofilm height ---------------------
function A = matrix_h(j, r, nR, dt, dr, gamma, h)
    % common coefficients for interior points
    cp = ((r(j+1) + r(j)).*(h(j+1) + h(j)).^3)./(16*r(j)*dr^2); cn = ((r(j) + r(j-1)).*(h(j) + h(j-1)).^3)./(16*r(j)*dr^2);
    rpp = r(j+2) + r(j+1); rp = r(j+1) + r(j); rn = r(j) + r(j-1); rnn = r(j-1) + r(j-2);
    dp = 2*r(j+1)*dr^2; d = 2*r(j)*dr^2; dn = 2*r(j-1)*dr^2;
    % coefficients for N-1 boundary
    Cp = ((r(nR) + r(nR-1))*(h(nR) + h(nR-1))^3)/(16*r(nR-1)*dr^2); Cn = ((r(nR-1) + r(nR-2))*(h(nR-1) + h(nR-2))^3)/(16*r(nR-1)*dr^2);
    Rpp = 2*r(nR) + dr; Rp = r(nR) + r(nR-1); Rn = r(nR-1) + r(nR-2); Rnn = r(nR-2) + r(nR-3);
    Dp = 2*r(nR)*dr^2; Dm = 2*r(nR-1)*dr^2; Dn = 2*r(nR-2)*dr^2;
    % construct matrix for linear system
    index1 = [ 1 ; 1 ; 1 ; 2 ; 2 ; 2 ; 2 ; 2 ; j ; j ; j ; j ; j ; nR-1 ; nR-1 ; nR-1 ; nR-1 ; nR ];
    index2 = [ 1 ; 2 ; 3 ; 1 ; 2 ; 3 ; 4 ; 5 ; j+2 ; j+1 ; j ; j-1 ; j-2 ; nR-3 ; nR-2 ; nR-1 ; nR ; nR ];
    entries = [ -3/(2*dr) ; 4/(2*dr) ; -1/(2*dr) ; ... % no-flux condition
       -5/(2*dr^3) ; 18/(2*dr^3) ; -24/(2*dr^3) ; 14/(2*dr^3) ; -3/(2*dr^3) ; ... % third derivative condition
       gamma*dt/6*cp.*rpp./dp ; ...
       gamma*dt/6*(-cp.*rpp./dp - cp.*rp./dp - cp.*rp./d - cn.*rp./d) ; ...
       1 + gamma*dt/6*(cp.*rp./dp + cp.*rp./d + cp.*rn./d + cn.*rp./d + cn.*rn./d + cn.*rn./dn) ; ...
       gamma*dt/6*(-cp.*rn./d - cn.*rn./d - cn.*rn./dn -cn.*rnn./dn) ; gamma*dt/6*cn.*rnn./dn ; ... % interior grid points
       gamma*dt/6*(Cn*Rnn/Dn) ; gamma*dt/6*(-Cp*Rn/Dm -Cn*Rn/Dm - Cn*Rn/Dn - Cn*Rnn/Dn) ; ...
       1 + gamma*dt/6*(Cp*Rpp/Dp + Cp*Rp/Dp + Cp*Rp/Dm + Cp*Rn/Dm + Cn*Rp/Dm + Cn*Rn/Dm + Cn*Rn/Dn) ; ...
       gamma*dt/6*(-Cp*Rpp/Dp - Cp*Rp/Dp - Cp*Rp/Dm - Cn*Rp/Dm) ; ... % ghost points for no-flux condition
       1 ]; % Dirichlet condition
    A = sparse(index1, index2, entries, nR, nR);
end

function b = rhs_h(j, r, nR, precursor, J, dt, dr, gamma, h)
    % pre-allocate
    b = nan(nR,1);
    % common coefficients for interior points
    cp = ((r(j+1) + r(j)).*(h(j+1) + h(j)).^3)./(16*r(j)*dr^2); cn = ((r(j) + r(j-1)).*(h(j) + h(j-1)).^3)./(16*r(j)*dr^2);
    rpp = r(j+2) + r(j+1); rp = r(j+1) + r(j); rn = r(j) + r(j-1); rnn = r(j-1) + r(j-2);
    dp = 2*r(j+1)*dr^2; d = 2*r(j)*dr^2; dn = 2*r(j-1)*dr^2;
    % coefficients for N-1 boundary
    Cp = ((r(nR) + r(nR-1))*(h(nR) + h(nR-1))^3)/(16*r(nR-1)*dr^2); Cn = ((r(nR-1) + r(nR-2))*(h(nR-1) + h(nR-2))^3)/(16*r(nR-1)*dr^2);
    Rpp = 2*r(nR) + dr; Rp = r(nR) + r(nR-1); Rn = r(nR-1) + r(nR-2); Rnn = r(nR-2) + r(nR-3);
    Dp = 2*r(nR)*dr^2; Dm = 2*r(nR-1)*dr^2; Dn = 2*r(nR-2)*dr^2;
    % construct RHS
    b(1) = 0; % no-flux condition
    b(2) = 0; % third derivative condition
    b(j) = h(j) + dt*J(j) ...
         + (-gamma*dt/6*cp.*rpp./dp).*h(j+2)...
         + (-gamma*dt/6*(-cp.*rpp./dp - cp.*rp./dp - cp.*rp./d - cn.*rp./d)).*h(j+1)...
         + (-gamma*dt/6*(cp.*rp./dp + cp.*rp./d + cp.*rn./d + cn.*rp./d + cn.*rn./d + cn.*rn./dn)).*h(j)...
         + (-gamma*dt/6*(-cp.*rn./d - cn.*rn./d - cn.*rn./dn -cn.*rnn./dn)).*h(j-1)...
         + (-gamma*dt/6*cn.*rnn./dn).*h(j-2);
    b(nR-1) = h(nR-1) + dt*J(nR-1) ...
                 + (-gamma*dt/6*(-Cp*Rpp/Dp - Cp*Rp/Dp - Cp*Rp/Dm - Cn*Rp/Dm))*h(nR)...
                 + (-gamma*dt/6*(Cp*Rpp/Dp + Cp*Rp/Dp + Cp*Rp/Dm + Cp*Rn/Dm + Cn*Rp/Dm + Cn*Rn/Dm + Cn*Rn/Dn))*h(nR-1)...
                 + (-gamma*dt/6*(-Cp*Rn/Dm - Cn*Rn/Dm - Cn*Rn/Dn - Cn*Rnn/Dn))*h(nR-2)...
                 + (-gamma*dt/6*(Cn*Rnn/Dn))*h(nR-3); % ghost points for no-flux condition
    b(nR) = precursor; % Dirichlet condition
end

%--------------- Linear system for cell volume fraction -------------------
function A = lhs_phi(nR, nZ, Index2D_1, Index2D_2, dxi, dzeta, dt, Xi, theta_mat, h_mat, a1, a2, a3, rl_z, ri_zb, ri_zb_rp, ri_zb_rn, ri_zi, ri_zi_rp, ri_zi_rn, rr_zi, rr_zi_rn, rr_zi_rnn, ri_zt, ri_zt_rp, ri_zt_rn, rr_zt, rr_zt_rn, rr_zt_rnn, rr_zb, rr_zb_rn, rr_zb_rnn)
    % reshape
    a1_vec = reshape(a1, nR*nZ, 1); a2_vec = reshape(a2, nR*nZ, 1); a3_vec = reshape(a3, nR*nZ, 1);
    xi_vec = reshape(Xi, nR*nZ, 1); theta_vec = reshape(theta_mat, nR*nZ, 1); h_vec = reshape(h_mat, nR*nZ,1);
    % Matrix entries
    Entries2D = [ -3/(2*dxi)*ones(size(rl_z)) ; 4/(2*dxi)*ones(size(rl_z)) ; -1/(2*dxi)*ones(size(rl_z)) ; ...
        1 - dt/2*a3_vec(ri_zb) - 3*dt/(4*dzeta)*a2_vec(ri_zb) ; dt/(4*dxi)*a1_vec(ri_zb)./xi_vec(ri_zb).*xi_vec(ri_zb_rp).*h_vec(ri_zb_rp).^2.*theta_vec(ri_zb_rp) ; -dt/(4*dxi)*a1_vec(ri_zb)./xi_vec(ri_zb).*xi_vec(ri_zb_rn).*h_vec(ri_zb_rn).^2.*theta_vec(ri_zb_rn) ; 4*dt/(4*dzeta)*a2_vec(ri_zb) ; -dt/(4*dzeta)*a2_vec(ri_zb) ; ...
        1 - dt/2*a3_vec(ri_zi) ; dt/(4*dxi)*a1_vec(ri_zi)./xi_vec(ri_zi).*xi_vec(ri_zi_rp).*h_vec(ri_zi_rp).^2.*theta_vec(ri_zi_rp) ; -dt/(4*dxi)*a1_vec(ri_zi)./xi_vec(ri_zi).*xi_vec(ri_zi_rn).*h_vec(ri_zi_rn).^2.*theta_vec(ri_zi_rn) ; dt/(4*dzeta)*a2_vec(ri_zi) ; -dt/(4*dzeta)*a2_vec(ri_zi) ; ...
        1 - dt/2*a3_vec(rr_zi) + 3*dt/(4*dxi)*a1_vec(rr_zi)./xi_vec(rr_zi).*xi_vec(rr_zi).*h_vec(rr_zi).^2.*theta_vec(rr_zi) ; -4*dt/(4*dxi)*a1_vec(rr_zi)./xi_vec(rr_zi).*xi_vec(rr_zi_rn).*h_vec(rr_zi_rn).^2.*theta_vec(rr_zi_rn) ; dt/(4*dxi)*a1_vec(rr_zi)./xi_vec(rr_zi).*xi_vec(rr_zi_rnn).*h_vec(rr_zi_rnn).^2.*theta_vec(rr_zi_rnn) ; dt/(4*dzeta)*a2_vec(rr_zi) ; -dt/(4*dzeta)*a2_vec(rr_zi) ; ...
        1 - dt/2*a3_vec(ri_zt) + 3*dt/(4*dzeta)*a2_vec(ri_zt) ; dt/(4*dxi)*a1_vec(ri_zt)./xi_vec(ri_zt).*xi_vec(ri_zt_rp).*h_vec(ri_zt_rp).^2.*theta_vec(ri_zt_rp) ; -dt/(4*dxi)*a1_vec(ri_zt)./xi_vec(ri_zt).*xi_vec(ri_zt_rn).*h_vec(ri_zt_rn).^2.*theta_vec(ri_zt_rn) ; -4*dt/(4*dzeta)*a2_vec(ri_zt) ; dt/(4*dzeta)*a2_vec(ri_zt) ; ...
        1 - dt/2*a3_vec(rr_zt) + 3*dt/(4*dxi)*a1_vec(rr_zt)./xi_vec(rr_zt).*xi_vec(rr_zt).*h_vec(rr_zt).^2.*theta_vec(rr_zt) + 3*dt/(4*dzeta)*a2_vec(rr_zt) ; -4*dt/(4*dxi)*a1_vec(rr_zt)./xi_vec(rr_zt).*xi_vec(rr_zt_rn).*h_vec(rr_zt_rn).^2.*theta_vec(rr_zt_rn) ; dt/(4*dxi)*a1_vec(rr_zt)./xi_vec(rr_zt).*xi_vec(rr_zt_rnn).*h_vec(rr_zt_rnn).^2.*theta_vec(rr_zt_rnn) ; -4*dt/(4*dzeta)*a2_vec(rr_zt) ; dt/(4*dzeta)*a2_vec(rr_zt) ; ...
        1 - dt/2*a3_vec(rr_zb) + 3*dt/(4*dxi)*a1_vec(rr_zb)./xi_vec(rr_zb).*xi_vec(rr_zb).*h_vec(rr_zb).^2.*theta_vec(rr_zb) - 3*dt/(4*dzeta)*a2_vec(rr_zb) ; -4*dt/(4*dxi)*a1_vec(rr_zb)./xi_vec(rr_zb).*xi_vec(rr_zb_rn).*h_vec(rr_zb_rn).^2.*theta_vec(rr_zb_rn) ; dt/(4*dxi)*a1_vec(rr_zb)./xi_vec(rr_zb).*xi_vec(rr_zb_rnn).*h_vec(rr_zb_rnn).^2.*theta_vec(rr_zb_rnn) ; 4*dt/(4*dzeta)*a2_vec(rr_zb) ; -dt/(4*dzeta)*a2_vec(rr_zb) ];
    A = sparse(Index2D_1, Index2D_2, Entries2D, nR*nZ, nR*nZ);
%     % Construct indices for non-corner points
%     A = sparse(nR*nZ, nR*nZ);
%     % Symmetry condition for phi_r
%     A = A + sparse(index_r_in, index_r_in, -3/(2*dxi), nR*nZ, nR*nZ)...
%         + sparse(index_r_in, index_r_in_p, 4/(2*dxi), nR*nZ, nR*nZ)...
%         + sparse(index_r_in, index_r_in_pp, -1/(2*dxi), nR*nZ, nR*nZ);
%     % Boundary condition for phi_z
%     A = A + sparse(index_z_in, index_z_in, -3/(2*dzeta), nR*nZ, nR*nZ)...
%         + sparse(index_z_in, index_z_in_p, 4/(2*dzeta), nR*nZ, nR*nZ)...
%         + sparse(index_z_in, index_z_in_pp, -1/(2*dzeta), nR*nZ, nR*nZ);
%     % Interior grid points
%     A = A + sparse(index_int, index_int, 1 - dt/2*a3_vec(index_int), nR*nZ, nR*nZ)...
%         + sparse(index_int, index_int_rp, dt/(4*dxi)*a1_vec(index_int)./xi_vec(index_int).*xi_vec(index_int_rp).*h_vec(index_int_rp).^2.*theta_vec(index_int_rp), nR*nZ, nR*nZ)...
%         + sparse(index_int, index_int_rn, -dt/(4*dxi)*a1_vec(index_int)./xi_vec(index_int).*xi_vec(index_int_rn).*h_vec(index_int_rn).^2.*theta_vec(index_int_rn), nR*nZ, nR*nZ)...
%         + sparse(index_int, index_int_zp, dt/(4*dzeta)*a2_vec(index_int), nR*nZ, nR*nZ)...
%         + sparse(index_int, index_int_zn, -dt/(4*dzeta)*a2_vec(index_int), nR*nZ, nR*nZ);
%     % Interior in z, boundary in r
%     A = A + sparse(index_r_out, index_r_out, 1 - dt/2*a3_vec(index_r_out) + 3*dt/(4*dxi)*a1_vec(index_r_out)./xi_vec(index_r_out).*xi_vec(index_r_out).*h_vec(index_r_out).^2.*theta_vec(index_r_out), nR*nZ, nR*nZ)...
%         + sparse(index_r_out, index_r_out_rn, -4*dt/(4*dxi)*a1_vec(index_r_out)./xi_vec(index_r_out).*xi_vec(index_r_out_rn).*h_vec(index_r_out_rn).^2.*theta_vec(index_r_out_rn), nR*nZ, nR*nZ)...
%         + sparse(index_r_out, index_r_out_rnn, dt/(4*dxi)*a1_vec(index_r_out)./xi_vec(index_r_out).*xi_vec(index_r_out_rnn).*h_vec(index_r_out_rnn).^2.*theta_vec(index_r_out_rnn), nR*nZ, nR*nZ)...
%         + sparse(index_r_out, index_r_out_zp, dt/(4*dzeta)*a2_vec(index_r_out), nR*nZ, nR*nZ)...
%         + sparse(index_r_out, index_r_out_zn, -dt/(4*dzeta)*a2_vec(index_r_out), nR*nZ, nR*nZ);
%     % Interior in r, boundary in z
%     A = A + sparse(index_z_out, index_z_out, 1 - dt/2*a3_vec(index_z_out) + 3*dt/(4*dzeta)*a2_vec(index_z_out), nR*nZ, nR*nZ)...
%         + sparse(index_z_out, index_z_out_rp, dt/(4*dxi)*a1_vec(index_z_out)./xi_vec(index_z_out).*xi_vec(index_z_out_rp).*h_vec(index_z_out_rp).^2.*theta_vec(index_z_out_rp), nR*nZ, nR*nZ)...
%         + sparse(index_z_out, index_z_out_rn, -dt/(4*dxi)*a1_vec(index_z_out)./xi_vec(index_z_out).*xi_vec(index_z_out_rn).*h_vec(index_z_out_rn).^2.*theta_vec(index_z_out_rn), nR*nZ, nR*nZ)...
%         + sparse(index_z_out, index_z_out_zn, -4*dt/(4*dzeta)*a2_vec(index_z_out), nR*nZ, nR*nZ)...
%         + sparse(index_z_out, index_z_out_znn, dt/(4*dzeta)*a2_vec(index_z_out), nR*nZ, nR*nZ);
%     % Outer corner
%     A = A + sparse(index_corner, index_corner, 1 - dt/2*a3_vec(index_corner) + 3*dt/(4*dxi)*a1_vec(index_corner)./xi_vec(index_corner).*xi_vec(index_corner).*h_vec(index_corner).^2.*theta_vec(index_corner) + 3*dt/(4*dzeta)*a2_vec(index_corner), nR*nZ, nR*nZ)...
%         + sparse(index_corner, index_corner_rn, -4*dt/(4*dxi)*a1_vec(index_corner)./xi_vec(index_corner).*xi_vec(index_corner_rn).*h_vec(index_corner_rn).^2.*theta_vec(index_corner_rn), nR*nZ, nR*nZ)...
%         + sparse(index_corner, index_corner_rnn, dt/(4*dxi)*a1_vec(index_corner)./xi_vec(index_corner).*xi_vec(index_corner_rnn).*h_vec(index_corner_rnn).^2.*theta_vec(index_corner_rnn), nR*nZ, nR*nZ)...
%         + sparse(index_corner, index_corner_zn, -4*dt/(4*dzeta)*a2_vec(index_corner), nR*nZ, nR*nZ)...
%         + sparse(index_corner, index_corner_znn, dt/(4*dzeta)*a2_vec(index_corner), nR*nZ, nR*nZ);
end

function b = rhs_phi(nR, nZ, jXi, jZeta, dxi, dzeta, dt, Xi, theta_mat, h_mat, volume_fraction, a1, a2, a3)
    Phi = reshape(volume_fraction, nR, nZ);
    dPdR = nan(nR, nZ); dPdZ = nan(nR, nZ);
    % Compute d/dR(XTH^2P)/R
    dPdR(1,:) = 2*h_mat(1,:).^2.*Phi(1,:).*(-3*theta_mat(1,:) + 4*theta_mat(2,:) - theta_mat(3,:))/(2*dxi); % Replace with L'Hop
    dPdR(jXi,:) = (Xi(jXi+1,:).*theta_mat(jXi+1,:).*h_mat(jXi+1,:).^2.*Phi(jXi+1,:) - Xi(jXi-1,:).*theta_mat(jXi-1,:).*h_mat(jXi-1,:).^2.*Phi(jXi-1,:))./(2*dxi*Xi(jXi,:));
    dPdR(nR,:) = (3*Xi(nR,:).*theta_mat(nR,:).*h_mat(nR,:).^2.*Phi(nR,:) - 4*Xi(nR-1,:).*theta_mat(nR-1,:).*h_mat(nR-1,:).^2.*Phi(nR-1,:) + Xi(nR-2,:).*theta_mat(nR-2,:).*h_mat(nR-2,:).^2.*Phi(nR-2,:))./(2*dxi*Xi(nR,:));
    % Compute dP/dZ
    dPdZ(:,1) = (-3*Phi(:,1) + 4*Phi(:,2) - Phi(:,3))/(2*dzeta);
    dPdZ(:,jZeta) = (Phi(:,jZeta+1) - Phi(:,jZeta-1))/(2*dzeta);
    dPdZ(:,nZ) = (3*Phi(:,nZ) - 4*Phi(:,nZ-1) + Phi(:,nZ-2))/(2*dzeta);
    % Construct Ax matrix vector product
    B = Phi - dt/2*a1.*dPdR - dt/2*a2.*dPdZ + dt/2*a3.*Phi; % Crank--Nicolson
    B(1,:) = 0; % replace term with explicit BC at r = 0
    b = reshape(B, nR*nZ, 1);
end

%---------------- Linear system for substratum nutrient -------------------
function A = matrix_gs(jd, r, dr, dt, D, Qs, nR, h, threshold)
    index1 = [ 1 ; 1 ; jd ; jd ; jd ; nR ; nR ];
    index2 = [ 1 ; 2 ; jd+1 ; jd ; jd-1 ; nR-1 ; nR ];
    entries = [ 1 + 2*dt*D/dr^2 + dt*D*Qs/2*(h(1) >= threshold) ; -2*dt*D/dr^2 ; ... % ghost points for no-flux condition
       -dt*D/(4*dr^2)*(r(jd+1) + r(jd))./r(jd) ; ones(nR-2, 1) + dt*D/(dr^2) + dt*D*Qs/2*(h(jd) >= threshold) ; -dt*D/(4*dr^2)*(r(jd) + r(jd-1))./r(jd) ; ... % interior grid points
       -dt*D/(dr^2) ; 1 + dt*D/(dr^2) + dt*D*Qs/2*(h(nR) >= threshold) ]; % ghost points for no-flux condition
    A = sparse(index1, index2, entries, nR, nR);
end

function b = rhs_gs(jd, r, dr, dt, D, Qs, nR, gs, gb, h, threshold)
    b = nan(nR,1);
    b(1) = gs(1)*(1 - 2*dt*D/dr^2 - dt*D*Qs/2*(h(1) >= threshold)) + gs(2)*(2*dt*D/dr^2) + dt*D*Qs*gb(1)*(h(1) >= threshold); % ghost points for no-flux condition
    b(jd) = gs(jd).*(1 - dt*D/(dr^2) - dt*D*Qs/2*(h(jd) >= threshold))...
        + gs(jd+1).*(dt*D/(4*dr^2)*(r(jd+1) + r(jd))./r(jd))...
        + gs(jd-1).*(dt*D/(4*dr^2)*(r(jd) + r(jd-1))./r(jd)) + dt*D*Qs*gb(jd).*(h(jd) >= threshold); % interior grid points
    b(nR) = gs(nR)*(1 - dt*D/(dr^2) - dt*D*Qs/2*(h(nR) >= threshold)) + gs(nR-1)*(dt*D/(dr^2)) + dt*D*Qs*gb(nR)*(h(nR) >= threshold); % ghost points for no-flux condition
end

% Fully regularised
%------------------ Linear system for biofilm nutrient --------------------
function A = matrix_gb(h, r, dr, nR, gamma, Pe, Qb, Upsilon, bar_phi, jd, dt, theta, threshold)
    % Construct linear system
    index1 = [ 1 ; 1 ; jd ; jd ; jd ; nR ; nR ];
    index2 = [ 1 ; 2 ; jd+1 ; jd ; jd-1 ; nR-1 ; nR ];
    entries = [ 1*(h(1) < threshold) + (Pe*h(1) + dt*gamma*Pe*(-3*theta(1) + 4*theta(2) - theta(3))/(6*dr) + 2*dt*h(1)/(dr^2) + dt*Qb/2 + dt*Upsilon*bar_phi(1)*h(1)/2)*(h(1) >= threshold) ; -2*dt*h(1)/(dr^2)*(h(1) >= threshold) ; ... % ghost points for no-flux condition
        (dt*gamma*Pe*(r(jd+1).*theta(jd+1))./(12*dr*r(jd)) - dt*(r(jd) + r(jd+1)).*(h(jd) + h(jd+1))./(8*r(jd)*dr^2)).*(h(jd) >= threshold) ; ...
        ones(nR-2,1).*(h(jd) < threshold) + (Pe*h(jd) + dt*(r(jd) + r(jd+1)).*(h(jd) + h(jd+1))./(8*r(jd)*dr^2) + dt*(r(jd) + r(jd-1)).*(h(jd) + h(jd-1))./(8*r(jd)*dr^2) + dt*Qb/2 + dt*Upsilon*bar_phi(jd).*h(jd)/2).*(h(jd) >= threshold) ; ...
        (-dt*gamma*Pe*(r(jd-1).*theta(jd-1))./(12*dr*r(jd)) - dt*(r(jd) + r(jd-1)).*(h(jd) + h(jd-1))./(8*r(jd)*dr^2)).*(h(jd) >= threshold) ; ... % interior grid points
        -dt*h(nR)/(dr^2)*(h(nR) >= threshold) ; 1*(h(nR) < threshold) + (Pe*h(nR) + dt*gamma*Pe*(3*theta(nR) - 4*theta(nR-1) + theta(nR-2))/(12*r(nR)*dr) + dt*h(nR)/dr^2 + dt*Qb/2 + dt*Upsilon*bar_phi(nR)*h(nR)/2)*(h(nR) >= threshold) ]; % ghost points for no-flux condition
    A = sparse(index1, index2, entries, nR, nR);
end

function b = rhs_gb(h, r, dr, nR, gamma, Pe, Qb, Upsilon, bar_phi, jd, dt, theta, gb, gs, threshold)
    b = nan(nR,1);
    b(1) = gb(1)*(h(1) < threshold) + gb(1)*(Pe*h(1) - dt*gamma*Pe*(-3*theta(1) + 4*theta(2) - theta(3))/(6*dr) - 2*dt*h(1)/(dr^2) - dt*Qb/2 - dt*Upsilon*bar_phi(1)*h(1)/2)*(h(1) >= threshold)...
        + gb(2)*(2*dt*h(1)/(dr^2))*(h(1) >= threshold) + dt*Qb*gs(1)*(h(1) >= threshold); % ghost points for no-flux condition
    b(jd) = gb(jd).*(h(jd) < threshold)  + gb(jd).*(Pe*h(jd) - dt*(r(jd) + r(jd+1)).*(h(jd) + h(jd+1))./(8*r(jd)*dr^2) - dt*(r(jd) + r(jd-1)).*(h(jd) + h(jd-1))./(8*r(jd)*dr^2) - dt*Qb/2 - dt*Upsilon*bar_phi(jd).*h(jd)/2).*(h(jd) >= threshold) ...
        + gb(jd+1).*(-dt*gamma*Pe*(r(jd+1).*theta(jd+1))./(12*dr*r(jd)) + dt*(r(jd) + r(jd+1)).*(h(jd) + h(jd+1))./(8*r(jd)*dr^2)).*(h(jd) >= threshold) ...
        + gb(jd-1).*( dt*gamma*Pe*(r(jd-1).*theta(jd-1))./(12*dr*r(jd)) + dt*(r(jd) + r(jd-1)).*(h(jd) + h(jd-1))./(8*r(jd)*dr^2)).*(h(jd) >= threshold)...
        + dt*Qb*gs(jd).*(h(jd) >= threshold); % interior grid points
    b(nR) = gb(nR)*(h(nR) < threshold) + gb(nR)*(Pe*h(nR) - dt*gamma*Pe*(3*theta(nR) - 4*theta(nR-1) + theta(nR-2))/(12*r(nR)*dr) - dt*h(nR)/dr^2 - dt*Qb/2 - dt*Upsilon*bar_phi(nR)*h(nR)/2)*(h(nR) >= threshold)...
        + gb(nR-1)*(dt*h(nR)/(dr^2))*(h(nR) >= threshold) + dt*Qb*gs(nR)*(h(nR) >= threshold); % ghost points for no-flux condition
end

% % Regularised, retaining diffusion
% %------------------ Linear system for biofilm nutrient --------------------
% function A = matrix_gb(h, r, dr, nR, gamma, Pe, Qb, Upsilon, bar_phi, jd, dt, theta, threshold)
%     % Construct linear system
%     index1 = [ 1 ; 1 ; jd ; jd ; jd ; nR ; nR ];
%     index2 = [ 1 ; 2 ; jd+1 ; jd ; jd-1 ; nR-1 ; nR ];
%     entries = [ Pe*h(1) + dt*gamma*Pe*(-3*theta(1) + 4*theta(2) - theta(3))/(6*dr)*(h(1) >= threshold) + 2*dt*h(1)/(dr^2) + dt*Qb/2*(h(1) >= threshold) + dt*Upsilon*bar_phi(1)*h(1)/2*(h(1) >= threshold) ; -2*dt*h(1)/(dr^2) ; ... % ghost points for no-flux condition
%         dt*gamma*Pe*(r(jd+1).*theta(jd+1))./(12*dr*r(jd)).*(h(jd) >= threshold) - dt*(r(jd) + r(jd+1)).*(h(jd) + h(jd+1))./(8*r(jd)*dr^2) ; ...
%         Pe*h(jd) + dt*(r(jd) + r(jd+1)).*(h(jd) + h(jd+1))./(8*r(jd)*dr^2) + dt*(r(jd) + r(jd-1)).*(h(jd) + h(jd-1))./(8*r(jd)*dr^2) + dt*Qb/2*(h(jd) >= threshold) + dt*Upsilon*bar_phi(jd).*h(jd)/2.*(h(jd) >= threshold) ; ...
%         -dt*gamma*Pe*(r(jd-1).*theta(jd-1))./(12*dr*r(jd)).*(h(jd) >= threshold) - dt*(r(jd) + r(jd-1)).*(h(jd) + h(jd-1))./(8*r(jd)*dr^2) ; ... % interior grid points
%         -dt*h(nR)/(dr^2) ; Pe*h(nR) + dt*gamma*Pe*(3*theta(nR) - 4*theta(nR-1) + theta(nR-2))/(12*r(nR)*dr)*(h(nR) >= threshold) + dt*h(nR)/dr^2 + dt*Qb/2*(h(nR) >= threshold) + dt*Upsilon*bar_phi(nR)*h(nR)/2*(h(nR) >= threshold) ]; % ghost points for no-flux condition
%     A = sparse(index1, index2, entries, nR, nR);
% end
% 
% function b = rhs_gb(h, r, dr, nR, gamma, Pe, Qb, Upsilon, bar_phi, jd, dt, theta, gb, gs, threshold)
%     b = nan(nR,1);
%     b(1) = gb(1)*(Pe*h(1) - dt*gamma*Pe*(-3*theta(1) + 4*theta(2) - theta(3))/(6*dr)*(h(1) >= threshold) - 2*dt*h(1)/(dr^2) - dt*Qb/2*(h(1) >= threshold) - dt*Upsilon*bar_phi(1)*h(1)/2*(h(1) >= threshold))...
%         + gb(2)*(2*dt*h(1)/(dr^2)) + dt*Qb*gs(1)*(h(1) >= threshold); % ghost points for no-flux condition
%     b(jd) = gb(jd).*(Pe*h(jd) - dt*(r(jd) + r(jd+1)).*(h(jd) + h(jd+1))./(8*r(jd)*dr^2) - dt*(r(jd) + r(jd-1)).*(h(jd) + h(jd-1))./(8*r(jd)*dr^2) - dt*Qb/2.*(h(jd) >= threshold) - dt*Upsilon*bar_phi(jd).*h(jd)/2.*(h(jd) >= threshold)) ...
%         + gb(jd+1).*(-dt*gamma*Pe*(r(jd+1).*theta(jd+1))./(12*dr*r(jd)).*(h(jd) >= threshold) + dt*(r(jd) + r(jd+1)).*(h(jd) + h(jd+1))./(8*r(jd)*dr^2)) ...
%         + gb(jd-1).*( dt*gamma*Pe*(r(jd-1).*theta(jd-1))./(12*dr*r(jd)).*(h(jd) >= threshold) + dt*(r(jd) + r(jd-1)).*(h(jd) + h(jd-1))./(8*r(jd)*dr^2))...
%         + dt*Qb*gs(jd).*(h(jd) >= threshold); % interior grid points
%     b(nR) = gb(nR)*(Pe*h(nR) - dt*gamma*Pe*(3*theta(nR) - 4*theta(nR-1) + theta(nR-2))/(12*r(nR)*dr)*(h(nR) >= threshold) - dt*h(nR)/dr^2 - dt*Qb/2*(h(nR) >= threshold) - dt*Upsilon*bar_phi(nR)*h(nR)/2*(h(nR) >= threshold))...
%         + gb(nR-1)*(dt*h(nR)/(dr^2)) + dt*Qb*gs(nR)*(h(nR) >= threshold); % ghost points for no-flux condition
% end

end