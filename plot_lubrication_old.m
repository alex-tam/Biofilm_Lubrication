function plot_lubrication_old
%--------------------------------------------------------------------------
%PLOT_LUBRICATION Plot solutions to lubrication model.
%   Alex Tam, 18/03/2019.
%--------------------------------------------------------------------------
%--------------------------- Import global data ---------------------------
r = importdata('r.csv');
nR = importdata('nR.csv'); nZ = importdata('nZ.csv');
R_dim = importdata('dish_size.csv');
t = importdata('t.csv');
R = importdata('R_mat.csv'); R = reshape(R, nR, nZ);
Zeta = importdata('Zeta_mat.csv'); Zeta = reshape(Zeta, nR, nZ);
contact_line = importdata('contact_line.csv');
thickness_index = importdata('thickness_index.csv');
times = 1:1000:10001;
gamma = 10; threshold = 1.1e-4; Psi_m = 1/9; Psi_d = 0;

%------------------------- Plot solution profiles -------------------------
for plots = 1:length(times)
    % Import solution data
    time_step = times(plots);
    h = importdata(['biofilm_height-',num2str(time_step),'.csv']);   
    gb = importdata(['biofilm_nutrient-',num2str(time_step),'.csv']);
    gs = importdata(['substratum_nutrient-',num2str(time_step),'.csv']);
    bar_phi = importdata(['depth-averaged_vol_frac-',num2str(time_step),'.csv']);
    u_z = importdata(['u_z-',num2str(time_step),'.csv']); Uz = reshape(u_z, nR, nZ); 
    vol_frac = importdata(['vol_frac-',num2str(time_step),'.csv']); Vol_Frac = reshape(vol_frac, nR, nZ);  
    Z = Zeta.*repmat(h, 1, nZ);
    % Plot variable profiles
    hold on
    set(gca, 'FontSize', 16) % change axis tick font size
    plot(r, h, 'LineWidth', 1.5); 
    plot(r, gb, 'LineWidth', 1.5);
    plot(r, gs, 'LineWidth', 1.5);
    plot(r, bar_phi, 'LineWidth', 1.5);
    legend({'\(h\)', '\(g_b\)', '\(g_s\)', '\(\bar{\phi}_n\)'}, 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast')
    xlabel('\(r\)', 'Interpreter', 'latex')
    ylabel('Numerical solution', 'Interpreter', 'latex')
    title(['t = ', num2str(t(times(plots))),'.'], 'Interpreter', 'latex')
    xlim([0, R_dim]); ylim([0 5]);
    print(gcf, '-depsc', ['radial_lub_sol-',num2str(time_step),'.eps'])
    figure
    surf(R, Z, Uz,'EdgeColor','none','LineStyle','none','FaceLighting','phong'); colorbar
    xlabel('\(r\)','FontSize',20,'Interpreter','latex')
    ylabel('\(z\)','FontSize',20,'Interpreter','latex')
    zlabel('\(u_z\)','Interpreter','latex','FontSize',20)
    view(0, 90); 
    print(gcf, '-depsc', ['u_z-',num2str(time_step),'.eps'])
    figure % (az, el)
    surf(R, Z, Vol_Frac,'EdgeColor','none','LineStyle','none','FaceLighting','phong'); colorbar
    xlabel('\(r\)','FontSize',20,'Interpreter','latex')
    ylabel('\(z\)','FontSize',20,'Interpreter','latex')
    zlabel('\(\phi_n\)','Interpreter','latex','FontSize',20)
    view(0, 90) % (az, el)
    figure
end

%----------------------- Plot contact line position -----------------------
format long
fprintf('The final biofilm size is: %f.\n', contact_line(end))
hold on; set(gca, 'FontSize', 16) % change axis tick font size
plot(t, contact_line, 'LineWidth', 1.5); xlim([0, t(end)]); ylim([0 r(end)]);
xlabel('\(t\)', 'Interpreter', 'latex', 'FontSize', 16); ylabel('\(S(t)\)', 'Interpreter', 'latex', 'FontSize', 16); figure

%-------------------------- Plot thickness index --------------------------
format long
hold on; set(gca, 'FontSize', 16) % change axis tick font size
plot(t, thickness_index, 'LineWidth', 1.5); xlim([0, t(end)]); ylim([0 1]);
xlabel('\(t\)', 'Interpreter', 'latex', 'FontSize', 16); ylabel('\(I_t\)', 'Interpreter', 'latex', 'FontSize', 16); figure

%-------------------------- Plot surface tension --------------------------
nPoints = length(r); nR = nPoints; dr = r(2) - r(1);
% h_r = [ (-3*h(1) + 4*h(2) - h(3))/(2*dr) ; (h(3:nPoints) - h(1:nPoints-2))/(2*dr) ; (3*h(nPoints) - 4*h(nPoints-1) + h(nPoints-2))/(2*dr) ];
% h_rr = [ (-3*h_r(1) + 4*h_r(2) - h_r(3))/(2*dr) ; (h_r(3:nPoints) - h_r(1:nPoints-2))/(2*dr); (3*h_r(nPoints) - 4*h_r(nPoints-1) + h_r(nPoints-2))/(2*dr) ];
% h_rrr = [ (-3*h_rr(1) + 4*h_rr(2) - h_rr(3))/(2*dr) ; (h_rr(3:nPoints) - h_rr(1:nPoints-2))/(2*dr); (3*h_rr(nPoints) - 4*h_rr(nPoints-1) + h_rr(nPoints-2))/(2*dr) ];
% theta = (h_rrr + h_rr./r - h_r./(r.^2)); theta(1) = 0; % surface tension term
% theta(1) = 0;

H = h;
theta = nan(nR,1);
theta(1) = 0;
theta(2) = (r(3)*(H(4) - H(2)) - r(2)*(H(3) - H(1)))/(dr^3*(r(3)+r(2))) - r(2)*(H(3) - H(1))/(r(2)*dr^3);
theta(3:nR-2) = ((r(5:nR) + r(4:nR-1)).*(H(5:nR) - H(4:nR-1)) - (r(4:nR-1) + r(3:nR-2)).*(H(4:nR-1) - H(3:nR-2)))./(4*r(4:nR-1)*dr^3) ...
    - ((r(3:nR-2) + r(2:nR-3)).*(H(3:nR-2) - H(2:nR-3)) - (r(2:nR-3) + r(1:nR-4)).*(H(2:nR-3) - H(1:nR-4)))./(4*r(2:nR-3)*dr^3);
theta(nR-1) = (r(nR)*(3*H(nR) - 4*H(nR-1) + H(nR-2)) - r(nR-1)*(H(nR) - H(nR-2)))/(dr^3*(r(nR) + r(nR-1)))...
    + (r(nR-2)*(H(nR-1) - H(nR-3)) - r(nR-1)*(H(nR) - H(nR-2)))/(dr^3*(r(nR-1) + r(nR-2)));
theta(nR) = 0;
theta_r = [ (-3*theta(1) + 4*theta(2) - theta(3))/(2*dr) ; (theta(3:nR) - theta(1:nR-2))/(2*dr) ; (3*theta(nR) - 4*theta(nR-1) + theta(nR-2))/(2*dr) ];

set(gca, 'FontSize', 16) % change axis tick font size
plot(r, h.*theta); xlabel('\(r\)', 'Interpreter', 'latex'); ylabel('\(h\Theta\)', 'Interpreter', 'latex'); figure
plot(r, theta_r); xlabel('\(r\)', 'Interpreter', 'latex'); ylabel('\(\Theta_r\)', 'Interpreter', 'latex'); figure

% jd = (2:nR-1)';
% theta_uz = Z.^2/2.*(Z/3 - repmat(h, 1, nZ)).*repmat(theta, 1, nZ);
% derivative_uz = nan(nR, nZ);
% derivative_uz(1,:) = (-3*theta_uz(1,:) + 4*theta_uz(2,:) - theta_uz(3,:))/dr;
% derivative_uz(jd,:) = (R(jd+1,:).*theta_uz(jd+1,:) - R(jd-1,:).*theta_uz(jd-1,:))./(2*R(jd,:)*dr);
% derivative_uz(nR,:) = (3*R(nR,:).*theta_uz(nR,:) - 4*R(nR-1,:).*theta_uz(nR-1,:) + R(nR-2,:).*theta_uz(nR-2,:))./(2*R(nR,:)*dr);
% 
% theta_uzeta = repmat(h, 1, nZ).^3.*repmat(theta, 1, nZ);
% derivative_uzeta = nan(nR, nZ); % pre-allocate (rz^2/2(z/3-h)Theta)_r/r
% derivative_uzeta(1,:) = (-3*theta_uzeta(1,:) + 4*theta_uzeta(2,:) - theta_uzeta(3,:))/dr;
% derivative_uzeta(jd,:) = (R(jd+1,:).*theta_uzeta(jd+1,:) - R(jd-1,:).*theta_uzeta(jd-1,:))./(2*R(jd,:)*dr);
% derivative_uzeta(nR,:) = (3*theta_uzeta(nR,:) - 4*theta_uzeta(nR-1,:) + theta_uzeta(nR-2,:))./(2*R(nR,:)*dr);
% 
% surf(R, Z, theta_uz,'EdgeColor','none','LineStyle','none','FaceLighting','phong'); colorbar; figure
% surf(R, Z, derivative_uz,'EdgeColor','none','LineStyle','none','FaceLighting','phong'); colorbar; figure

%--- Plot coefficeints
% Construct matrices
    Xi = R; dxi = dr;
    jXi = (2:nR-1)'; jZeta = (2:nZ-1)'; dzeta = Zeta(1,2) - Zeta(1,1);
    h_mat = repmat(h, 1, nZ);
    gb_mat = repmat(gb, 1, nZ);
    bar_phi_mat = repmat(bar_phi, 1, nZ);
    theta_mat = repmat(theta, 1, nZ);
    % Compute dUdZ
    dUdZ = nan(nR, nZ);
    dUdZ(:,1) = (-3*Uz(:,1) + 4*Uz(:,2) - Uz(:,3))/(2*dzeta);
    dUdZ(:,jZeta) = (Uz(:,jZeta+1) - Uz(:,jZeta-1))/(2*dzeta);
    dUdZ(:,nZ) = (3*Uz(:,nZ) - 4*Uz(:,nZ-1) + Uz(:,nZ-2))/(2*dzeta);
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
    a2 = Uz./h_mat - Zeta./h_mat.*dHdt + gamma*Zeta.^2.*(Zeta/2 - ones(nR, nZ)).*h_mat.*theta_mat.*dHdX;
    a3 = gb - Psi_d - dUdZ./h_mat - gamma*Zeta.*(Zeta - ones(nR, nZ)).*h_mat.*theta_mat.*dHdX;
    a1 = a1.*(h_mat >= threshold); a2 = a2.*(h_mat >= threshold); a3 = a3.*(h_mat >= threshold);

    surf(R, Z, a1,'EdgeColor','none','LineStyle','none','FaceLighting','phong'); colorbar; 
    xlabel('\(r\)','FontSize',20,'Interpreter','latex')
    ylabel('\(z\)','FontSize',20,'Interpreter','latex')
    zlabel('\(A_1\)','Interpreter','latex','FontSize',20); figure
    surf(R, Z, a2,'EdgeColor','none','LineStyle','none','FaceLighting','phong'); colorbar; 
    xlabel('\(r\)','FontSize',20,'Interpreter','latex')
    ylabel('\(z\)','FontSize',20,'Interpreter','latex')
    zlabel('\(A_2\)','Interpreter','latex','FontSize',20); figure
    surf(R, Z, a3,'EdgeColor','none','LineStyle','none','FaceLighting','phong'); colorbar; 
    xlabel('\(r\)','FontSize',20,'Interpreter','latex')
    ylabel('\(z\)','FontSize',20,'Interpreter','latex')
    zlabel('\(A_3\)','Interpreter','latex','FontSize',20); figure

%----------------------- Check boundary conditions ------------------------
dr = r(2) - r(1); dz = Zeta(1,2) - Zeta(1,1); z = Zeta(1,:);
ddR = (-3*Vol_Frac(1,:) + 4*Vol_Frac(2,:) - Vol_Frac(3,:))/(2*dr);
ddZ = (-3*Vol_Frac(:,1) + 4*Vol_Frac(:,2) - Vol_Frac(:,3))/(2*dz);
plot(z, ddR); xlabel('\(\zeta\)', 'Interpreter', 'latex', 'FontSize', 16); ylabel('\(\phi_\xi\)', 'Interpreter', 'latex', 'FontSize', 16); figure
plot(r, ddZ); xlabel('\(\xi\)', 'Interpreter', 'latex', 'FontSize', 16); ylabel('\(\phi_\zeta\)', 'Interpreter', 'latex', 'FontSize', 16);

%--------------------------------------------------------------------------
end