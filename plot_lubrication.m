function plot_lubrication
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
times = 1:1000:10001;

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
% Integral_uz = nan(nR,nZ);
%     for row = 1:nR
%         Integral_uz(row, :) = cumtrapz(Z(row,:), Vol_Frac(row,:));
%     end
%     theta_uz = Z.^2/2.*(Z/3 - repmat(h, 1, nZ)).*repmat(theta, 1, nZ);
%     derivative_uz = nan(nR, nZ);
%     derivative_uz(1,:) = (-3*theta_uz(1,:) + 4*theta_uz(2,:) - theta_uz(3,:))/dr;
%     derivative_uz(jd,:) = (R(jd+1,:).*theta_uz(jd+1,:) - R(jd-1,:).*theta_uz(jd-1,:))./(2*R(jd,:)*dr);
%     derivative_uz(nR,:) = (3*theta_uz(nR,:) - 4*theta_uz(nR-1,:) + theta_uz(nR-2,:))./(2*R(nR,:)*dr);
%     surf(R, Z, theta_uz,'EdgeColor','none','LineStyle','none','FaceLighting','phong'); colorbar; figure
%     surf(R, Z, derivative_uz,'EdgeColor','none','LineStyle','none','FaceLighting','phong'); colorbar; figure
%----------------------- Check boundary conditions ------------------------
dr = r(2) - r(1); dz = Zeta(1,2) - Zeta(1,1); z = Zeta(1,:);
ddR = (-3*Vol_Frac(1,:) + 4*Vol_Frac(2,:) - Vol_Frac(3,:))/(2*dr);
ddZ = (-3*Vol_Frac(:,1) + 4*Vol_Frac(:,2) - Vol_Frac(:,3))/(2*dz);
plot(z, ddR); xlabel('\(\zeta\)', 'Interpreter', 'latex', 'FontSize', 16); ylabel('\(\phi_\xi\)', 'Interpreter', 'latex', 'FontSize', 16); figure
plot(r, ddZ); xlabel('\(\xi\)', 'Interpreter', 'latex', 'FontSize', 16); ylabel('\(\phi_\zeta\)', 'Interpreter', 'latex', 'FontSize', 16);

%--------------------------------------------------------------------------
end