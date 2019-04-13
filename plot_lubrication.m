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
    phi = importdata(['vol_frac-',num2str(time_step),'.csv']); Vol_Frac = reshape(phi, nR, nZ);  
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
    xlim([0, R_dim]); ylim([0 2]);
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
nB = find(h < 1.1e-4, 1)-1;
h_r = [ (-3*h(1) + 4*h(2) - h(3))/(2*dr) ; (h(3:nPoints) - h(1:nPoints-2))/(2*dr) ; (3*h(nPoints) - 4*h(nPoints-1) + h(nPoints-2))/(2*dr) ];
h_rr = [ (-3*h_r(1) + 4*h_r(2) - h_r(3))/(2*dr) ; (h_r(3:nPoints) - h_r(1:nPoints-2))/(2*dr); (3*h_r(nPoints) - 4*h_r(nPoints-1) + h_r(nPoints-2))/(2*dr) ];
h_rrr = [ (-3*h_rr(1) + 4*h_rr(2) - h_rr(3))/(2*dr) ; (h_rr(3:nPoints) - h_rr(1:nPoints-2))/(2*dr); (3*h_rr(nPoints) - 4*h_rr(nPoints-1) + h_rr(nPoints-2))/(2*dr) ];
theta = (h_rrr + h_rr./r - h_r./(r.^2)); theta(1) = 0; % surface tension term
% h_r = [ (-3*h(1) + 4*h(2) - h(3))/(2*dr) ; (h(3:nB) - h(1:nB-2))/(2*dr) ; (3*h(nB) - 4*h(nB-1) + h(nB-2))/(2*dr) ];
% h_rr = [ (-3*h_r(1) + 4*h_r(2) - h_r(3))/(2*dr) ; (h_r(3:nB) - h_r(1:nB-2))/(2*dr) ; (3*h_r(nB) - 4*h_r(nB-1) + h_r(nB-2))/(2*dr) ];
% h_rrr = [ (-3*h_rr(1) + 4*h_rr(2) - h_rr(3))/(2*dr) ; (h_rr(3:nB) - h_rr(1:nB-2))/(2*dr) ; (3*h_rr(nB) - 4*h_rr(nB-1) + h_rr(nB-2))/(2*dr) ];
% theta = nan(nR, 1);
% theta(1:nB) = h_rrr + h_rr./r(1:nB) - h_r./(r(1:nB).^2);
% theta(nB+1:nR) = 0;
theta(1) = 0;
theta_r = [ (-3*theta(1) + 4*theta(2) - theta(3))/(2*dr) ; (theta(3:nR) - theta(1:nR-2))/(2*dr) ; (3*theta(nR) - 4*theta(nR-1) + theta(nR-2))/(2*dr) ];
set(gca, 'FontSize', 16) % change axis tick font size
plot(r, h.^2.*theta); xlabel('\(r\)', 'Interpreter', 'latex'); ylabel('\(h\Theta\)', 'Interpreter', 'latex'); figure
plot(r, theta_r); xlabel('\(r\)', 'Interpreter', 'latex'); ylabel('\(\Theta_r\)', 'Interpreter', 'latex'); figure

%----------------------- Check boundary conditions ------------------------
dr = r(2) - r(1); dz = Zeta(1,2) - Zeta(1,1); z = Zeta(1,:);
ddR = (-3*Vol_Frac(1,:) + 4*Vol_Frac(2,:) - Vol_Frac(3,:))/(2*dr);
ddZ = (-3*Vol_Frac(:,1) + 4*Vol_Frac(:,2) - Vol_Frac(:,3))/(2*dz);
plot(z, ddR); xlabel('\(\zeta\)', 'Interpreter', 'latex', 'FontSize', 16); ylabel('\(\phi_\xi\)', 'Interpreter', 'latex', 'FontSize', 16); figure
plot(r, ddZ); xlabel('\(\xi\)', 'Interpreter', 'latex', 'FontSize', 16); ylabel('\(\phi_\zeta\)', 'Interpreter', 'latex', 'FontSize', 16);

%--------------------------------------------------------------------------
end