function plot_lubrication_simplified
%--------------------------------------------------------------------------
%PLOT_LUBRICATION Plot solutions to lubrication model.
%   Alex Tam, 18/03/2019.
%--------------------------------------------------------------------------
%--------------------------- Import global data ---------------------------
r = importdata('r.csv'); R_dim = importdata('dish_size.csv');
t = importdata('t.csv'); contact_line = importdata('contact_line.csv');
thickness_index = importdata('thickness_index.csv');
times = 1:10000:100001;

%------------------------- Plot solution profiles -------------------------
for plots = 1:length(times)
    % Import solution data
    time_step = times(plots);
    h = importdata(['biofilm_height-',num2str(time_step),'.csv']);   
    gb = importdata(['biofilm_nutrient-',num2str(time_step),'.csv']);
    gs = importdata(['substratum_nutrient-',num2str(time_step),'.csv']);
    % Plot variable profiles
    hold on
    set(gca, 'FontSize', 16) % change axis tick font size
    plot(r, h, 'LineWidth', 1.5); 
    plot(r, gb, 'LineWidth', 1.5);
    plot(r, gs, 'LineWidth', 1.5);
    legend({'\(h\)', '\(g_b\)', '\(g_s\)'}, 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast')
    xlabel('\(r\)', 'Interpreter', 'latex')
    ylabel('Numerical solution', 'Interpreter', 'latex')
    title(['t = ', num2str(t(times(plots))),'.'], 'Interpreter', 'latex')
    xlim([0, R_dim]); ylim([0 2]);
    print(gcf, '-depsc', ['radial_ext_sol-',num2str(time_step),'.eps'])
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
nPoints = length(r); dr = r(2) - r(1); 
h_r = [ (-3*h(1) + 4*h(2) - h(3))/(2*dr) ; (h(3:nPoints) - h(1:nPoints-2))/(2*dr) ; (3*h(nPoints) - 4*h(nPoints-1) + h(nPoints-2))/(2*dr) ];
h_rr = [ (-3*h_r(1) + 4*h_r(2) - h_r(3))/(2*dr) ; (h_r(3:nPoints) - h_r(1:nPoints-2))/(2*dr); (3*h_r(nPoints) - 4*h_r(nPoints-1) + h_r(nPoints-2))/(2*dr) ];
h_rrr = [ (-3*h_rr(1) + 4*h_rr(2) - h_rr(3))/(2*dr) ; (h_rr(3:nPoints) - h_rr(1:nPoints-2))/(2*dr); (3*h_rr(nPoints) - 4*h_rr(nPoints-1) + h_rr(nPoints-2))/(2*dr) ];
theta = h.^3.*(h_rrr + h_rr./r - h_r./(r.^2)); theta(1) = 0; % surface tension term
theta_r = [ (-3*theta(1) + 4*theta(2) - theta(3))/(2*dr) ; (theta(3:nPoints) - theta(1:nPoints-2))/(2*dr) ; (3*theta(nPoints) - 4*theta(nPoints-1) + theta(nPoints-2))/(2*dr) ];
set(gca, 'FontSize', 16) % change axis tick font size
plot(r, theta); xlabel('\(r\)', 'Interpreter', 'latex'); ylabel('\(\Theta\)', 'Interpreter', 'latex'); figure
plot(r, theta_r); xlabel('\(r\)', 'Interpreter', 'latex'); ylabel('\(\Theta_r\)', 'Interpreter', 'latex')
%--------------------------------------------------------------------------
end