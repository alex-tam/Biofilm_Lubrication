function lubrication_simplified_sensitivity
%--------------------------------------------------------------------------
% LUBRICATION_SIMPLIFIED Solve the simplified lubrication model.
%   Solves the lubrication model for biofilm expansion, assuming constant
%   cell volume fraction.
%   Alex Tam, 20/05/2019.
%--------------------------------------------------------------------------
%-------------------------- Model assumptions -----------------------------
epsilon = 0.1; % [-] biofilm aspect ratio
H0 = 0.1; % [-] maximum initial biofilm height
phi = 0.9; % [-] cell volume fraction
precursor = 1e-4; % [-] precursor film height
h_star = 0.002; % [-] source term threshold

%------------------------ Dimensional parameters --------------------------
R_b = 2.875; % [mm] Initial biofilm size
R_p = 41.5; % [mm] Petri dish radius
D_s = (1 - 0.023*0.3)*4.04e-2; % [mm^2/min] Diffusivity of glucose in substratum
D_b = 0.24*4.04e-2; % [mm^2/min] Diffusivity of glucose in biofilm
Q = 2.92e-3; % [mm/min] Mass transfer coefficient
G = 0.5/(pi*R_p^2); % [g/mm^2] Initial nutrient concentration (0.5g)
psi_n = 50; % [mm^2/g/min] cell production rate
eta = 3.7e-3; % [/min] nutrient consumption rate

%----------------------- Dimensionless parameters -------------------------
R = R_p/R_b; % [-] dimensionless Petri dish radius
Upsilon = eta*R_b^2/D_b; % [-] nutrient consumption rate
Qs = Q*R_b/(epsilon*D_s); % [-] substratum permeability
Qb = Q*R_b/(epsilon*D_b); % [-] biofilm permeability
D = D_s/(psi_n*G*R_b^2); % [-] nutrient diffusion coefficient
Pe = (psi_n*G*R_b^2)/D_b; % [-] Peclet number
gamma = 1; % [-] surface tension coefficient
T = 7200*psi_n*G; % [-] dimensionless time in the experiment

%------------------------- Numerical parameters ---------------------------
nR = 1001; j = (3:nR-2)'; jd = (2:nR-1)';
nTimes = 100001;
r = linspace(0, R, nR)'; dr = r(2) - r(1);
t = linspace(0, T, nTimes); dt = t(2) - t(1);
output_files = 10; % number of files to output
dlmwrite('r.csv', r, 'precision', 20);
dlmwrite('dish_size.csv', R, 'precision', 20);
dlmwrite('t_constant_vf.csv', t, 'precision', 20);

%---------------------------- Vary parameters -----------------------------
VAR = [ 1 ; 5 ]';
dlmwrite('var.csv', VAR);

% Pre-allocate test-dependent quantities
final_size = nan(size(VAR));
ridge = nan(size(VAR));
thickness = nan(size(VAR));

for tests = 1:length(VAR)
    var = VAR(tests);
    Qs = var;

%-------------------------- Initial conditons -----------------------------
    H = precursor + (H0-precursor)*(ones(size(r)) - (r/1).^2).^4.*(r <= 1); % initial condition
    Gs = ones(nR, 1); % [-] initial substratum nutrient concentration within the biofilm
    Gb = 0*ones(nR, 1); % [-] initial biofilm nutrient concentration
    contact_line = r(find(H <= h_star, 1))*ones(nTimes,1); % [-] initial contact line position
    thickness_index = max(H)/contact_line(1)*ones(nTimes, 1); % [-] initial thickness index
    ridge_index = max(H)/H(1)*ones(nTimes, 1); % [-] initial ridge index

%------------------------------ Solve PDEs --------------------------------
    for i = 1:nTimes-1
        % Store variables
        h = H; gs = Gs; gb = Gb; J = source(gb, h, h_star);
        % Write data
        if sum(i == 1:(nTimes-1)/output_files:nTimes) == 1
            dlmwrite(['biofilm_height-',num2str(i),'-var-',num2str(var),'.csv'], H, 'precision', '%.20f');
            dlmwrite(['biofilm_nutrient-',num2str(i),'-var-',num2str(var),'.csv'], Gb, 'precision', '%.20f');
            dlmwrite(['substratum_nutrient-',num2str(i),'-var-',num2str(var),'.csv'], Gs, 'precision', '%.20f');
            dlmwrite(['source-',num2str(i),'-var-',num2str(var),'.csv'], J, 'precision', '%.20f');
        end
        % 1. Biofilm height
        Ah = matrix_h(j, r, nR, dt, dr, gamma, h); bh = rhs_h(j, r, nR, precursor, J, dt, dr, gamma, h);
        H = Ah\bh;
        % 2. Nutrient concentration (substratum)
        Ags = matrix_gs(jd, r, dr, dt, D, Qs, nR, h, h_star); bgs = rhs_gs(jd, r, dr, dt, D, Qs, nR, gs, gb, h, h_star);
        Gs = Ags\bgs;
        % Compute surface tension term
        theta = h.^3.*surface_tension(h, r, dr, nR, h_star);
        % 3. Nutrient concentration (biofilm)
        Agb = matrix_gb(h, r, dr, nR, gamma, Pe, Qb, Upsilon, phi, jd, dt, theta, h_star); bgb = rhs_gb(h, r, dr, nR, gamma, Pe, Qb, Upsilon, phi, jd, dt, theta, gb, gs, h_star);
        Gb = Agb\bgb;
        % Store contact line position
        contact_line(i+1) = r(find(H <= h_star, 1)); % [-] contact line position
        thickness_index(i+1) = max(H)/contact_line(i+1); % [-] thickness index
        ridge_index(i+1) = max(H)/H(1); % [-] ridge index
    end

    dlmwrite(['biofilm_height-',num2str(i+1),'-var-',num2str(var),'.csv'], H, 'precision', '%.20f');
    dlmwrite(['biofilm_nutrient-',num2str(i+1),'-var-',num2str(var),'.csv'], Gb, 'precision', '%.20f');
    dlmwrite(['substratum_nutrient-',num2str(i+1),'-var-',num2str(var),'.csv'], Gs, 'precision', '%.20f');
    dlmwrite(['source-',num2str(i+1),'-var-',num2str(var),'.csv'], source(Gb, H, h_star), 'precision', '%.20f');
    dlmwrite(['contact_line','-var-',num2str(var),'.csv'], contact_line, 'precision', '%.5f');
    dlmwrite(['It_lubrication_constant_vf','-var-',num2str(var),'.csv'], thickness_index, 'precision', '%.5f');
    dlmwrite(['Ir_lubrication_constant_vf','-var-',num2str(var),'.csv'], ridge_index, 'precision', '%.5f');

    % Store test-dependent quantities
    final_size(tests) = contact_line(i+1);
    ridge(tests) = ridge_index(i+1);
    thickness(tests) = thickness_index(i+1);
end

dlmwrite('final_size.csv', final_size)
dlmwrite('ridge.csv', ridge);
dlmwrite('thickness.csv', thickness);

%% Functions
%--------------------------- Surface tension ------------------------------
function theta = surface_tension(h, r, dr, nR, h_star)
    theta = nan(nR,1);
    theta(1) = 0;
    theta(2) = (r(3)*(h(4) - h(2)) - r(2)*(h(3) - h(1)))/(dr^3*(r(3)+r(2))) - r(2)*(h(3) - h(1))/(r(2)*dr^3);
    theta(3:nR-2) = ((r(5:nR) + r(4:nR-1)).*(h(5:nR) - h(4:nR-1)) - (r(4:nR-1) + r(3:nR-2)).*(h(4:nR-1) - h(3:nR-2)))./(4*r(4:nR-1)*dr^3) ...
        - ((r(3:nR-2) + r(2:nR-3)).*(h(3:nR-2) - h(2:nR-3)) - (r(2:nR-3) + r(1:nR-4)).*(h(2:nR-3) - h(1:nR-4)))./(4*r(2:nR-3)*dr^3);
    theta(nR-1) = (r(nR)*(3*h(nR) - 4*h(nR-1) + h(nR-2)) - r(nR-1)*(h(nR) - h(nR-2)))/(dr^3*(r(nR) + r(nR-1)))...
        + (r(nR-2)*(h(nR-1) - h(nR-3)) - r(nR-1)*(h(nR) - h(nR-2)))/(dr^3*(r(nR-1) + r(nR-2)));
    theta(nR) = 3/r(nR)*(9*r(nR)*h(nR) - 12*r(nR)*h(nR-1) + 3*r(nR)*h(nR-2) - 4*(r(nR) + r(nR-1))*(h(nR) - h(nR-1)) + r(nR-1)*(h(nR) - h(nR-2)))/(2*dr^3) ...
        - 4/(r(nR) + r(nR-1))*(3*r(nR)*h(nR) - 4*r(nR)*h(nR-1) + r(nR)*h(nR-2) - r(nR-1)*(h(nR) - h(nR-2)))/(dr^3) ...
        + 1/(r(nR-1))*((r(nR) + r(nR-1))*(h(nR) - h(nR-1)) - (r(nR-1) + r(nR-2))*(h(nR-1) - h(nR-2)))/(2*dr^3);
    theta = theta.*(h >= h_star);
end

%----------------------------- Source term --------------------------------
function J = source(Gb, h, threshold)
    J = Gb.*h.*(h >= threshold); % modified source term
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
       gamma*dt/6*cp.*rpp./dp ; gamma*dt/6*(-cp.*rpp./dp - cp.*rp./dp - cp.*rp./d - cn.*rp./d) ; 1 + gamma*dt/6*(cp.*rp./dp + cp.*rp./d + cp.*rn./d + cn.*rp./d + cn.*rn./d + cn.*rn./dn) ; gamma*dt/6*(-cp.*rn./d - cn.*rn./d - cn.*rn./dn -cn.*rnn./dn) ; gamma*dt/6*cn.*rnn./dn ; ... % interior grid points
       gamma*dt/6*(Cn*Rnn/Dn) ; gamma*dt/6*(-Cp*Rn/Dm - Cn*Rn/Dm - Cn*Rn/Dn - Cn*Rnn/Dn) ; 1 + gamma*dt/6*(Cp*Rpp/Dp + Cp*Rp/Dp + Cp*Rp/Dm + Cp*Rn/Dm + Cn*Rp/Dm + Cn*Rn/Dm + Cn*Rn/Dn) ; gamma*dt/6*(-Cp*Rpp/Dp - Cp*Rp/Dp - Cp*Rp/Dm - Cn*Rp/Dm) ; ... % ghost points for no-flux condition
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

%---------------- Linear system for substratum nutrient -------------------
function A = matrix_gs(jd, r, dr, dt, D, Qs, nR, h, h_star)
    index1 = [ 1 ; 1 ; jd ; jd ; jd ; nR ; nR ];
    index2 = [ 1 ; 2 ; jd+1 ; jd ; jd-1 ; nR-1 ; nR ];
    entries = [ 1 + 2*dt*D/dr^2 + dt*D*Qs/2*(h(1) >= h_star) ; -2*dt*D/dr^2 ; ... % ghost points for no-flux condition
       -dt*D/(4*dr^2)*(r(jd+1) + r(jd))./r(jd) ; ones(nR-2, 1) + dt*D/(dr^2) + dt*D*Qs/2*(h(jd) >= h_star) ; -dt*D/(4*dr^2)*(r(jd) + r(jd-1))./r(jd) ; ... % interior grid points
       -dt*D/(dr^2) ; 1 + dt*D/(dr^2) + dt*D*Qs/2*(h(nR) >= h_star) ]; % ghost points for no-flux condition
    A = sparse(index1, index2, entries, nR, nR);
end

function b = rhs_gs(jd, r, dr, dt, D, Qs, nR, gs, gb, h, h_star)
    b = nan(nR,1);
    b(1) = gs(1)*(1 - 2*dt*D/dr^2 - dt*D*Qs/2*(h(1) >= h_star)) + gs(2)*(2*dt*D/dr^2) + dt*D*Qs*gb(1)*(h(1) >= h_star); % ghost points for no-flux condition
    b(jd) = gs(jd).*(1 - dt*D/(dr^2) - dt*D*Qs/2*(h(jd) >= h_star))...
        + gs(jd+1).*(dt*D/(4*dr^2)*(r(jd+1) + r(jd))./r(jd))...
        + gs(jd-1).*(dt*D/(4*dr^2)*(r(jd) + r(jd-1))./r(jd)) + dt*D*Qs*gb(jd).*(h(jd) >= h_star); % interior grid points
    b(nR) = gs(nR)*(1 - dt*D/(dr^2) - dt*D*Qs/2*(h(nR) >= h_star)) + gs(nR-1)*(dt*D/(dr^2)) + dt*D*Qs*gb(nR)*(h(nR) >= h_star); % ghost points for no-flux condition
end

%------------------ Linear system for biofilm nutrient --------------------
function A = matrix_gb(h, r, dr, nR, gamma, Pe, Qb, Upsilon, phi, jd, dt, theta, h_star)
    index1 = [ 1 ; 1 ; 1 ; jd ; jd ; jd ; nR ; nR ];
    index2 = [ 1 ; 2 ; 3 ; jd+1 ; jd ; jd-1 ; nR-1 ; nR ];
    entries = [ -3/(2*dr) ; 4/(2*dr) ; -1/(2*dr) ; ... % Dirichlet condition g_r(0) = 0
        dt*Pe*gamma*(1-phi)*(r(jd+1).*h(jd+1).^3.*theta(jd+1))./(12*r(jd)*dr) - (h(jd) >= h_star).*(dt*(r(jd) + r(jd+1)).*(h(jd) + h(jd+1))./(8*r(jd)*dr^2)) ; ...
        Pe*h(jd) + (h(jd) >= h_star).*(dt*(r(jd) + r(jd+1)).*(h(jd) + h(jd+1))./(8*r(jd)*dr^2) + dt*(r(jd) + r(jd-1)).*(h(jd) + h(jd-1))./(8*r(jd)*dr^2) + dt*Qb/2 + dt*Upsilon*phi*h(jd)/2) ; ...
        -dt*Pe*gamma*(1-phi)*(r(jd-1).*h(jd-1).^3.*theta(jd-1))./(12*r(jd)*dr) - (h(jd) >= h_star).*(dt*(r(jd) + r(jd-1)).*(h(jd) + h(jd-1))./(8*r(jd)*dr^2)) ; ... % interior grid points
        -dt*h(nR)/(dr^2)*(h(nR) >= h_star) ; Pe*h(nR) + dt*gamma*Pe*(1-phi)*(3*h(nR)^3*theta(nR) - 4*h(nR-1)^3*theta(nR-1) + h(nR-2)^3*theta(nR-2))/(12*r(nR)*dr) + (h(nR) >= h_star)*(dt*h(nR)/dr^2 + dt*Qb/2 + dt*Upsilon*phi*h(nR)/2) ]; % ghost point for no-flux condition
    A = sparse(index1, index2, entries, nR, nR);
end

function b = rhs_gb(h, r, dr, nR, gamma, Pe, Qb, Upsilon, phi, jd, dt, theta, gb, gs, h_star)
    b = nan(nR,1);
    b(1) = 0; % Dirichlet condition g_r(0) = 0
    b(jd) = gb(jd).*( Pe*h(jd) - (h(jd) >= h_star).*(dt*(r(jd) + r(jd+1)).*(h(jd) + h(jd+1))./(8*r(jd)*dr^2) + dt*(r(jd) + r(jd-1)).*(h(jd) + h(jd-1))./(8*r(jd)*dr^2) + dt*Qb/2 + dt*Upsilon*phi*h(jd)/2) )...
          + gb(jd+1).*( -dt*Pe*gamma*(1-phi)*(r(jd+1).*h(jd+1).^3.*theta(jd+1))./(12*r(jd)*dr) + (h(jd) >= h_star).*(dt*(r(jd) + r(jd+1)).*(h(jd) + h(jd+1))./(8*r(jd)*dr^2)) )...
          + gb(jd-1).*( dt*Pe*gamma*(1-phi)*(r(jd-1).*h(jd-1).^3.*theta(jd-1))./(12*r(jd)*dr) + (h(jd) >= h_star).*(dt*(r(jd) + r(jd-1)).*(h(jd) + h(jd-1))./(8*r(jd)*dr^2)) )...
          + dt*Qb*gs(jd).*(h(jd) >= h_star); % interior grid points
    b(nR) = gb(nR)*( Pe*h(nR) - dt*gamma*Pe*(1-phi)*(3*h(nR)^3*theta(nR) - 4*h(nR-1)^3*theta(nR-1) + h(nR-2)^3*theta(nR-2))/(12*r(nR)*dr) - (h(nR) >= h_star)*(dt*h(nR)/dr^2 + dt*Qb/2 + dt*Upsilon*phi*h(nR)/2) )...
          + gb(nR-1)*( dt*h(nR)/(dr^2)*(h(nR) >= h_star) )...
          + dt*Qb*gs(nR)*(h(nR) >= h_star); % ghost point for no-flux condition
end
end