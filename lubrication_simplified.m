function lubrication_simplified
%--------------------------------------------------------------------------
% LUBRICATION_SIMPLIFIED Solve the simplified lubrication model.
%   Solves the lubrication model for biofilm expansion, assuming constant
%   cell volume fraction.
%   Alex Tam, 13/03/2019.
%--------------------------------------------------------------------------
%-------------------------- Model assumptions -----------------------------
epsilon = 0.1; % [-] biofilm aspect ratio
H0 = 0.1; % [-] maximum initial biofilm height
phi = 0.9; % [-] cell volume fraction
precursor = 1e-4; % [-] precursor film height
threshold = 1.1*precursor; % [-] source term threshold
dlmwrite('threshold.csv', threshold);

%------------------------ Dimensional parameters --------------------------
R_b = 2.875; % [mm] Initial biofilm size
R_p = 41.5; % [mm] Petri dish radius
D_s = (1 - 0.023*0.3)*4.04e-2; % [mm^2/min] Diffusivity of glucose in substratum
D_b = 0.24*4.04e-2; % [mm^2/min] Diffusivity of glucose in biofilm
Q = 2.92e-3; % [mm/min] Mass transfer coefficient
G = 0.5/(pi*R_p^2); % [g/mm^2] Initial nutrient concentration (0.5g)

%-------------------------- Fitted parameters -----------------------------
psi_n = 12.1; % [mm^2/g/min] cell production rate
eta = 3.7e-3; % [/min] nutrient consumption rate

%----------------------- Dimensionless parameters -------------------------
R = R_p/R_b; % [-] dimensionless Petri dish radius
Upsilon = eta*R_b^2/D_b; % [-] nutrient consumption rate
Qs = Q*R_b/(epsilon*D_s); % [-] substratum permeability
Qb = Q*R_b/(epsilon*D_b); % [-] biofilm permeability
D = D_s/(psi_n*G*R_b^2); % [-] nutrient diffusion coefficient
Pe = (psi_n*G*R_b^2)/D_b; % [-] Peclet number
gamma = 1; % [-] surface tension coefficient
T = 14208*psi_n*G; % [-] dimensionless time in the experiment

%------------------------- Numerical parameters ---------------------------
nPoints = 1001; j = (3:nPoints-2)'; jd = (2:nPoints-1)';
nTimes = 20001;
r = linspace(0, R, nPoints)'; dr = r(2) - r(1);
t = linspace(0, T, nTimes); dt = t(2) - t(1);
output_files = 10; % number of files to output
dlmwrite('r.csv', r, 'precision', 20);
dlmwrite('dish_size.csv', R, 'precision', 20);
dlmwrite('t.csv', t, 'precision', 20);

%-------------------------- Initial conditons -----------------------------
H = precursor + (H0-precursor)*(ones(size(r)) - r.^2).*(r < 1); % initial condition
Gs = ones(nPoints, 1); % [-] initial substratum nutrient concentration within the biofilm
Gb = 0*ones(nPoints, 1); % [-] initial biofilm nutrient concentration

%------------------------------ Solve PDEs --------------------------------
for i = 1:nTimes-1
    % Extract data
    if sum(i == 1:(nTimes-1)/output_files:nTimes) == 1
        dlmwrite(['biofilm_height-',num2str(i),'.csv'], H, 'precision', '%.20f');
        dlmwrite(['biofilm_nutrient-',num2str(i),'.csv'], Gb, 'precision', '%.20f');
        dlmwrite(['substratum_nutrient-',num2str(i),'.csv'], Gs, 'precision', '%.20f');
    end
    % Store variables
    h = H; gs = Gs; gb = Gb; J = source(gb, h, threshold);
    % 1. Biofilm height
    Ah = matrix_h(j, r, nPoints, dt, dr, gamma, h); bh = rhs_h(j, r, nPoints, precursor, J, dt, dr, gamma, h);
    H = Ah\bh;
    % 2. Nutrient concentration (substratum)
    Ags = matrix_gs(jd, r, dr, dt, D, Qs, nPoints); bgs = rhs_gs(jd, r, dr, dt, D, Qs, nPoints, gs, gb);
    Gs = Ags\bgs;
    % Compute surface tension term: expanded form, sequential differentiation (2nd order)
    h_r = [ (-3*h(1) + 4*h(2) - h(3))/(2*dr) ; (h(3:nPoints) - h(1:nPoints-2))/(2*dr) ; (3*h(nPoints) - 4*h(nPoints-1) + h(nPoints-2))/(2*dr) ];
    h_rr = [ (-3*h_r(1) + 4*h_r(2) - h_r(3))/(2*dr) ; (h_r(3:nPoints) - h_r(1:nPoints-2))/(2*dr); (3*h_r(nPoints) - 4*h_r(nPoints-1) + h_r(nPoints-2))/(2*dr) ];
    h_rrr = [ (-3*h_rr(1) + 4*h_rr(2) - h_rr(3))/(2*dr) ; (h_rr(3:nPoints) - h_rr(1:nPoints-2))/(2*dr); (3*h_rr(nPoints) - 4*h_rr(nPoints-1) + h_rr(nPoints-2))/(2*dr) ];
    theta = h.^3.*(h_rrr + h_rr./r - h_r./(r.^2));
    theta(1) = 0; % enforced by velocity condition
    % 3. Nutrient concentration (biofilm)
    Agb = matrix_gb(h, r, dr, nPoints, gamma, Pe, Qb, Upsilon, phi, jd, dt, theta); bgb = rhs_gb(h, r, dr, nPoints, gamma, Pe, Qb, Upsilon, phi, jd, dt, theta, gb, gs);
    Gb = Agb\bgb;
end

dlmwrite(['biofilm_height-',num2str(i+1),'.csv'], H, 'precision', '%.20f');
dlmwrite(['biofilm_nutrient-',num2str(i+1),'.csv'], Gb, 'precision', '%.20f');
dlmwrite(['substratum_nutrient-',num2str(i+1),'.csv'], Gs, 'precision', '%.20f');
        
%% Functions
%----------------------------- Source term --------------------------------
function J = source(Gb, h, threshold)
    J = Gb.*(h >= threshold); % modified source term
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
function A = matrix_gb(h, r, dr, nPoints, gamma, Pe, Qb, Upsilon, phi, jd, dt, theta)
    % Construct linear system
    index1 = [ 1 ; 1 ; jd ; jd ; jd ; nPoints ; nPoints ];
    index2 = [ 1 ; 2 ; jd+1 ; jd ; jd-1 ; nPoints-1 ; nPoints ];
    entries = [ Pe*h(1) + dt*gamma*(1-phi)*Pe*(-3*theta(1) + 4*theta(2) - theta(3))/(6*dr) + 2*dt*h(1)/(dr^2) + dt*Qb/2 + dt*Upsilon*phi*h(1)/2 ; -2*dt*h(1)/(dr^2) ; ... % ghost points for no-flux condition
        dt*gamma*(1-phi)*Pe*(r(jd+1).*theta(jd+1))./(12*dr*r(jd)) - dt*(r(jd) + r(jd+1)).*(h(jd) + h(jd+1))./(8*r(jd)*dr^2) ; Pe*h(jd) + dt*(r(jd) + r(jd+1)).*(h(jd) + h(jd+1))./(8*r(jd)*dr^2) + dt*(r(jd) + r(jd-1)).*(h(jd) + h(jd-1))./(8*r(jd)*dr^2) + dt*Qb/2 + dt*Upsilon*phi*h(jd)/2 ; -dt*gamma*(1-phi)*Pe*(r(jd-1).*theta(jd-1))./(12*dr*r(jd)) - dt*(r(jd) + r(jd-1)).*(h(jd) + h(jd-1))./(8*r(jd)*dr^2) ; ... % interior grid points
        -dt*h(nPoints)/(dr^2) ; Pe*h(nPoints) + dt*gamma*(1-phi)*Pe*(3*theta(nPoints) - 4*theta(nPoints-1) + theta(nPoints-2))/(12*r(nPoints)*dr) + dt*h(nPoints)/dr^2 + dt*Qb/2 + dt*Upsilon*phi*h(nPoints)/2 ]; % ghost points for no-flux condition
    A = sparse(index1, index2, entries, nPoints, nPoints);
end

function b = rhs_gb(h, r, dr, nPoints, gamma, Pe, Qb, Upsilon, phi, jd, dt, theta, gb, gs)
    b = nan(nPoints,1);
    b(1) = gb(1)*(Pe*h(1) - dt*gamma*(1-phi)*Pe*(-3*theta(1) + 4*theta(2) - theta(3))/(6*dr) - 2*dt*h(1)/(dr^2) - dt*Qb/2 - dt*Upsilon*phi*h(1)/2)...
        + gb(2)*(2*dt*h(1)/(dr^2)) + dt*Qb*gs(1); % ghost points for no-flux condition
    b(jd) = gb(jd).*(Pe*h(jd) - dt*(r(jd) + r(jd+1)).*(h(jd) + h(jd+1))./(8*r(jd)*dr^2) - dt*(r(jd) + r(jd-1)).*(h(jd) + h(jd-1))./(8*r(jd)*dr^2) - dt*Qb/2 - dt*Upsilon*phi*h(jd)/2) ...
        + gb(jd+1).*(-dt*gamma*(1-phi)*Pe*(r(jd+1).*theta(jd+1))./(12*dr*r(jd)) + dt*(r(jd) + r(jd+1)).*(h(jd) + h(jd+1))./(8*r(jd)*dr^2)) ...
        + gb(jd-1).*( dt*gamma*(1-phi)*Pe*(r(jd-1).*theta(jd-1))./(12*dr*r(jd)) + dt*(r(jd) + r(jd-1)).*(h(jd) + h(jd-1))./(8*r(jd)*dr^2))...
        + dt*Qb*gs(jd); % interior grid points
    b(nPoints) = gb(nPoints)*(Pe*h(nPoints) - dt*gamma*(1-phi)*Pe*(3*theta(nPoints) - 4*theta(nPoints-1) + theta(nPoints-2))/(12*r(nPoints)*dr) - dt*h(nPoints)/dr^2 - dt*Qb/2 - dt*Upsilon*phi*h(nPoints)/2)...
        + gb(nPoints-1)*(dt*h(nPoints)/(dr^2)) + dt*Qb*gs(nPoints); % ghost points for no-flux condition
end
%--------------------------------------------------------------------------
end