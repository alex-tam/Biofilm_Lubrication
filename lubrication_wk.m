function lubrication_wk
%--------------------------------------------------------------------------
% LUBRICATION_WK Solve the Ward & King (2012) no-slip model.
%   Solves the lubrication model for biofilm expansion, assuming constant
%   cell volume fraction and no nutrient depletion.
%   Alex Tam, 16/05/2019.
%--------------------------------------------------------------------------
%---------------------------- Model parameters ----------------------------
gamma = 1; % [-] surface tension coefficient
sigma = 24.995; rho = 0.2; % [-] source term parameters
precursor = 1e-4; % [-] precursor flim height
h_star = 0.02; % [-] source term threshold
H0 = 1; % [-] initial biofilm height

%-------------------------- Numerical parameters --------------------------
nPoints = 1001; nTimes = 13001; % numerical grid sizes
j = (3:nPoints-2)'; % indices of interior grid points
r = linspace(0, 10, nPoints)'; dr = r(2) - r(1); % spatial domain
t = linspace(0, 130, nTimes); dt = t(2) - t(1); % time domain
output_files = 13; % number of files to output
dlmwrite('r.csv', r, 'precision', '%.5f');

%--------------------------- Initial condition ----------------------------
h = precursor + (H0-precursor)*(ones(size(r)) - (r/1.5).^2).^4.*(r <= 1.5); % Ward & King

%------------------------------- Solve PDE --------------------------------
for i = 1:nTimes-1
    % Store old variables
    h_old = h; J = source(h_old, sigma, rho, h_star);
    % Write data
    if sum(i == 1:(nTimes-1)/output_files:nTimes) == 1
        dlmwrite(['biofilm_height-',num2str(i),'.csv'], h, 'precision', '%.5f');
        dlmwrite(['source-',num2str(i),'.csv'], J, 'precision', '%.5f');
    end
    % Construct and solve linear system
    A = matrix_h(j, r, nPoints, dt, dr, gamma, h_old); b = rhs_h(j, r, nPoints, precursor, J, dt, dr, gamma, h_old);
    h = A\b;
end

dlmwrite(['biofilm_height-',num2str(i+1),'.csv'], h, 'precision', '%.5f');
dlmwrite(['source-',num2str(i+1),'.csv'], J, 'precision', '%.5f');

%------------------------------ Source term -------------------------------
function J = source(h, sigma, rho, h_star)
    J = (1/(sqrt(sigma))*tanh(sqrt(sigma)*h) - rho*h).*(h >= h_star); % [-] regularised source term
end

%----------------------------- Linear system ------------------------------
function A = matrix_h(j, r, nPoints, dt, dr, gamma, h)
    % common coefficients for interior points
    cp = ((r(j+1) + r(j)).*(h(j+1) + h(j)).^3)./(16*r(j)*dr^2); cn = ((r(j) + r(j-1)).*(h(j) + h(j-1)).^3)./(16*r(j)*dr^2);
    rpp = r(j+2) + r(j+1); rp = r(j+1) + r(j); rn = r(j) + r(j-1); rnn = r(j-1) + r(j-2);
    dp = 2*r(j+1)*dr^2; d = 2*r(j)*dr^2; dn = 2*r(j-1)*dr^2;
    % coefficients for N-1 boundary
    Cp = ((r(nPoints) + r(nPoints-1))*(h(nPoints) + h(nPoints-1))^3)/(16*r(nPoints-1)*dr^2); Cn = ((r(nPoints-1) + r(nPoints-2))*(h(nPoints-1) + h(nPoints-2))^3)/(16*r(nPoints-1)*dr^2);
    Rpp = 2*r(nPoints) + dr; Rp = r(nPoints) + r(nPoints-1); Rn = r(nPoints-1) + r(nPoints-2); Rnn = r(nPoints-2) + r(nPoints-3);
    Dp = 2*r(nPoints)*dr^2; D = 2*r(nPoints-1)*dr^2; Dn = 2*r(nPoints-2)*dr^2;
    % construct matrix for linear system
    index_1 = [ 1 ; 1 ; 1 ; 2 ; 2 ; 2 ; 2 ; 2 ; j ; j ; j ; j ; j ; nPoints-1 ; nPoints-1 ; nPoints-1 ; nPoints-1 ; nPoints ];
    index_2 = [ 1 ; 2 ; 3 ; 1 ; 2 ; 3 ; 4 ; 5 ; j+2 ; j+1 ; j ; j-1 ; j-2 ; nPoints ; nPoints-1 ; nPoints-2 ; nPoints-3 ; nPoints ];
    entries = [ -3/(2*dr) ; 4/(2*dr) ; -1/(2*dr) ; ... % no-flux condition
        -5/(2*dr^3) ; 18/(2*dr^3) ; -24/(2*dr^3) ; 14/(2*dr^3) ; -3/(2*dr^3) ; ... % third derivative condition
        gamma*dt/6*cp.*rpp./dp ; gamma*dt/6*(-cp.*rpp./dp - cp.*rp./dp - cp.*rp./d - cn.*rp./d) ; 1 + gamma*dt/6*(cp.*rp./dp + cp.*rp./d + cp.*rn./d + cn.*rp./d + cn.*rn./d + cn.*rn./dn) ; gamma*dt/6*(-cp.*rn./d - cn.*rn./d - cn.*rn./dn -cn.*rnn./dn) ; gamma*dt/6*cn.*rnn./dn ; ... % interior grid points
        gamma*dt/6*(-Cp*Rpp/Dp - Cp*Rp/Dp - Cp*Rp/D - Cn*Rp/D) ; 1 + gamma*dt/6*(Cp*Rpp/Dp + Cp*Rp/Dp + Cp*Rp/D + Cp*Rn/D + Cn*Rp/D + Cn*Rn/D + Cn*Rn/Dn) ; gamma*dt/6*(-Cp*Rn/D -Cn*Rn/D - Cn*Rn/Dn - Cn*Rnn/Dn) ; gamma*dt/6*(Cn*Rnn/Dn) ; ... % ghost point method for no-flux condition
        1 ]; % Dirichlet condition at r = R
    A = sparse(index_1, index_2, entries, nPoints, nPoints);
%     A = sparse(1,1,-3/(2*dr),nPoints,nPoints)...
%       + sparse(1,2,4/(2*dr),nPoints,nPoints)...
%       + sparse(1,3,-1/(2*dr),nPoints,nPoints)... % no-flux condition
%       + sparse(2,1,-5/(2*dr^3),nPoints,nPoints)...
%       + sparse(2,2,18/(2*dr^3),nPoints,nPoints)...
%       + sparse(2,3,-24/(2*dr^3),nPoints,nPoints)...
%       + sparse(2,4,14/(2*dr^3),nPoints,nPoints)...
%       + sparse(2,5,-3/(2*dr^3),nPoints,nPoints)... % third derivative condition
%       + sparse(j, j+2, gamma*dt/6*cp.*rpp./dp, nPoints, nPoints)...
%       + sparse(j, j+1, gamma*dt/6*(-cp.*rpp./dp - cp.*rp./dp - cp.*rp./d - cn.*rp./d), nPoints, nPoints)...
%       + sparse(j, j, 1 + gamma*dt/6*(cp.*rp./dp + cp.*rp./d + cp.*rn./d + cn.*rp./d + cn.*rn./d + cn.*rn./dn), nPoints, nPoints)...
%       + sparse(j, j-1, gamma*dt/6*(-cp.*rn./d - cn.*rn./d - cn.*rn./dn -cn.*rnn./dn), nPoints, nPoints)...
%       + sparse(j, j-2, gamma*dt/6*cn.*rnn./dn, nPoints, nPoints)... % interior grid points
%       + sparse(nPoints-1, nPoints, gamma*dt/6*(-Cp*Rpp/Dp - Cp*Rp/Dp - Cp*Rp/D - Cn*Rp/D), nPoints, nPoints)...
%       + sparse(nPoints-1, nPoints-1, 1 + gamma*dt/6*(Cp*Rpp/Dp + Cp*Rp/Dp + Cp*Rp/D + Cp*Rn/D + Cn*Rp/D + Cn*Rn/D + Cn*Rn/Dn), nPoints, nPoints)...
%       + sparse(nPoints-1, nPoints-2, gamma*dt/6*(-Cp*Rn/D -Cn*Rn/D - Cn*Rn/Dn - Cn*Rnn/Dn), nPoints, nPoints)...
%       + sparse(nPoints-1, nPoints-3, gamma*dt/6*(Cn*Rnn/Dn), nPoints, nPoints)... % ghost points for no-flux
%       + sparse(nPoints, nPoints, 1, nPoints, nPoints); % Dirichlet condition
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
    Dp = 2*r(nPoints)*dr^2; D = 2*r(nPoints-1)*dr^2; Dn = 2*r(nPoints-2)*dr^2;
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
                 + (-gamma*dt/6*(-Cp*Rpp/Dp - Cp*Rp/Dp - Cp*Rp/D - Cn*Rp/D))*h(nPoints)...
                 + (-gamma*dt/6*(Cp*Rpp/Dp + Cp*Rp/Dp + Cp*Rp/D + Cp*Rn/D + Cn*Rp/D + Cn*Rn/D + Cn*Rn/Dn))*h(nPoints-1)...
                 + (-gamma*dt/6*(-Cp*Rn/D - Cn*Rn/D - Cn*Rn/Dn - Cn*Rnn/Dn))*h(nPoints-2)...
                 + (-gamma*dt/6*(Cn*Rnn/Dn))*h(nPoints-3); % ghost points for no-flux
    b(nPoints) = precursor; % Dirichlet condition
end
%--------------------------------------------------------------------------
end