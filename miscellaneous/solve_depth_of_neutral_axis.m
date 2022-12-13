clear; clc; close all;

%--------------------------------------------------------------------------
% This function is not finished, can be developed in future work.
%--------------------------------------------------------------------------

n = 10000;
h = 240;

M = def_c(35, h, 360, [16, 16, 16], [5, 2, 5], [145]);
P = strength_parameters(M);
x = linspace(-2*h*1E-3, 2*h*1E-3, n);

P.A_sr = M.n_L .* pi .* (M.o_L/2).^2;
P.index_sc = 1:1;
P.index_s = 2:3;

df = zeros(n, 1);
for i = 1:n

    df(i) = soe(M, P, x(i)).lhs;

end

plot(x*1E3, df, 'b.')
ylim([-1E6, 1E6])
grid on


function R = soe(M, P, x)

% Strains in the reinforcement rows

    % Compression                                               %[-]
    if P.index_sc ~= 0
        R.varepsilon_sc = P.varepsilon_cu3 .* (x - M.d_L(P.index_sc)) / x;
    else
        R.varepsilon_sc = 0;
    end

    % Tension                                                   %[-]
    if P.index_s ~= 0
        R.varepsilon_s = P.varepsilon_cu3 .* (M.d_L(P.index_s) - x) / x;
    else 
        R.varepsilon_s = 0;
    end


% Stresses in rows
% (The strains can not be greater than the yeild strength of steel)

    % Compression
    R.sigma_sc = min(R.varepsilon_sc * M.E_s, M.f_yd);          %[Pa]

    % Tension
    R.sigma_s = min(R.varepsilon_s * M.E_s, M.f_yd);            %[Pa]


% Resulting force from the cross section components

    % Square stress distribution from concrete
    R.F_c = P.lambda*x*M.w*M.f_cd;                              %[N]

    % Compression steel
    R.F_sc = P.A_sr(P.index_sc) .* R.sigma_sc;                  %[N]

    % Tension steel
    R.F_s = P.A_sr(P.index_s) .* R.sigma_s;                     %[N]

% Left hand side of state of equilibrium equation
R.lhs = R.F_c + sum(R.F_sc) - sum(R.F_s) - M.N_Ed;              %[N]

end