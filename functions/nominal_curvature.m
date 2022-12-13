function M_Ed = nominal_curvature(info_lev, M)

%--------------------------------------------------------------------------
% Calculated the design value of the bending moment via the method of
% nominal curvature.
%
% {Curvataure - 5.8.8.3 (1) - p.73}
% The method is currently valid for symmetrical cross section. The 
% effective depth, d, is different if the reinforcement is not
% concentrated on opposite sides, see {5.8.8.3 (2)}.
%
% {Method based on nominal curvature  - 5.8.8.1 (1) - p.72}
% The method of nominal stiffness is well suited for isolated members
% with a constant axial force and a defined effective length, L_0.
%
% Input:    Information about the member (structure array)
% Output:   Design value of the bending moment
%--------------------------------------------------------------------------


% Factor depending on the curvature distribution
% (c_0 is solely dependent on the first order effects. c is dependent on
% the overall curvature. pi^2 (approx. 10) is usally used, which
% corresponds to a sinusoidal curvature distrubution, where 8 is the lower
% limit which is a constant curvature.)
c = pi^2;

% The relative axial force at maximum moment resistance
% (this parameter can be assumed to be the following value)
n_bal = 0.4;                                        %[-]    %{5.8.8.3 (3)}  

% Correction factor depending on axial load
n_u = 1 + M.omega;
K_r = min((n_u - M.n) / (n_u - n_bal), 1);          %[-]    %{5.8.8.3 (3)}

% Factor for taking account of creep (f_ck in [MPa])
beta = 0.35 + (M.f_ck*10^(-6))/200 - M.lambda/150;  %[-]    %{5.8.8.3 (4)}  
K_varphi = max(1 + beta*M.varphi_ef, 1);            %[-]    %{5.8.8.3 (4)}

% Curvature for members with constant symmetrical cross section
r_0 = 1 / (M.varepsilon_yd / (0.45*M.d));           %[m]            
r = 1 / (K_r*K_varphi * 1/r_0);                     %[m]    %{5.8.8.3 (1)}

% The deflection
e_2 = (1/r) * M.L_0^2 / c;                          %[m]    %{5.8.8.2 (3)}

% The nominal second order moment
M_2 = M.N_Ed*e_2;                                   %[Nm]   %{5.8.8.2 (3)}

% The total design bending moment
if M.r_m == 1
    M_Ed = M.M_0Ed + M_2;                           %[Nm]   %{5.8.8.2 (1)}
else
    M_Ed = M.M_0e + M_2;                            %[Nm]   %{5.8.8.2 (1)}
end


%--------------------------------------------------------------------------
%% Outputting requested information


if info_lev == 1
    fprintf(['(e_2, M_2, M_Ed_nc) = (%.2f, %.2f, %.2f) ' ...
        '[mm, kNm, kNm] | '], e_2*10^(3), M_2*10^(-3), M_Ed*10^(-3))

elseif info_lev == 2
    fprintf(['The deflection of the column is: e_2 = %.2f [mm] \n' ...
        'The nominal second order moment is: M_2 = %.2f [kNm] \n' ...
        'The total second order moment is: M_ed_nc = %.2f [kNm] \n'] ...
        , e_2*10^(3), M_2*10^(-3), M_Ed*10^(-3))

end