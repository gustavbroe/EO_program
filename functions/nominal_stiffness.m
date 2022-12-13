function [M_Ed, nr] = nominal_stiffness(info_lev, cut_off, M)

%--------------------------------------------------------------------------
% Calculated the design value of the bending moment via the method of
% nominal stiffness.
%
% The procedure is varies depanding on the geometric reinforcement ratio.
%
% {Nominal stiffness - 5.8.7.2 (1) - p.70}
% The following technique can be used to estimate the nominal stiffness of
% slender compression memebers with arbitrary cross section.
%
% Input:    Information about the member (structure array)
% Output:   Design value of the bending moment
%--------------------------------------------------------------------------

% Assuming that the cross section is NOT non-reinforced
nr = 0;

% High geometric reinforcement ratio
if M.rho >= 0.01

    K_s = 0;                                %[-]            %{(5.26)}
    K_c = 0.3 / (1 + 0.5*M.varphi_ef);      %[-]


% Medium (normal) geometric reinforcement ratio 
% rho >= 0.002
else

    % Adds output if the cross section is considered non-reinforced
    if M.rho < 0.002

        % Changing assumption
        nr = 1;

    end

    % Factor accounting for the conrete strength class
    k_1 = sqrt( (M.f_ck*10^(-6))/20 );      %[-]            %{(5.23)}

    % Factor accounting for the axial force and slenderness
    k_2 = min(M.n * M.lambda/170, 0.20);    %[-]            %{(5.24)}

    % K_s determines if the reinforcement steel has influence over the result
    K_s = 1;                                %[-]            %{(5.22)}

    % K_c accounts for the influence of cracking, creep, etc.
    K_c = (k_1*k_2) / (1 + M.varphi_ef);    %[-]            %{(5.22)}

end

% Calculating the nominal stiffness
EI = K_c*M.E_cd*M.I_c + K_s*M.E_s*M.I_s;    %[Nm^2]         %{(5.21)}

% Buckling load based on nominal stiffness
N_B = (pi^2*EI)/M.L_0^2;                    %[N]            %{BK (7.7)}

% Factoring distrubution of 1st and 2nd moment
    % (Triangular moment distribution, simply supported, moment in one end)
    if M.r_m == 1
        c_0 = 9;                            %[-]            %{5.8.7.3 (2)}
    else
        c_0 = 8;
    end
    % (sine-wave shaped distribution of the second order bending moment)
    beta_ns = (pi^2) / c_0;                 %[-]            %{(5.29)}


% Checking if the buckling load is greater than the applied load
% (With a reduction of the applied load)
if N_B > cut_off*M.N_Ed

    % Total design moment (including second ordre effects)
    if M.r_m == 1
        M_Ed = M.M_0Ed * ( 1 + beta_ns/( (N_B/M.N_Ed)-1 ) );%{5.8.7.3 (1)}
    else
        M_Ed = M.M_0e * ( 1 + beta_ns/( (N_B/M.N_Ed)-1 ) ); %{5.8.7.3 (1)}
    end
else

    % The function outputs NaN
    M_Ed = NaN;
end


% % The initial modulus of elasticity (idealised stress strain curve)
% E_c0m = 1.05*E_cm;                    %[Pa]               %{7.9 - p.206}
% 
% % Buckling tenson for concrete columns (Ritter's equation) [Pa]
% sigma_cr = f_cd / (1 + f_ck/(pi^2 * E_c0m) * (L_0/i)^2 ); %{7.12 - p.206}
% 
% % Critical tenson for reinforcement steel
% sigma_sc = E_s*sigma_cr/E_cm;         %[Pa]               %{7.16 - p.207}
% 
% % Buckling load for reinforced concrete
% N_B = sigma_cr*A_c + sigma_sc*A_s;    %[N]                %{7.15 - p.207}



%--------------------------------------------------------------------------
%% Outputting requested information

if info_lev >= 1
    if nr == 1
        % Printing error message
        fprintf(['NON-REINFORCED | s = %.1f [mm] | o_L = %.1f [mm] | ' ...
        'rho = %.4f <= 0.002 [-] | '] ...
        , M.h*10^3, M.o_L(end)*10^3, M.rho);
    end
end


if info_lev == 1
    fprintf(['(EI, N_B, M_Ed_ns) = (%.2f, %.2f, %.2f) ' ...
        '[kNm^2, kN, kNm] | '], EI*10^(-3), N_B*10^(-3), M_Ed*10^(-3));

elseif info_lev == 2
    fprintf(['(I_c, I_s, K_c, K_s, EI, N_B, M_Ed_ns) = (%.2e, %.2e, ' ...
        '%.2f, %.2f, %.2f, %.2f, %.2f) ' ...
        '[mm^4, mm^4, -, -, kNm^2, kN, kNm] | '], M.I_c*10^12, M.I_s*10^12, K_c, K_s, ...
        EI*10^(-3), N_B*10^(-3), M_Ed*10^(-3));

elseif info_lev == 3
    fprintf(['The nominal stifness of the column is: ' ...
        'EI = %.2f [kNm^2] \n' ...
        'The buckling load with nominal stiffness is: ' ...
        'N_B = %.2f [kN] \n' ...
        'The design value of the second order bending moment is: ' ...
        'M_Rd_ns = %.2f [kNm] \n'], EI*10^(-3), N_B*10^(-3), M_Ed*10^(-3))

end