function M = bending_moment(info_lev, M)

%--------------------------------------------------------------------------
% This function uses information about load eccentricity, axial load, and 
% geometric imperfections to calculate the bending moment.
%
% Assumptions
%  - The moment distribution is parabolic.
%
% Input:    Information about the member (structure array)
% Output:   New information to the struct
%--------------------------------------------------------------------------

% CONTRIBUTION FROM TRANSVERSE LOAD
    % Maximum bending moment from line load
    % (standard solution)
    M_tv = 1/8 * M.p_tv * M.L^2;                    %[Nm]

% FIRST ORDER BENDING MOMENT
    % In design load combination (ULS)
    M.M_0Ed = M.N_Ed * (M.e_0 + M.e_1) + M_tv;      %[Nm]

    % In quasi-permanent load combination (SLS) 
    % (The load combination factor in quasi permanent, psi_2 combinations 
    % for industry and warehouse is 0.7, this is applied the service load, 
    % which approximatly accounts for 25% for both cases)
    M.M_0Eqp = M.N_Ed*(1 - 0.3*0.25) * (M.e_0 + M.e_1);           %[Nm]

% MOMENT RATIO BETWEEN THE ENDS
    % The moment through the members decrease linearly. If it is simply
    % supported and subject to a point moment.
    % (OBS: |M_02| >= |M_01| ) 
    M.M_01 = 0;                                     %[kNm]
    M.M_02 = M.M_0Ed;                               %[kNm]
    M.r_m = M.M_01/M.M_02;                          %[-]

% DIFFERING END MOMENT
    % If the end moment are differnt a equivalent can be used with the
    % assumption of a constant moment distribution. 
    M.M_0e = max(0.6*M.M_02 + 0.4*M.M_01, 0.4*M.M_02);


%--------------------------------------------------------------------------
%% Outputting requested information

if info_lev == 1
    fprintf('M_0Ed = %.2f [kNm] |', M.M_0Ed*10^(-3))

elseif info_lev == 2
    fprintf(['M_0Ed = %.2f [kNm] | M_0Eqp = %.2f [kNm] | ' ...
        'r_m =  %.1f [-] | '], M.M_0Ed*10^(-3), M.M_0Eqp*10^(-3) ...
        , M.r_m*10^(-3))

elseif info_lev == 3
    fprintf(['First order bending moment in design load combination ' ...
        '(ULS): M_0Ed = %.2f [kNm] \n'], M.M_0Ed*10^(-3));

    fprintf(['First order bending moment in quasi-permanent load ' ...
        'combination (SLS): M_0Eqp = %.2f [kNm] \n'], M.M_0Eqp*10^(-3));

    fprintf(['The ratio between the end moments (ULS): ' ...
        'M_01/M_02 = %.2f/%.2f = %.1f [kNm] \n'] ...
        , M.M_01*10^(-3), M.M_02*10^(-3), M.M_01*10^(-3)/(M.M_02*10^(-3)));

end