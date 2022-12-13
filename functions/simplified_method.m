function M_Ed = simplified_method(info_lev, M)

%--------------------------------------------------------------------------
% Calculated the design value of the bending moment via the simplified
% method of nominal curvature. 
%
% This method estimates the maximum deformation of the column, be assuming
% breaking strain, varepsilon_cu, in the concrete and yield strain in the
% reinforcement steel, which is the worst case scenario.
%
% The method is desribed in the book:
% Beton konstruktion efter DS/EN 1992-1-1, section 7.3.2, page 29
%
% Input:    Information about the member (structure array)
% Output:   Approximated design value of the bending moment
%--------------------------------------------------------------------------



% The maximum deformation
% (1/10 is familiar to the c and c_0 parameter seen in the previous two
% methods.)
u_max = 1/10 * (M.varepsilon_cu + M.varepsilon_yd)/M.d * M.L_0^2;   %[m]

% The desgin value of the bending moment
M_Ed = M.M_0Ed + M.N_Ed*u_max;                                      %[Nm]


%--------------------------------------------------------------------------
%% Outputting requested information


if info_lev == 1
    fprintf('u_max = %.2f [mm] | M_Ed_sm = %.2f [kNm] | ' ...
        , u_max*10^(3), M_Ed*10^(-3))

elseif info_lev == 2
    fprintf(['The largest possible deformation is: ' ...
        'u_max = %.4f [mm] \n' ...
        'The design value of the bending moment is: ' ...
        'M_Ed_sm = %.4f [kNm] \n'] ...
        ,u_max*10^(3) ,M_Ed*10^(-3))

end