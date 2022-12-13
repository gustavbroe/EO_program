function e_0 = geo_imperfections(info_lev, M)

%--------------------------------------------------------------------------
% The geo_imperfections calculates the geometric imperfections, which is
% calucalted as an excentricity. These imperfections are an consequnce of
% casting, production, and construction of the concrete members.
%
% {Geometric imperfections - 5.2 (1) - p.54}
% If the column is affected by unfavoirable effetcs of deviations in
% geometry and the position of loads, these could be considered.
%
% {5.2 (2-3)}
% Geometric imperfections is to be considered in ULS (both accidental
% and persistent design situations) and SLS (Serviceability Limit States)
%
% {5.2 (4)}
% If following rules are valid for members with axial compression and
% structures with vertical load.
%
% Input:    Information about the member (structure array)
% Output:   An excentricity, e_0, in meters.
%--------------------------------------------------------------------------


% ECCENTRICITY CAUSED BY GEOMETRIC IMPERFECTION             %{5.2 - p.55}
    % Basic value of inclination (recommended value)
    theta_0 = 1/200;                            %[-]        %{5.2 (5)}

    % Reduction factor for length and height
    % (The length, L, is meters)
    alpha_h = max(2/3, min(2/sqrt(M.L) ,1));    %[-]        %{5.2 (5)}

        % Number of elements
        % (Isolated member)
        m = 1;                                  %[-]        %{5.2 (6)}

    % Reduction factor for number of elements
    alpha_m = sqrt(0.5*(1+1/m));                %[-]        %{5.2 (5)}

    % Represented as an inclination
    theta_i = theta_0*alpha_h*alpha_m;          %[-]        %{5.2 (5)}

    
% The final eccentricity from geometric imperfection
% (Isolated member)
% (NA states that e_0 has to be geater than L_0/200 {5.2(1)P - p.29})
e_0 = min(theta_i*M.L_0/2 ,M.L_0*1/200);          %[m]        %{5.2 (7.a)}



%--------------------------------------------------------------------------
%% Outputting requested information

if info_lev == 1
    fprintf('e_0 = %.2f [mm] | ', e_0*10^(3));

elseif info_lev == 2
    fprintf(['The geometric imperfections represented as an ' ...
        'eccentricity is: e_0 = %.2f [mm] \n'], e_0*10^(3));

end