function [total, concrete, steel] = calc_GWP(info_lev, EPD_types, dimensions)

%--------------------------------------------------------------------------
% Calculates the Global Warming Potential (GWP) based on climate
% information from the 'climate-info'-function.
%
%
% Input:    EPD_types (String array)
%           Information about the member (structure array) 
% Output:   GWP for column (total) and the concrete and steel components 
%           respectively. [kg CO2-ep]
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Translates function inputs

% EPD types
EPD_concrete = EPD_types{1};
EPD_steel = EPD_types{2};

% Dimensions
M = dimensions;


%--------------------------------------------------------------------------
% Converting moment capacity to Global Warming Potential (GWP)

% Imports climate information
MI = material_info(0);


% The total volume of the column for each combination (square cross section)
V.tot = M.A_c * M.L;        %[m^3]

% Volumes of ...

    % Steel (both vertical and horizontal)
    V.s = M.A_s*M.L + M.V_sT;      %[m^3]

    % Concrete
    V.c = V.tot - V.s;      %[m^3]

% Weight of ...

    % Steel
    mass.s = V.s * MI.(EPD_steel).rho;

    % Concrete
    mass.c = V.c * MI.(EPD_concrete).rho;

% GWP of ...

    % Steel
    steel = V.s*MI.(EPD_steel).GWP;  %[kg CO2-eq]
    
    % Concrete
    concrete = V.c*MI.(EPD_concrete).GWP; %[kg CO2-eq]

% Total GWP
total = steel + concrete;


%--------------------------------------------------------------------------
% Printing requested information

if info_lev == 1

    % Only the total result
    fprintf('GWP = %.1f [kg CO2-eq.] | ', total)

elseif info_lev == 2

    % The total and individual results
    fprintf('GWP_(tot, s, c) = (%.1f, %.1f, %.1f) [kg CO2-eq.] | ' ...
        ,total, steel, concrete);

elseif info_lev == 3

    % Supplementary information and the results
    fprintf(['The total global warming potential from cradle to cradle' ...
        ' for the member is %.2f [kg CO2-eq.]. \n'], total)

    fprintf(['Steel accounts for %.1f%% (%.1f [kg CO2-eq.]) ' ...
        'of the total GWP, with a volumen of %.1e [m^3], ' ...
        'a mass of %.1f [kg]. (%.2f [kg CO2-eq. pr. m^3] | %s) \n']...
        , steel/total*100, steel, V.s, mass.s, ...
        MI.(EPD_steel).GWP, EPD_steel)

    fprintf(['Concrete accounts for %.1f%% (%.1f [kg CO2-eq.]) ' ...
        'of the total GWP, with a volumen of %.1e [m^3], ' ...
        'a mass of %.1f [kg]. (%.2f [kg CO2-eq. pr. m^3] | %s) \n']...
        , concrete/total*100, concrete, V.c, mass.c, ...
        MI.(EPD_concrete).GWP, EPD_concrete)

end


% Ending the function
end