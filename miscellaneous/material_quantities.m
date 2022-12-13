%--------------------------------------------------------------------------
% Introduction to script:
%
%   The script calculates and output the volume and mass of materials 
%   needed to construct both walls and columns, with dimensions just above 
%   the GWP-optimal.

%--------------------------------------------------------------------------
% Initializing MATLAB
clear; clc; close all hidden;

% Adding function folder to path
addpath('functions')

% Loading material information
MI = material_info(0);


%--------------------------------------------------------------------------
% Estimating the total quantities for walls

    % Total amount of bearing interior wall
    %  (Emilie, 2022, p. 8)
    A_w = (183 + 16 + 1137);            %[m^2]

    % Characteristic compressive strength of concrete
    % (16, 20, and 35 is currently available)
    f_ck = 16;   EPD_c = 'c16';   EPD_s = 'rp';

    % Wall thickness and spacing between reinforcement bars
    w = 3600;                           %[mm]
    s = 150;                            %[mm]

    % Definitions
    M = def_w(f_ck, 100, w, [10 10], [w/s w/s], []);

    % Steel (both vertical and horizontal)

        % Volumen
        V.s = (M.A_s*M.L) + M.V_sT;     %[m^3]

        % Mass
        mass.s = V.s * MI.(EPD_s).rho;  %[kg]

    % Concrete

        % Volumen
        V.c = (M.A_c*M.L) - V.s;        %[m^3]

        % Mass
        mass.c = V.c * MI.(EPD_c).rho;  %[kg]

    % Converting quantities to (x)/m^2
    q_s = mass.s/(M.L*w*10^(-3));       %[kg/m^2]
    q_c = V.c/(M.L*w*10^(-3));          %[m^3/m^2]

    % Printing result
    fprintf(['Material quantities needed to construct all bearing ' ...
        'interior wall are: (%.1f [m^2])\n\n'], A_w);

    fprintf(' - %.1f [m^3] of concrete (%.1f [m^3/m^2] / %.1f [ton])\n' ...
        ,q_c*A_w , q_c, q_c*A_w * MI.(EPD_c).rho *10^(-3));

    fprintf(' - %.1f [ton] of steel (%.1f [kg/m^2] / %.1f [m^3])\n\n' ...
        ,q_s*A_w *10^(-3) , q_s, q_s*A_w/MI.(EPD_s).rho);


%--------------------------------------------------------------------------
% Estimating the total quantities for columns

    % Total number of columns
    % (Yasemin, 2022, p. 8)
    n_c = 151;                          %[-]

    % Characteristic compressive strength of concrete
    % (16, 20, and 35 is currently available)
    f_ck = 16;   EPD_c = 'c16';   EPD_s = 'rp';

    % Cross sectional side length and diameter of longitunal reinforcement
    s = 220;                            %[mm]
    o_L = 20;                                     %[mm]

    % Definitions
    M = def_w(f_ck, s, s, [o_L o_L], [2 2], []);

    % Steel (both vertical and horizontal)

        % Volumen
        V.s = (M.A_s*M.L) + M.V_sT;     %[m^3]

        % Mass
        mass.s = V.s * MI.(EPD_s).rho;  %[kg]

    % Concrete

        % Volumen
        V.c = (M.A_c*M.L) - V.s;        %[m^3]

        % Mass
        mass.c = V.c * MI.(EPD_c).rho;  %[kg]

    % Printing result
    fprintf(['Material quantities needed to construct all columns ' ...
        'are: (%.f) \n\n'], n_c);

    fprintf(' - %.1f [m^3] of concrete (%.1f [m^3/column] / %.1f [ton])\n' ...
        ,V.c*n_c , V.c, mass.c*n_c *10^(-3));

    fprintf(' - %.1f [ton] of steel (%.1f [kg/m^2] / %.1f [m^3])\n\n' ...
        ,mass.s*n_c *10^(-3) , mass.s, V.s*n_c);