%% Initializing MATLAB
clear; clc; close all hidden;

% Adding function folder to path
addpath('functions')

% Loading climate information and plot colors
MI = material_info(0);   PC = plot_colors();

%-------------------------------------------------------------------------
% Number of points
np = 100;

% Applied axial load
N_Ed = linspace(0, 4000, np);

% Characteristic compressive strength of concrete
% (16, 20, 25, 30 [F/R]  |  35, 40 [L/R])
f_ck = 35;   EPD_c = 'c35R';   EPD_s = 'rp';

% Number of longitunal reinforcement rows
n_L = [2, 2];

% Reinforcement and cross section dimensions
h = 300;   w = 300;   o_L = [16, 16];

%--------------------------------------------------------------------------
%% Calculating applied moment, moment capacity and GWP for each combination

% Rreallocating arrays and starting counter
ss = zeros(np, 1);   oo=ss;   MM=ss;   rr=ss;   M0 = ss;
SM = ss;   NC = ss;   NS = ss;   nn = ss;   GG = ss;   PP = ss;

% Loop for side length
for x = 1:size(N_Ed, 2)

    % Definitions
    M = def_N(N_Ed(x), f_ck, h, w, o_L, n_L, []);

    % Geometrical imperfection
    M.e_0 = geo_imperfections(0, M);
    
    % First order bending moment
    M = bending_moment(0, M);

    % Effects of creep
    M.varphi_ef = effective_creep(0, M);

    % Moment capacity
    [MM(x), rr(x)] = moment_capacity(0, M);

    % Applied second order moment
        
        %Simplified method II
        SM(x) = simplified_method(0, M);

        % Nominal curvature
        NC(x) = nominal_curvature(0, M);

        % Nominal stiffness
        [NS(x), nn(x)] = nominal_stiffness(0, 1.01, M);

    % Global warming potential
    GG(x) = calc_GWP(0, {EPD_c, EPD_s}, M);

    % Price
    PP(x) = calc_Price(0, M, EPD_c, EPD_s);

    % First order moment
    M0(x) = M.M_0Ed;

end

%--------------------------------------------------------------------------
% Post processing

% Line thickness
LT = 2;

% Creating new figure for the new plot
p1 = figure('Name', 'Force', 'PaperSize', [21 19]);

% Making it the current figure and insuring that future enteries will be
% added to it.
figure(p1); hold on;

% Plotting the results (moment as a function of applied force)
    % Moment capacity
    plot(N_Ed, MM*1E-3, 'k-', 'LineWidth', LT, ...
        'DisplayName', 'Moment capacity, M_R_d');

    % First order moment
    plot(N_Ed, M0*1E-3, 'k--', 'LineWidth', LT, ...
        'DisplayName', 'First order moment, M_0_E_d');

    % The three estimation methods
    plot(N_Ed, NS*1E-3, 'r-', 'LineWidth', LT, ...
        'DisplayName', 'Nominal stiffness, M_E_d_n_s');
    plot(N_Ed, NC*1E-3, 'b-', 'LineWidth', LT, ...
        'DisplayName', 'Nominal curvature, M_E_d_n_c');
    plot(N_Ed, SM*1E-3, 'b--', 'LineWidth', LT, ...
        'DisplayName', 'Simplified method, M_E_d_s_m');


% Setting axis limits
ylim([0, max(MM)*1E-3*1.1]);

% Adding plot title
title('', ...
    sprintf(['(f_c_k = %.f [MPa] / h = w = %.f [mm] / L_0 = %.2f [m]' ...
    ' / c+o_T = %.f [mm])'] ...
    , M.f_ck*10^(-6), h, M.L_0, (M.c+M.o_T)*10^(3) ));

% Changing plot legends
xlabel('Applied axial load,{\it N_E_d} [N]')
ylabel('Bending moment,{\it M} [Nm]')

% Concluding plot
conclude_plot(p1, 'right')

%--------------------------------------------------------------------------
% Concluding plot

function conclude_plot(plot_identifier, screen_side)

% Adding major and minor grid lines
grid on; 

    % Finding the 'Axes'-objects
    ax = findall(plot_identifier, 'Type', 'Axes');

    % Adding minor grid lines for both axis
    set(ax, 'XMinorGrid', 'on');
    set(ax, 'YMinorGrid', 'on');


% Gathering display information
    % Current screen size
    si = get(0,'screensize');
    % Height of toolbar
    th = si(4)*0.06;

% The height/width-ratio of the figure PaperSize
hwr = plot_identifier.PaperSize(2)/plot_identifier.PaperSize(1);

% Changing the position of the plot
if strcmp(screen_side, 'left')

    % To the left side
    set(plot_identifier, 'Position', ...
        [0 th si(3)*0.5 si(3)*0.5*hwr])

elseif strcmp(screen_side, 'right')

    % To the right side
    set(plot_identifier, 'Position', ...
        [si(3)*0.5 th si(3)*0.5 si(3)*0.5*hwr])

end


% Inserting legends
legend;

% Concluding creation of the plot
hold off

end

%--------------------------------------------------------------------------
% Nested definition function

function M = def_N(N_Ed, f_ck, h, w, o_L, n_L, d_L)

%--------------------------------------------------------------------------
% This fanction is defined be the user and converts the input to a
% structual array carrying all the information.
%
% Input:    Height and width of the cross section.
% Output:   A structure array with geometrical and mechanical properties.
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Input parameters [SI-Units]

% CROSS SECTION
    % Height and width
    M.h = h *10^(-3);                           %[mm]
    % Width
    M.w = w *10^(-3);                           %[mm]
    % Area
    M.A_c = M.h * M.w;                          %[m^2]
    % Moment of inertia (rectangular)
    M.I_c = 1/12*M.w*M.h^3;                     %[m^4]

% REINFORCEMENT
    % Thickness of concret cover
    M.c = 20 *10^(-3);                          %[mm]
    % Diameter (longitudinal)
    M.o_L = o_L *10^(-3);                           %[mm]
    % Number of bars in each row [no less than 2]
    M.n_L = n_L;                              %[-]
    % Diameter (transversal)
    M.o_T = 6 *10^(-3);                         %[mm]
    % Depth to first row of longitudinal bars
    if size(M.n_L, 2) > 1
        M.d_sc = M.c + M.o_T + M.o_L(1)/2;      %[m]
    else
        M.d_sc = [];
    end
    % Effective depth
    M.d = M.h - (M.c + M.o_T + M.o_L(end)/2);    %[m]
    % Depth of the bars center of each row
    M.d_L = [M.d_sc, (d_L*10^(-3)), M.d];      %[m]

% COLUMN
    % Length
    M.L = 4.18;                                 %[m]
    %Effective length factor
    M.L_0 = 1*M.L;                              %[m]

% PARTIAL FACTORS
    % Control class (Normal)
    M.gamma_3 = 1.00;                           %[-]

% CONCRETE
    % Characteristic compressive strength
    M.f_ck = f_ck *10^6;                          %[MPa]
    % Secant modulus at 40% of the average compressive strength
    M.E_cm = 33 *10^9;                          %[GPa]
    % Ultimate strain
    M.varepsilon_cu = 0.35 /100;                %[%]
    % Partial factor (Reinforced and under compression)
    M.gamma_c = 1.45 * M.gamma_3;               %[-]

% STEEL
    % Modolus of elesticity
    M.E_s = 200 *10^9;                          %[GPa]
    % Yeild strength
    M.f_yk = 500 *10^6;                         %[MPa]
    % Partial factor (Slap armering)
    M.gamma_s = 1.20 * M.gamma_3;               %[-]

% LOAD
    % Design value of axial load (ULS)
    M.N_Ed = N_Ed *10^3;                         %[kN]
    % Transverse line load (from wind, etc.)
    % (NOTE: A transversal load would change the shape of the moment 
    % distribution, meaning that c_0 regarding nominal stiffness method 
    % has to be changed accordingly.) 
    M.p_tv = 0 *10^3;                          %[kN/m]


%--------------------------------------------------------------------------    
% ASSUMPTIONS
    % Concrete age at loading
    % (In-situ concrete has to remain in the cast for at less 3 days.
    % Lisbet Sand, DTU Campusservice, said during the tour of B112, that
    % the each level was constructed over 1-2 weeks.
    % In-situ concrete has to remain in the cast for at less 3 days,
    % then the rest of building is built in the following 2-4 month.
    % 30 days is assumed, as load is progressively added.)
    M.t0 = 30;                         %[days]

    % Eccentricity of load
    % (A good assumption, advised by PG)
    M.e_1 = 0 *10^(-3);                        %[mm]

    % The characterustuc strain at maximum force
    % (For Class B reinforcement steel, smalles value selected)
    M.varepsilon_uk = 5.0 /100;

    % Spacing between horizontal reinforcement
    % (This is assumes to be minimum value stated by EC2)
    M.s_T = min([20*min(M.o_L), min(M.h, M.w), 400*10^(-3)]);   %[m]


%--------------------------------------------------------------------------
% Preliminary calculations

% AREAS
    % Total area of vertical steel
    M.A_s = 0;
    for i = 1:size(M.o_L, 2)
        M.A_s = M.A_s + M.n_L(i) * pi * (M.o_L(i)/2)^2;         %[m^2]
    end
    
    % Second moment of area (about center)
    M.I_s = 0;
    for i = 1:size(M.o_L, 2)
        M.I_s = M.I_s + M.n_L(i) * (pi/64*M.o_L(i)^4 + ...
            pi * (M.o_L(i)/2)^2 * (M.h/2 - M.d_L(i))^2);     %[m^4]
    end

% VOLUMES
    % Total volumen of horizontal reinforcement steel
    % (A loop is assumed to consistent of four rods, pairs of two with 
    % equal length)                     %[m^3]
    M.V_sT = M.L/M.s_T * pi * (M.o_T/2)^2 * ...
        ( 2*(M.h - (2*M.c+M.o_T)) + 2*(M.w - (2*M.c+M.o_T)) );  %[m^3]                                               

% CONVERTING CHARACTERISTIC VALUES TO DESIGN VALUES
    % Yeild strength of steel
    M.f_yd = M.f_yk / M.gamma_s;                %[Pa]
    % Compressive strength of concrete
    M.f_cd = M.f_ck / M.gamma_c;                %[Pa]
    % Modulus of concrete
    M.E_cd = M.E_cm / M.gamma_c;                %[Pa]
    % Yield strain of steel
    M.varepsilon_yd = M.f_yd/M.E_s;             %[-]
    % strain at maximum force of steel
    M.varepsilon_ud = M.varepsilon_uk/M.gamma_s;%[-]

% REINFORCEMENT RATIO
    % Geometric reinforcement ratio
    M.rho = M.A_s / M.A_c;                      %[-]
    % Mechanical reinforcement ratio
    M.omega = (M.A_s*M.f_yd) / (M.A_c*M.f_cd);  %[-]

% ADDITIONAL GEOMETRIC PARAMETERS
    % Radius of gyration of the uncracked concrete section
    M.i_c = sqrt(M.I_c / M.A_c);                %[m]
    % Slenderness ratio
    M.lambda = M.L_0 / M.i_c;                   %[-]

% ADDITIONAL GEOMETRIC PARAMETERS
    % Calculating the relative axial force
    M.n = M.N_Ed / (M.A_c*M.f_cd);              %[-]        %{5.8.7.2 (2)}

end