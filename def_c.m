function M = def_c(f_ck, h, w, o_L, n_L, d_L)

%--------------------------------------------------------------------------
% This fanction is defined be the user and converts the input to a
% structual array carrying all the information throughout the program.
%
% Input:    Characteristic strength and cross section dimensions in [mm]
%               - Height
%               - Width
%               - Diamter of longitudinal reinforcement rows
%               - Number of bars in each row
%               - Depth of middle rows (left blank if only top and bottom)
%
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
    M.L = 4.180;                                 %[m]
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
    M.f_yk = 550 *10^6;                         %[MPa]
    % Partial factor (Slap armering)
    M.gamma_s = 1.20 * M.gamma_3;               %[-]

% LOAD
    % Design value of axial load (ULS)
    M.N_Ed = 455.1 *10^3;                         %[kN]
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
    M.e_1 = 55 *10^(-3);                        %[mm]

    % The characterustuc strain at maximum force
    % (For Class B reinforcement steel, smalles value selected)
    M.varepsilon_uk = 5.0 /100;

    % Spacing between horizontal reinforcement
    % (This is assumes to be the following, assumed to around the mean OF)
    M.s_T = 175 *1E-3;   %[m]


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