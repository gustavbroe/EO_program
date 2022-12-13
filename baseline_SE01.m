%% Initializing MATLAB
clear; clc; close all;

% Adding function folder to path
addpath('functions', 'use')

% Loading materials information
MI = material_info(1);

% Input control?
IC = false;


% ------------------- Inputs equal to column type SE01 -------------------
% (Information not stated below is defined in the def_c.m function.)

% STANDARD

% Characteristic strength of concrete in compression
f_ck = 35;                  %[MPa]

% Choice of EPDs (Relevant for density)
EPD_c = 'c35R';
EPD_s = 'rp';

% Height
h = 240;                    %[mm]

% Width
w = 360;                    %[mm]

% Longitudinal reinforcement

    % Diameter
    o_L = [16, 16, 16];     %[mm]

    % Number of bars in each row
    n_L = [5, 2, 5];        %[-]

    % Depth of middle of row
    d_L = 145;              %[mm]

% Preliminary calculations with inputs
M = def_c(f_ck, h, w, o_L, n_L, d_L);


% ADDITIONAL

% Transverse reinforcement

    % Diamter
    o_T = 6;                %[mm]

    % Width of middle hoop
    w_T = 200;              %[mm]

    % Spacing between group of hoops
    s_T = 150;              %[mm]


% --------------------------- Areas and masses ---------------------------

% Outputs area of longitudinal reinforcement
fprintf('Longitudinal reinforcement: Total area = %.4f [mm^2] \n', ...
    M.A_s*1E6);

% Calculates and outputs volumes

    % Longitudinal
    V_L = M.A_s * M.L;                                      %[m^3]
    fprintf('Longitudinal reinforcement: Volume = %.4e [m^3] \n', V_L);

    % Density
    fprintf('Longitudinal reinforcement: Density = %.2f [kg/m^3] \n', ...
        MI.(EPD_s).rho);

    % Mass
    m_L = V_L * MI.(EPD_s).rho;                             %[kg]
    fprintf('Longitudinal reinforcement: Mass = %.2f [kg] \n', ...
        m_L);


    % Transverse
    % (Assuming that the hoops are square)

        % Calculates the area of one transverse bars
        A_1T = pi * (o_T*1E-3/2)^2;                         %[m^2]
        fprintf(['Transverse reinforcement: Area of one bar = %.4f ' ...
            '[mm^2] \n'], A_1T*1E6);

        % Length of hoops per set
        l_T = 4*(M.h-2*M.c) + 2*(M.w-2*M.c) + 2*w_T*1E-3;   %[m]
        fprintf(['Transverse reinforcement: Total length of one group ' ...
            'of hoops = %.4f [mm] \n'], l_T*1E3);

        % Number of hoops groups
        n_T = M.L/(s_T*1E-3);                               %[-]
        fprintf(['Transverse reinforcement: Number of hoops = %.1f ' ...
            '[-] \n'], n_T);

        % Volumes
        V_T = A_1T * l_T * n_T;                             %[m^3]
        fprintf('Transverse reinforcement: Volume = %.4e [m^3] \n', V_T);
    
        % Density
        fprintf('Transverse reinforcement: Density = %.2f [kg/m^3] \n', ...
            MI.(EPD_s).rho);
    
        % Mass
        m_T = V_T * MI.(EPD_s).rho;                         %[kg]
        fprintf('Transverse reinforcement: Mass = %.2f [kg] \n', ...
            m_T);

    % Total reinforcement mass
    m_s = m_L + m_T;                                        %[kg]
    fprintf('Total reinforcement: Mass = %.2f [kg] \n', ...
        m_s);

    % GWP-GHG
    fprintf(['Reinforcement products: GWP-GHG = %.2f ' ...
        '[kg CO2-eq / m^3] \n'], MI.(EPD_s).GWP);
    

    % Concrete
    V_c = M.w*M.h*M.L - (V_T + V_L);                        %[m^3]
    fprintf('Concrete: Volume (excluding steel) = %.4f [m^3] \n', V_c);
    
    % Density
    fprintf('Concrete: Density = %.2f [kg/m^3] \n', MI.(EPD_c).rho);

    % Mass
    m_c = V_c * MI.(EPD_c).rho;                             %[kg]
    fprintf('Concrete: Mass = %.2f [kg] \n', m_c);

    % GWP-GHG
    fprintf('Concrete: GWP-GHG = %.2f [kg CO2-eq / m^3] \n', ...
        MI.(EPD_c).GWP);


% -------------------- Potential environmental impact --------------------

% Concrete
Q_c = V_c*MI.(EPD_c).GWP;                                   %[kg CO2-eq]
fprintf('Concrete: GWP-GHG = %.2f [kg CO2-eq] \n', Q_c);

% Concrete
Q_s = (V_L + V_T)*MI.(EPD_s).GWP;                           %[kg CO2-eq]
fprintf('Steel: GWP-GHG = %.2f [kg CO2-eq] \n\n', Q_s);

% Concrete
Q_SE01 = Q_c + Q_s;                                         %[kg CO2-eq]
fprintf('Column SE01: GWP-GHG = %.2f [kg CO2-eq] \n\n', Q_SE01);


% ----------------------- Geometric imperfections -----------------------

M.e_0 = geo_imperfections(0, M);                            %[m]
fprintf('Eccentricity of the axial: e_0 = %.2f [mm] \n\n', M.e_0*1E3);


% ---------------------- Applied first order moment ----------------------

M = bending_moment(0, M);                                   %[Nm]
fprintf('Applied first order moment: M_0Ed = %.2f [kNm] \n\n', ...
    M.M_0Ed*1E-3);


% --------------------- Effective creep coefficient ---------------------

M.varphi_ef = effective_creep(1, M); fprintf('\n\n')


% ------------------ Checking for second order effects ------------------

soe = check_soe(1, M); fprintf('\n\n')

% --------------------------- Moment capacity ---------------------------

M.M_Rd = moment_capacity(1, M); fprintf('\n\n')

% Printing result
fprintf('SE01 utilisation factor: %.3f \n\n', M.M_0Ed/M.M_Rd)

% ---------------------------- Input control ----------------------------

% Controlling input if desired
if IC == true
    input_control(0.03, M);
end