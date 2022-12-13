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
h = 220;                    %[mm]

% Width
w = 3600;                    %[mm]

% Longitudinal reinforcement

    % Diameter
    o_L = [6, 6];           %[mm]
    
    % Spacing
    s_L = 150;              %[mm]

    % Number of bars in each row
    n_L = [w/s_L, w/s_L];   %[-]

% Preliminary calculations with inputs
M = def_w(f_ck, h, w, o_L, n_L, []);


% ADDITIONAL

% Transverse reinforcement

    % Diamter
    o_T = 6;                %[mm]

    % Spacing between group of hoops
    s_T = 150;              %[mm]


% --------------------------- Areas and masses ---------------------------

% Outputs area of longitudinal reinforcement
fprintf('Longitudinal reinforcement: Total area = %.4f [mm^2] \n', ...
    M.A_s*1E6);

% Calculates and outputs volumes

    % Longitudinal

        % Number of bars
        n_L = w/s_L;                                            %[-]
        fprintf('Longitudinal reinforcement: Number = %.2f [m^3] \n', n_L);
    
        % Area
        A_L = M.A_s;                                            %[m^2]
        fprintf('Longitudinal reinforcement: Area = %.4e [m^3] \n', A_L);
    
        % Volumen
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

        % Number of bars
        n_T = M.L*1E3/s_T;                                      %[-]
        fprintf('Transverse reinforcement: Number = %.2f [m^3] \n', n_L);
    
        % Area
        % (Remember there is two members, multiply by two)
        A_T = 2 * (n_T * pi * (o_T*1E-3/2)^2);                  %[m^2]
        fprintf('Transverse reinforcement: Area = %.4e [m^3] \n', A_T);
    
        % Volumen
        V_T = A_T * M.w;            %[m^3]
        fprintf('Transverse reinforcement: Volume = %.4e [m^3] \n', V_T);
    
        % Density
        fprintf('Transverse reinforcement: Density = %.2f [kg/m^3] \n', ...
            MI.(EPD_s).rho);
    
        % Mass
        m_T = V_T * MI.(EPD_s).rho;                             %[kg]
        fprintf('Transverse reinforcement: Mass = %.2f [kg] \n', ...
            m_T);

    % Total reinforcement voluen
    V_s = V_L + V_T;                                        %[kg]
    fprintf('Total reinforcement: Volumen = %.2e [m^3] \n', ...
        V_s);

    % Total reinforcement mass
    m_s = m_L + m_T;                                        %[kg]
    fprintf('Total reinforcement: Mass = %.2f [kg] \n', ...
        m_s);
    
    % Concrete
    V_c = M.w*M.h*M.L - (V_T + V_L);                        %[m^3]
    fprintf('Concrete: Volume (excluding steel) = %.4f [m^3] \n', V_c);
    
    % Density
    fprintf('Concrete: Density = %.2f [kg/m^3] \n', MI.(EPD_c).rho);

    % Mass
    m_c = V_c * MI.(EPD_c).rho;                             %[kg]
    fprintf('Concrete: Mass = %.2f [kg] \n', m_c);

    % Potential environmental impact per volumen (GWP-GHG)

        % Reinforcement (Steel)
        fprintf(['Reinforcement products: GWP-GHG per volumen  = %.2f ' ...
            '[kg CO2-eq / m^3] \n'], MI.(EPD_s).GWP);

        % Concrete
        fprintf('Concrete: GWP-GHG per volumen = %.2f [kg CO2-eq / m^3] \n', ...
            MI.(EPD_c).GWP);

% -------------------- Potential environmental impact --------------------

% Concrete
Q_c = V_c*MI.(EPD_c).GWP;                                   %[kg CO2-eq]
fprintf('Concrete: GWP-GHG = %.2f [kg CO2-eq] \n', Q_c);

% Concrete
Q_s = (V_L + V_T)*MI.(EPD_s).GWP;                           %[kg CO2-eq]
fprintf('Steel: GWP-GHG = %.2f [kg CO2-eq] \n\n', Q_s);

% Concrete
Q_VES01 = Q_c + Q_s;                                        %[kg CO2-eq]
fprintf('Wall VES01: GWP-GHG = %.2f [kg CO2-eq] \n\n', Q_VES01);


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


% ------------------- Calculating second order effects -------------------

% M_Ed_sm = simplified_method(1, M); fprintf('\n')
% M_Ed_nc = nominal_curvature(1, M); fprintf('\n')
% M_Ed_ns = nominal_stiffness(1, 1.5, M); fprintf('\n\n')

% --------------------------- Moment capacity ---------------------------

M.M_Rd = moment_capacity(1, M); fprintf('\n\n')

% Printing result
fprintf('VES01 utilisation factor: %.3f \n\n', M.M_0Ed/M.M_Rd)

% ---------------------------- Input control ----------------------------

% Controlling input if desired
if IC == true
    input_control(0.03, M);
end