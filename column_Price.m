%% Initializing MATLAB
clear; clc; close all;

% Adding function folder to path
addpath('functions')

% Loading climate information and plot colors
MI = material_info(1);   PC = plot_colors();


%% User-input

% Number of points (Only even numbers)
np = 200;

% Limits for side length and longitunal reinforcement diameter
o_L_min = 4;                            o_L_max = 26;   %[mm]
s_min = max(150, 2*(o_L_max+15+6));     s_max = 300;    %[mm]

% Characteristic compressive strength of concrete
% (16, 20, 25, 30 [F/R]  |  35, 40 [L/R])
f_ck = 16;   EPD_c = 'c16F';   EPD_s = 'rp';

% Number of longitunal reinforcement rows
n_L = [2, 2];

% Printing
printing = 1;


%% Analysis of calc_Price function

% Points in intervals
o_L = linspace(o_L_min, o_L_max, np);
s = linspace(s_min, s_max, np);

% Rreallocating arrays and starting counter
ss = zeros(np, np);   oo=ss;   PP=ss;

% Loop for diameter of longitunal reinforcement
for x = 1:size(o_L, 2)

    % Loop for side length
    for y = 1:size(s, 2)

        % Definitions
        M = def_c(f_ck, s(y), s(y), [o_L(x), o_L(x)], n_L, []);

        % Price
        PP(x,y) = calc_Price(0, M, EPD_c, EPD_s);

        % Storing result
        ss(x,y) = s(y);   oo(x,y) = o_L(x);

    end

end


%--------------------------------------------------------------------------
%% Plotting the price for every combination

% Creating figure for the future plot
p1 = figure('Name', 'Price', 'PaperUnits', 'centimeters', ...
    'PaperSize', [21 19]);

% Making it the current figure and insuring that future enteries will be
% added to it.
figure(p1);   hold on

% Plotting the price as a surface plot
surf(ss, oo, PP ,"EdgeColor", "k", "EdgeAlpha", 0.5, ...
    "HandleVisibility", "off");

% Adding plot title
title('', ...
    sprintf('(L = %.f [mm], f_c_k = %.f, rho_s = %.f [kg/m^3])', ...
    M.L*10^3, M.f_ck*10^(-6), MI.(EPD_s).rho));

% Changing plot legends
xlabel('Side length of cross section, {\it s} [mm]')
ylabel('Diameter of longitunal reinforcement, {\it o_L} [mm]')
zlabel('Price of column, [DKK]')

% Concluding plot
conclude_plot(p1, 'right')


%--------------------------------------------------------------------------
%% Plotting prices as a function of side length

% Creating new figure for the new plot
p2 = figure('Name', 'Side length', 'PaperUnits', 'centimeters', ...
    'PaperSize', [21 8]);

% Making it the current figure and insuring that future enteries will be
% added to it.
figure(p2); hold on;

% Plotting the results
plot(ss(1,:)', PP(size(PP,2)/2, :)', "HandleVisibility", "off", ...
    'Color', PC.gwp, 'Marker', '.', 'LineStyle', 'none', ...
    'MarkerSize', 10);

% Adding plot title
title('', ...
    sprintf(['(o_L = %.2f [mm], L = %.f [mm], f_c_k = %.f, ' ...
    'rho_s = %.f [kg/m^3])'], ...
    oo(size(PP,2)/2, 1), M.L*10^3, M.f_ck*10^(-6), MI.(EPD_s).rho));

% Changing plot legends
xlabel('Side length of cross section,{\it s} [mm]')
ylabel('Price of column, [DKK]')

% Changing plot position
plot_position(p2, 'left');

% Adding grid lines
grid on

% Simply concluding plot
hold off


%--------------------------------------------------------------------------
%% Plotting prices as a function of diameter

% Creating new figure for the new plot
p3 = figure('Name', 'Diameter', 'PaperUnits', 'centimeters', ...
    'PaperSize', [21 8]);

% Making it the current figure and insuring that future enteries will be
% added to it.
figure(p3); hold on;

% Plotting the results
plot(oo(:,1), PP(:, size(PP,2)/2), "HandleVisibility", "off", ...
    'Color', PC.gwp, 'Marker', '.', 'LineStyle', 'none', ...
    'MarkerSize', 10);

% Adding plot title
title('', ...
    sprintf(['(s = %.2f [mm], L = %.f [mm], f_c_k = %.f, ' ...
    'rho_s = %.f [kg/m^3])'], ...
    ss(1, size(PP,2)/2), M.L*10^3, M.f_ck*10^(-6), MI.(EPD_s).rho));

% Changing plot legends
xlabel('Diameter of longitunal reinforcement, {\it o_L} [mm]')
ylabel('Price of column, [DKK]')

% Changing plot position
plot_position(p3, 'top_left');

% Adding grid lines
grid on

% Simply concluding plot
hold off


%--------------------------------------------------------------------------
% Printing (if wished for)

if printing == 1

% Path to the results folder
path_folder = '.\results\Pricing\Column';

% Figure array of currently open figures
fa = findall(groot, 'Type', 'Figure', '-not', 'Name', 'Input check');

% Sortint the figure array
[~, s_idx] = sort([fa.Number]);
fa = fa(s_idx);

% Number of figures
number_figure = 1:size(fa ,1);

% Name of the figures
name_file = append('L', string(M.L*10^(3)), '_nL', join(string(M.n_L),'_'), ...
    '_oL', sprintf('%.f',M.o_L(end)*10^(3)), ...
    '_C', sprintf('%.f_',M.f_ck*10^(-6)) , string(get(fa, 'Name')));

% Printing both a PDF and EPS
print_figures(path_folder, name_file, number_figure);

end


%--------------------------------------------------------------------------
% Concluding plot

function conclude_plot(plot_identifier, screen_side)

% Adding colorbar and grid lines (+ chaning the colormap)
colorbar;   grid on;   colormap gray;

% Changing rotation
view(300, 25);

% Changing plot position
plot_position(plot_identifier, screen_side);

% Concluding creation of the plot
hold off

end


%--------------------------------------------------------------------------
% Plot position

function plot_position(plot_identifier, screen_side)

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

elseif strcmp(screen_side, 'top_left')

    % To the left side
    set(plot_identifier, 'Position', ...
        [0 th+si(3)*0.5*hwr si(3)*0.5 si(3)*0.5*hwr])

end

end


%--------------------------------------------------------------------------
% Changing figure properties

function change_prop(fig)

% Loading plot colors
PC = plot_colors();

% Finding the 'Axes'-objects
ax = findall(fig, 'Type', 'Axes');

    % Changing its properties
    set(ax, 'Position', [0.1, 0.1, 1, 1])
    set(ax, 'OuterPosition', [0, 0.1, 0.9, 0.9])


% Finding the 'ColorBar'-objects
cb = findall(fig, 'Type', 'ColorBar');

    % Changing its properties
    set(cb, 'Position', [0.9, 0.1, 0.03, 0.8])


% Finding the legend object
lo = findall(fig, 'Type', 'Text');

    % Changing properties
    set(lo(2), "Units", "normalized")
    set(lo(2), "Position", [0.5, 1.025, 0.5])
    set(lo(2), "FontSize", 11)
    set(lo(2), "Color", [0.15,0.15,0.15])


end