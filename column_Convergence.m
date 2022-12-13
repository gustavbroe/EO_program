%% Initializing MATLAB
clear; clc; close all hidden;

% Adding function folder to path
addpath('functions')

% Loading climate information and plot colors
CI = material_info(1);   PC = plot_colors();

%--------------------------------------------------------------------------
% CATION
% The execution of this file takes a considerable amount of time:
% Approx. 10-20 minuttes.

%--------------------------------------------------------------------------
% User-input

% Parameter for convergence analysis (1 = GWP, 0 = Price)
CP = 1;

% Array with number of points
np_min = 20;                            np_max = 320;
np_points = np_min:10:np_max;

% Characteristic compressive strength of concrete
% (16, 20, 25, 30 [F/R]  |  35, 40 [L/R])
f_ck = 16;   EPD_c = 'c16F';   EPD_s = 'rp';

% Number of longitunal reinforcement rows
n_L = [2, 2];

% Limits for side length and longitunal reinforcement diameter
o_L_min = 6;                            o_L_max = 26;   %[mm]
s_min = max(200, 2*(o_L_max+20+6));     s_max = 350;    %[mm]

% Printing of figure
printing = 1;


%--------------------------------------------------------------------------
% Calculating applied moment, moment capacity and GWP for each combination

% Allocating space for results
res = zeros(size(np_points,2), 3);

% Starting counter
counter = 1;

% Looping through diffrent number of points on each axes
for np = np_points

% Points in intervals
o_L = linspace(o_L_min, o_L_max, np);
s = linspace(s_min, s_max, np);

% Rreallocating arrays and starting counter
ss = zeros(np, np);   oo=ss;   MM=ss;   rr=ss;
SM = ss;   NC = ss;   NS = ss;   nn = ss;   GG = ss;   PP = ss;

% Loop for diameter of longitunal reinforcement
for x = 1:size(o_L, 2)

    % Loop for side length
    for y = 1:size(s, 2)

        % Definitions
        M = def_c(f_ck, s(y), s(y), [o_L(x), o_L(x)], n_L, []);

        % Geometrical imperfection
        M.e_0 = geo_imperfections(0, M);
        
        % First order bending moment
        M = bending_moment(0, M);
    
        % Effects of creep
        M.varphi_ef = effective_creep(0, M);

        % Moment capacity
        [MM(x,y), rr(x,y)] = moment_capacity(0, M);

        % Applied second order moment
            
            %Simplified method II
            SM(x,y) = simplified_method(0, M);

            % Nominal curvature
            NC(x,y) = nominal_curvature(0, M);

            % Nominal stiffness
            [NS(x,y), nn(x,y)] = nominal_stiffness(0, 1.2, M);

        % Global warming potential
        GG(x,y) = calc_GWP(0, {EPD_c, EPD_s}, M);

        % Price
        PP(x,y) = calc_Price(0, M, EPD_c, EPD_s);

        % Storing result
        ss(x,y) = s(y);   oo(x,y) = o_L(x);

    end

end


% Finding the planes and intersection line between the three planes and MM
d = {SM, NC, NS};   vr = {0,0,0};
for k = 1:3

    % Logical array of valid cross sections with adequate moment cap.
    vcs = MM > d{k};

    % Making two copies
    vcsR = vcs;   vcsC = vcs;

    % Looping through the COLUMNS of the logical array
    for i = 1:size(vcs, 2)

        % Checking if there are more than 1 1 in the column
        if ~isempty(find(vcs(:,i) == 1, 1))

            % Index vector for the rows equal 1
            j = find(vcs(:,i) == 1);
    
            % Changen all but the first 1 to 0
            vcsC(j(2:end),i) = false;

        end

        % Checking if there are more than 1 1 in the column
        if ~isempty(find(vcs(i,:) == 1, 1))

            % Index vector for the rows equal 1
            j = find(vcs(i,:) == 1);
    
            % Changen all but the first 1 to 0
            vcsR(i,j(2:end)) = false;

        end

    % Ending loop
    end

    % Combining the two arrays
    vcs = logical(vcsR + vcsC);

    % Removing the that follow the edge and not the intersection
    j = find(vcs(1,:)==1);   vcs(1,j(2:end)) = false;
    j = find(vcs(:,1)==1);   vcs(j(2:end),1) = false;

    % Storing the results
    vr{k} = vcs;

% Ending for-loop
end


% Defining text for each method
tx = {' - Simplified method', ' - Nonimal curvature', ...
    ' - Nominal stiffness'};

% Post processing for GWP
if CP == 1

    % Storring the results (outside the loop)
    res(counter, :) = [min(GG(vr{1})), min(GG(vr{2})), min(GG(vr{3}))];
    
    %--------------------------------------------------------------------------
    % Outputting results for the current loop
    
    % Printing header
    fprintf('\nThe lowest GWP cross section: (%g/%g points) \n', np, np^2)
    
    % Finding index from highest to lowest method
    [~, ii] = sort([res(counter, :)]);
    for i = ii
    
        % Outputting the lowest possible GWP
        [row, col] = find(GG == min(GG(vr{i})));
        fprintf('%s = %.2f [kg CO2-ep] ', tx{i}, GG(row,col));
        fprintf('| s = %.2f [mm] / o_L = %.2f [mm] \n', ...
            ss(row,col), oo(row,col));
    
    end


% Post processing for Price
elseif CP == 0

    % Storring the results (outside the loop)
    res(counter, :) = [min(PP(vr{1})), min(PP(vr{2})), min(PP(vr{3}))];
    
    %--------------------------------------------------------------------------
    % Outputting results for the current loop
    
    % Printing header
    fprintf('\nThe lowest Price cross section: (%g/%g points) \n', np, np^2)
    
    % Finding index from highest to lowest method
    [~, ii] = sort([min(PP(vr{1})), min(PP(vr{2})), min(PP(vr{3}))]);
    for i = ii
    
        % Outputting the lowest possible GWP
        [row, col] = find(PP == min(PP(vr{i})));
        fprintf('%s = %.2f [DKK] / %.2f [kg CO2-eq] ', tx{i}, ...
            PP(row(end),col(end)), GG(row(end),col(end)));
        fprintf('| s = %.2f [mm] / o_L = %.2f [mm] \n', ...
            ss(row(end),col(end)), oo(row(end),col(end)));
    
    end

% Ending post processing
end


% Adding to the counter
counter = counter + 1;

% Ending the number of points for-loop
end


%--------------------------------------------------------------------------
% Plotting the relative error as a function of number of points

% Creating figure for the future plot
p1 = figure('Name', 'Absolute error', "PaperSize", [21, 12]);

% Making it the current figure and insuring that future enteries will be
% added to it.
figure(p1);   hold on

% Plot colors
plot_colors = [PC.sm, PC.nc, PC.ns];

% Calculating the absolute error
% (Further work: Calculate the ANALYTICAL solution for each method)
abs_error = abs(res - res(end, :));             %[kg CO2-ep]

% The relative error
rel_error = abs_error./res(end, :);             %[-]

% Loop through each method
for i = 1:size(tx,2)

    % Plotting the relative error for each method
    plot(np_points.^2, rel_error(:,i)*100, "Color", plot_colors(i), ...
        "DisplayName", tx{i}, "Marker", "o", "MarkerSize", 9, ...
        "LineStyle", "-", "LineWidth", 2);

end

% Changing plot legends
xlabel('Number of points, {\it n_{points}} [-]')
ylabel('Relative error, {\it error_{rel}} [%]')

% Concluding creation of the plot
conclude_plot(p1, 'right');


%--------------------------------------------------------------------------
% Plotting the error as a function of number of points

% Creating figure for the future plot
p2 = figure('Name', 'Error', "PaperSize", [21, 12]);

% Making it the current figure and insuring that future enteries will be
% added to it.
figure(p2);   hold on

% Calculating the error
% (Further work: Calculate the ANALYTICAL solution for each method)
error = res - res(end, :);             %[kg CO2-ep]

% Loop through each method
for i = 1:size(tx,2)

    % Plotting the relative error for each method
    plot(np_points.^2, error(:,i), "Color", plot_colors(i), ...
        "DisplayName", tx{i}, "Marker", "o", "MarkerSize", 9, ...
        "LineStyle", "-", "LineWidth", 2);

end

% Changing plot legends
xlabel('Number of points, {\it n_{points}} [-]')
ylabel('Error, {\it error} [kg CO_{2}-eq.]')

% Concluding creation of the plot
conclude_plot(p2, 'left');

%--------------------------------------------------------------------------
% Printing

if printing == 1

% Path to the results folder
path_folder = '.\results\Column\Convergence';

% Number of figures
number_figures = 1;

% Name of the figures
name_file = append('np', string(np_min), '_', string(np_max), ...
    '_nL', join(string([2,2]),'_'), '_C', sprintf('%.f_',f_ck), ...
    'Convergence');

% Printing both a PDF and EPS
print_figures(path_folder, name_file, number_figures);

end


%--------------------------------------------------------------------------
% Concluding plot

function conclude_plot(plot_identifier, screen_side)

% Adding grid lines
grid on

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

% Finding the 'Axes'-objects
ax = findall(plot_identifier, 'Type', 'Axes');

    % Adding minor grid lines for both axis
    set(ax, 'XMinorGrid', 'on');
    set(ax, 'YMinorGrid', 'on');

% Inserting legends
legend;

% Concluding creation of the plot
hold off

end
