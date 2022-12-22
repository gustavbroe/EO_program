%% Initializing MATLAB
clear; clc; close all hidden;

% Adding function folder to path
addpath('functions')

% Loading climate information, plot colors, and baseline values
MI = material_info(1);   PC = plot_colors();   BL = baseline();

%--------------------------------------------------------------------------
%% User-input

% Number of points
np = 60;

% Characteristic compressive strength of concrete
% (16, 20, 25, 30 [F/R]  |  35, 40 [L/R])
f_ck = 16;   EPD_c = 'c16F';   EPD_s = 'rp';            %[MPa]

% Width of cross section
w = 3600;                                               %[mm]

% Spacing between longitudinal reinforcement
s_L = 150;                                              %[mm]

% Limits for thickness and longitunal reinforcement diameter
o_L_min = 6;                o_L_max = 16;               %[mm]
h_min = 80;                 h_max = 250;                %[mm]

% Show first order and second order moments (transparent planes)
show_planes = 0;

% Printing of figures
printing = 1;

% Input control of lowest GWP solution?
control_gwp = 1;

% Consider pricing
con_price = 1;

    % Input control of cheapest wall?
    control_price = 0;


%--------------------------------------------------------------------------
%% Calculating applied moment, moment capacity and GWP for each combination

% Number of longitunal reinforcement rows
n_L = [w/s_L, w/s_L];

% Points in intervals
o_L = linspace(o_L_min, o_L_max, np);
h = linspace(h_min, h_max, np);

% Rreallocating arrays and starting counter
hh = zeros(np, np);   oo = hh;   MM = hh;   rr = hh;
SM = hh;   NC = hh;   NS = hh;   nn = hh;   GG = hh;   PP = hh;

% Loop for diameter of longitunal reinforcement
for x = 1:size(o_L, 2)

    % Loop for thickness
    for y = 1:size(h, 2)

        % Definitions
        M = def_w(f_ck, h(y), w, [o_L(x), o_L(x)], n_L, []);

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
        hh(x,y) = h(y);   oo(x,y) = o_L(x);

    end

end


%--------------------------------------------------------------------------
%% Plotting moment related results in right plot

% Printing/Paper properties for the figures
fig_prop = {'PaperUnits', 'centimeters', 'PaperSize', [21 19]};

% Creating figure for the future plot
p1 = figure('Name', 'Moment', fig_prop{:});

% Making it the current figure and insuring that future enteries will be
% added to it.
figure(p1);   hold on

% Plotting the moment capacity as a surface plot
surf(hh, oo, MM ,"EdgeColor", "k", "EdgeAlpha", 0.5, ...
    "HandleVisibility", "off");

% Plotting first order moment, M_0Ed
if show_planes == 1
    surf(hh, oo, repmat(M.M_0Ed,np,np) ,"EdgeColor", "k", ...
        "EdgeAlpha", 0.2, "FaceAlpha", 0.3 ,"FaceColor", PC.gray, ...
        "DisplayName", "{\it M_0_E_d}", "HandleVisibility", 'on');
end

% Creating logical array and plotting OVER REINFORCED points
plot3(hh(rr == 1), oo(rr == 1), MM(rr == 1), "MarkerSize", 4, ...
    "Marker", "*", "Color", PC.red, "LineStyle", "none", ...
    "DisplayName", "Over-reinforced");

% Creating logical array and plotting NON-REINFORCED points
plot3(hh(nn == 1), oo(nn == 1), MM(nn == 1), "MarkerSize", 7, ...
    "Marker", "o", "Color", PC.red, "LineStyle", "none", ...
    "DisplayName", "Non-reinforced");


% Defining text for each method
tx = {'Simplified method', 'Nonimal curvature', 'Nominal stiffness'};

% Plotting the planes and intersection line between the three planes and MM
d = {SM, NC, NS};    c = {PC.sm, PC.nc, PC.ns};   vr = {0,0,0};
for k = 1:3

    % Plotting the planes
    if show_planes == 1
        surf(hh, oo, d{k} ,"EdgeColor", "k", "EdgeAlpha", 0.2, ...
            "FaceAlpha", 0.3 ,"FaceColor", c{k}, ...
            "HandleVisibility", 'off');
    end

    % Logical array of valid cross sections with adequate moment cap.
    vcs = MM > d{k};

    % Making two copies
    vcsR = vcs;   vcsC = vcs;

    % Looping through the wallS of the logical array
    for i = 1:size(vcs, 2)

        % Checking if there are more than 1 1 in the wall
        if ~isempty(find(vcs(:,i) == 1, 1))

            % Index vector for the rows equal 1
            j = find(vcs(:,i) == 1);
    
            % Changen all but the first 1 to 0
            vcsC(j(2:end),i) = false;

        end

        % Checking if there are more than 1 1 in the wall
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

% Plotting the valid cross section just within the limit
plot3(hh(vcs), oo(vcs), MM(vcs), "MarkerSize", 16, "Marker", ".", ...
    "Color", c{k}, "LineStyle", "none", "DisplayName", tx{k});

% Ending for-loop
end


% Adding plot title
title('', ...
    sprintf(['(N_E_d = %.f [kN] / f_c_k = %.f [MPa] / L_0 = %.2f [m]' ...
    ' / c+o_T = %.f [mm])'] ...
    , M.N_Ed*10^(-3), M.f_ck*10^(-6), M.L_0, (M.c+M.o_T)*10^(3) ));

% Changing plot legends
xlabel('Thickness of wall, {\it h} [mm]')
ylabel('Diameter of longitunal reinforcement, {\it o_L} [mm]')
zlabel('Design value of the bending moment, {\it M} [Nm]')

% Concluding plot
conclude_plot(p1, 'right')


%--------------------------------------------------------------------------
%% Plotting GWP results in left plot

% Creating new figure for the new plot
p2 = figure('Name', 'GWP', fig_prop{:});

% Making it the current figure and insuring that future enteries will be
% added to it.
figure(p2); hold on;

% Plotting the results
surf(hh, oo, GG ,"EdgeColor", "k", "EdgeAlpha", 0.5, ...
    "HandleVisibility", "off");

% Adding plot title
title('', ...
    sprintf(['(N_E_d = %.f [kN] / f_c_k = %.f [MPa] / L_0 = %.2f [m]' ...
    ' / c+o_T = %.f [mm])'] ...
    , M.N_Ed*10^(-3), M.f_ck*10^(-6), M.L_0, (M.c+M.o_T)*10^(3) ));

% Changing plot legends
xlabel('Thickness of wall,{\it h} [mm]')
ylabel('Diameter of longitunal reinforcement,{\it o_L} [mm]')
zlabel('Global Warming Potential,{\it GWP} [kg CO_2-eq]')

% Conclude plot
conclude_plot(p2, 'left');


%--------------------------------------------------------------------------
%% Plotting PP results in another left plot

if con_price == 1

% Creating new figure for the new plot
p3 = figure('Name', 'Price', fig_prop{:});

% Making it the current figure and insuring that future enteries will be
% added to it.
figure(p3); hold on;

% Plotting the results
surf(hh, oo, PP ,"EdgeColor", "k", "EdgeAlpha", 0.5, ...
    "HandleVisibility", "off");

% Adding plot title
title('', ...
    sprintf(['(N_E_d = %.f [kN] / f_c_k = %.f [MPa] / L_0 = %.2f [m]' ...
    ' / c+o_T = %.f [mm])'] ...
    , M.N_Ed*10^(-3), M.f_ck*10^(-6), M.L_0, (M.c+M.o_T)*10^(3) ));

% Changing plot legends
xlabel('Thickness of wall,{\it h} [mm]')
ylabel('Diameter of longitunal reinforcement,{\it o_L} [mm]')
zlabel('Price of wall, [DKK]')

% Setting axes equal to the limits
axis([h_min h_max o_L_min o_L_max])

% Conclude plot
conclude_plot(p3, 'left');

end

%--------------------------------------------------------------------------
%% Outputting and plotting the lowest GWP solution (all plots)

% Printing header
fprintf('The lowest possible GWP for each method is: \n\n')

% Finding index from highest to lowest method
[~, ii] = sort([min(GG(vr{1})), min(GG(vr{2})), min(GG(vr{3}))]);
for i = ii

    % Plotting the valid cross section just within the limit in price plot
    if con_price == 1
    figure(p3); hold on;
    plot3(hh(vr{i}), oo(vr{i}), PP(vr{i}), "MarkerSize", 16, "Marker", ...
        ".", "Color", c{i}, "LineStyle", "none", "DisplayName", tx{i});
    figure(p3); hold off;   
    end
    
    % Plotting the valid cross section just within the limit in GWP plot
    figure(p2); hold on;
    plot3(hh(vr{i}), oo(vr{i}), GG(vr{i}), "MarkerSize", 16, "Marker", ...
        ".", "Color", c{i}, "LineStyle", "none", "DisplayName", tx{i});
    figure(p2); hold off;   
    
    if con_price == 1
    figure(p3); hold on;
    end

    % Outputting the lowest possible GWP
    [row, col] = find(GG == min(GG(vr{i})));
    fprintf('%s = %.2f [kg CO2-ep] / %.2f [DKK] / %.2f [%%] ', tx{i}, ...
        GG(row,col), PP(row,col), ...
        (BL.wall.GWP - GG(row(end),col(end)))/BL.wall.GWP*100);
    fprintf('| h = %.2f [mm] / o_L = %.2f [mm] \n', ...
        hh(row,col), oo(row,col));

    % Marking the point in price figure
    if con_price == 1
    plot3(hh(row,col), oo(row,col), PP(row,col), "MarkerSize", 18, ... 
        "Marker", "o", "Color", PC.green, "LineStyle", "none", ...
        "HandleVisibility","off", 'LineWidth', 2);
    figure(p3); hold off;  
    end
    
    % Marking the point in GWP figure
    figure(p2); hold on; 
    plot3(hh(row,col), oo(row,col), GG(row,col), "MarkerSize", 18, ... 
        "Marker", "o", "Color", PC.green, "LineStyle", "none", ...
        "HandleVisibility","off", 'LineWidth', 2);
    figure(p2); hold off;   figure(p1); hold on; 
    plot3(hh(row,col), oo(row,col), MM(row,col), "MarkerSize", 18, ...
        "Marker", "o", "Color", PC.green, "LineStyle", "none", ...
        "HandleVisibility","off", 'LineWidth', 2);
    figure(p1); hold off;   
    
    if con_price == 1
    figure(p3); hold on;
    end


    % Calling imput control if wished for
    if control_gwp == 1

        % Only doing it ones
        while control_gwp == 1

            % Calling function
            fig_IC = input_control(0.03, def_w(f_ck, ...
                hh(row(end),col(end)), w, ...
                [oo(row(end),col(end)), oo(row(end),col(end))], n_L, []));

            control_gwp = 0;

        end

    % Ending if-loop
    end

end


%--------------------------------------------------------------------------
%% Outputting cheapest wall and plotting it (all plots)

if con_price == 1

% Printing header
fprintf('\n\nThe cheapest wall for each method is: \n\n')

% Finding index from highest to lowest method
[~, ii] = sort([min(PP(vr{1})), min(PP(vr{2})), min(PP(vr{3}))]);
for i = ii

    % Outputting the lowest possible GWP
    [row, col] = find(PP == min(PP(vr{i})));
    fprintf('%s = %.2f [DKK] / %.2f [kg CO2-eq] / %.2f [%%] ', tx{i}, ...
        PP(row(end),col(end)), GG(row(end),col(end)), ...
        (BL.wall.GWP - GG(row(end),col(end)))/BL.wall.GWP*100);
    fprintf('| h = %.2f [mm] / o_L = %.2f [mm] \n', ...
        hh(row(end),col(end)), oo(row(end),col(end)));

    % Marking the point in both figures
    plot3(hh(row(end),col(end)), oo(row(end),col(end)), ...
        PP(row(end),col(end)), "MarkerSize", 18, ... 
        "Marker", "o", "Color", PC.yellow, "LineStyle", "none", ...
        "HandleVisibility","off", 'LineWidth', 2);
    figure(p3); hold off;   figure(p2); hold on; 
    plot3(hh(row(end),col(end)), oo(row(end),col(end)), ...
        GG(row(end),col(end)), "MarkerSize", 18, ... 
        "Marker", "o", "Color", PC.yellow, "LineStyle", "none", ...
        "HandleVisibility","off", 'LineWidth', 2);
    figure(p2); hold off;   figure(p1); hold on; 
    plot3(hh(row(end),col(end)), oo(row(end),col(end)), ...
        MM(row(end),col(end)), "MarkerSize", 18, ...
        "Marker", "o", "Color", PC.yellow, "LineStyle", "none", ...
        "HandleVisibility","off", 'LineWidth', 2);
    figure(p1); hold off;   figure(p3); hold on; 


    % Calling imput control if wished for
    if control_price == 1

        % Only doing it ones
        while control_price == 1

            % Calling function
            fig_IC = input_control(0.03, def_w(f_ck, ...
                hh(row(end),col(end)), w, ...
                [oo(row(end),col(end)), oo(row(end),col(end))], n_L, []));

            control_price = 0;

        end

    % Ending if-loop
    end

end

end

% Changing figure properties
change_prop(p1);   change_prop(p2);   

if con_price == 1
change_prop(p3);
end

%--------------------------------------------------------------------------
% Printing (if wished for)

if printing == 1

% Path to the results folder
path_folder = '.\results\Wall\Three_dimensional';

% Figure array of currently open figures
fa = findall(groot, 'Type', 'Figure', '-not', 'Name', 'Input Check: Wall');

% Sortint the figure array
[~, s_idx] = sort([fa.Number]);
fa = fa(s_idx);

% Number of figures
number_figure = 1:size(fa ,1);

% Name of the figures
name_file = append('np', string(np), ...
    '_h', sprintf('%.f_%.f', h_min, h_max), ...
    '_nL', join(string(M.n_L),'_'), ...
    '_oL', sprintf('%.f_%.f',o_L_min, o_L_max), ...
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


% Finding all the green points
gp = findall(fig, 'Color', PC.green);

    % Changing its properties
    set(gp(1), "DisplayName", "Lowest GWP for each method")
    set(gp(1), "HandleVisibility", "on")

% Finding all the yellow points
gp = findall(fig, 'Color', PC.yellow);

    % Changing its properties
    if ~isempty(gp)
    set(gp(1), "DisplayName", "Lowest price for each method")
    set(gp(1), "HandleVisibility", "on")
    end


% Finding the legend object
lo = findall(fig, 'Type', 'Legend');

    % Changing properties
    set(lo, "Orientation", "horizontal")
    set(lo, "Position", [0.5, 0.04, 0.01, 0.01])
    set(lo, "NumColumns", 3)


% Finding the legend object
lo = findall(fig, 'Type', 'Text');

    % Changing properties
    set(lo(2), "Units", "normalized")
    set(lo(2), "Position", [0.5, 1.025, 0.5])
    set(lo(2), "FontSize", 11)
    set(lo(2), "Color", [0.15,0.15,0.15])


end