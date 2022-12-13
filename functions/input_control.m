function fig = input_control(margin, definitions)

%--------------------------------------------------------------------------
% This function check if the input from the definition function (def_1...)
% is reasonable, by highlighting parameters and variables with color 
% and plotting the cross section.
%
%    (Window's position in parentheses are for walls, w > 4h)
%
%  - Top (Top left) window compares the input to regulations from EC2.
%
%  - Bottom left (Top right) window uses the 'plot_cross_section'-function 
%    to plot the cross section.
%
%  - Bottom right (Bottom) window displays steel, exposure, and control 
%    class from partial factors and mechanical properties, also 
%    highlighting values outside the limits given in EC2.
%
%
% Input:    A parent object (Only, Axes, UIAxes, etc.)
%           Information about the member (structure array) 
%           A number for the desired distance around the concrete edge.
% Output:   The generated plot as an object.
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Translates user input

    % Margin in figure 
    % (precent of figure width allocated to left and right margin)
    mar = margin;

    % Importing definitions
    M = definitions;

%--------------------------------------------------------------------------
%% Program

% Checking if the member is a wall
if M.w > 4*M.h
    is_wall = logical(1);
else
    is_wall = logical(0);
end

% Determining the figure name
if is_wall
    fig_name = 'Input Check: Wall';
else
    fig_name = 'Input Check: Column';
end

% Creating new figure
fig = uifigure('Name', fig_name, 'Resize', 'on'); 

% Gathering display information

    % Current screen size
    si = get(0,'screensize');

    % Height of toolbar
    th = si(4)*0.06;

    %Height plotting window header
    wh = si(4)*0.18;

% Changing the position of the plot

    % Height and witdh of the figure
    hof = si(4)-(th+wh);
    wof = si(3)*0.5;
    if is_wall 
        wof = si(3);
    end

    % Changing from default
    fig.Position = [0 th wof hof];


% Creating TOP panel

    % Creating uipanel in figure
    tp = uipanel(fig);

    % Size of margin in figure
    som = mar*hof;

    % Changing its position
    tp.Position = [som hof/2+som wof-som*2 hof/2-som*2];
    if is_wall
        % Top left
        tp.Position = [som hof/2+som wof*0.6-som*2 hof/2-som*2];
    end

    % Giving the panel a title

        % Inserting top panel text label
        tptl = uilabel(tp, "FontSize", 14);

        % Changing size
        tptl.Position = [som hof/2-som*4 wof-som*2 som*2];

        % Specifying text
        tptl.Text = ['<strong>Requirements given by the Eurocode 2' ...
            '</strong>'];

        % Setting interpreter
        tptl.Interpreter = 'html';


    % Creating table
    
        % Specifying column names
        vars = {'Name', 'Symbol', 'Unit', 'Value', 'Input'};
    
        % Gathering data from function
        data = table_data(M, is_wall);
    
        % Comparing values and assigning colors
        rc = table_style(data);
    
        % Formatting data
        data(:, end-1:end) = arrayfun(@(x)sprintf('%.3e',x), ...
            cell2mat(data(:, end-1:end)), 'UniformOutput', false);
    
        % Creating uitable in panel
        ttab = uitable(tp, "Data", data, ...
            "ColumnName", vars, "RowName", '');

        % Changing position
        ttab.Position = [som, som, wof-som*4, hof/2-som*5];
        if is_wall
            ttab.Position = [som, som, wof*0.6-som*4, hof/2-som*5];
        end

    
        % Changing row colors
        for i = 1:size(data, 1)
            addStyle(ttab, uistyle('BackgroundColor', rc{i}), ...
                'row', i);
        end


%--------------------------------------------------------------------------


% Creating LEFT BOTTOM panel

    % Creating uipanel in figure
    blp = uipanel(fig);

    % Panel dimensions
    wop = wof/2-som*1.5;
    hop = hof/2-som*1;
    if is_wall
        % Bottom
        wop = wof - som*2;
        hop = hof/2 - som*1;
    end

    % Changing its position
    blp.Position = [som som wop hop];


    % Giving the panel a title

        % Inserting top panel text label
        bptl = uilabel(blp, "FontSize", 14);

        % Changing size
        bptl.Position = [som hof/2-som*3 wof-som*2 som*2];

        % Specifying text
        bptl.Text = ['<strong>Visual representation of cross section' ...
            '</strong>'];

        % Setting interpreter
        bptl.Interpreter = 'html';


    % Inserting cross section plot

        % Creating uiaxes
        cp = uiaxes(blp);
        
        % Using 'plot_cross_section'-function
        plot_cross_section(cp, M, 0.1)

        % Changing its position
        cp.Position = [som/2 som/2 wop-som/2 hof/2-som*3.5];
        

%--------------------------------------------------------------------------


% Creating RIGHT BOTTOM panel

    % Creating uipanel in figure
    brp = uipanel(fig);

    % Changing its position
    brp.Position = [wop+som*2 som wop hop];
    if is_wall
        % Top right
        brp.Position = [wof*0.6 hof/2+som wof*0.4-som hof/2-som*2];
    end


    % Giving the panel a title

        % Inserting top panel text label
        blptl = uilabel(brp);

        % Changing size
        blptl.Position = [som som wop-som hop-som*1.5];
        if is_wall
            blptl.Position = [som som wof*0.4-som*3 hof/2-som*3.5];
        end

        % Specifying text (<br> works as extra linespacing)
        blptl.Text = html_text(M);

        % Setting interpreter
        blptl.Interpreter = 'html';

        % Changing vertical alignment from 'center' to 'top'
        blptl.VerticalAlignment = 'top';



% Ending main function
end

%--------------------------------------------------------------------------
%% Style of table rows


function rc = table_style(data)

% Allocating for size
rc = cell(size(data, 1),1); 

% Checking all given values and assign colors accordingly
for i = 1:size(data, 1)

    % If it is a minimum value
    if ~isempty(strfind(data{i,2}, 'min'))
        colors = {plot_colors().red, plot_colors().yellow, ...
            plot_colors().green};
    % Elseif it is a maximum value
    elseif  ~isempty(strfind(data{i,2}, 'max'))
        colors = {plot_colors().green, plot_colors().yellow, ...
            plot_colors().red};
    end

    % If the value is greater than the minimum value
    if cell2mat(data(i,4)) > cell2mat(data(i,5))
        rc{i} = colors{1};

    % If the to values are equal
    elseif cell2mat(data(i,4)) == cell2mat(data(i,5))
        rc{i} = colors{2};

    % If the value is less than the minimum value
    else
        rc{i} = colors{3};
    end

end

end



%--------------------------------------------------------------------------
% Rows of table data

function data = table_data(M, is_wall)

% Initializing counter
i = 1;

% If it is a column
if ~is_wall

% {5.3.1 - Structural models for overall analysis - p.57}
% (7)
data(i,:) = {'Minimum column length', 'L_min', ...
    '[m]', 3*max(M.h, M.w), M.L}; i = i+1;

% {9.5.1 - General - p.162}
% (1)
data(i,:) = {'Maximum height of cross section', 'h_max', ...
    '[m]', 4*min(M.h, M.w), max(M.h, M.w)}; i = i+1;

% {12.9.1 - Structural members - 198}
% (1)
data(i,:) = {'Minimum width of cross section', 'w_min', ...
'[mm]', 120, min(M.h, M.w)*10^3}; i = i+1;

% {9.5.2 - Longitudinal reinforcement - p.162}
% (1)
data(i,:) = {'Minimum diameter of longitudinal bars', 'ø_L_min', ...
    '[mm]', 8, min(M.o_L)*10^3}; i = i+1;

% (2)
data(i,:) = {'Minimum area of longitudinal bars', 'A_s_min', ...
    '[m^2]', min(0.1*M.N_Ed/M.f_yd, 0.002*M.A_c), M.A_s}; i = i+1;

% (3)
data(i,:) = {'Maximum area of longitudinal bars', 'A_s_max', ...
    '[m^2]', 0.08*M.A_c, M.A_s}; i = i+1;

% (4)
data(i,:) = {'Minimum number of longitudinal bars', 'n_L_min', ...
    '[-]', 4, sum(M.n_L)}; i = i+1;


% {9.5.3 - Transverse reinforcement - p.162}
%(1)
data(i,:) = {'Minimum diameter of transversal bars', 'o_T_min', ...
    '[mm]', max(6, 1/4*max(M.o_L)*10^3), M.o_T*10^3}; i = i+1;

%(3)
data(i,:) = {'Maximum spacing of transversal bars', 's_T_max', ...
    '[mm]', min([20*min(M.o_L), min(M.h, M.w), 400*10^(-3)]) *10^3, ...
    M.s_T*10^3}; i = i+1;

% {9.5.3 - Transverse reinforcement - p.162}
data(i,:) = {'Minimum geometric reinforcement ratio (RC)', 'rho_min', ...
    '[-]', 0.002, M.rho}; i = i+1;

% {BK - 3.3 Armeringsafstande - Fig. 3.13}
data(i,:) = {['Minimum distance between diameter circumference and ' ...
    'casting edge'], 'c_1_min', ...
    '[mm]', max([max(M.o_L)+5*10^(-3), M.c]) *10^3, (M.c + M.o_T)*10^3}; 
    i = i+1;

data(i,:) = {'Minimum distance between diameter circumferences', ...
    'a_min', '[mm]', max([max(M.o_L)*10^3, 32+5, 20]), ...
    min( (M.h-2*(M.c+M.o_T)-max(M.o_L))/(max(M.n_L)-1) - max(M.o_L), ...
         (M.w-2*(M.c+M.o_T)-max(M.o_L))/(max(M.n_L)-1) - max(M.o_L)) ...
     *10^3}; i = i+1;

end

%--------------------------------------------------------------------------

% If it is a wall
if is_wall

% {5.3.1 - Structural models for overall analysis - p.57}
% (7)
data(i,:) = {'Minimum wall height', 'L_min', ...
    '[mm]', 3*min(M.h, M.w)*10^3, M.L*10^3}; i = i+1;

% {9.6.1 - General - p.163}
% (1)
data(i,:) = {'Minimum width of cross section', 'w_min', ...
'[mm]', 4*min(M.h, M.w)*10^3, max(M.h, M.w)*10^3}; i = i+1;

% {12.9.1 - Structural members - 198}
% (1)
data(i,:) = {'Minimum thickness of wall', 'h_min', ...
'[mm]', 120, min(M.h, M.w)*10^3}; i = i+1;


% {9.6.2 - (Vertical) Longitudinal reinforcement - p.163}
% (1)
data(i,:) = {'Minimum area of longitudinal bars', 'A_L_min', ...
    '[m^2]', 0.002*M.A_c, M.A_s}; i = i+1;

% (2)
data(i,:) = {'Maximum area of longitudinal bars', 'A_L_max', ...
    '[m^2]', 0.04*M.A_c, M.A_s}; i = i+1;

% (3)
data(i,:) = {'Maximum spacing between longitudinal bars', 's_L_max', ...
    '[mm]', min(3*min(M.h, M.w)*10^3, 400), M.w/M.n_L(1)*10^3}; i = i+1;


% {9.6.3 - Horizontal reinforcement- p.163}
%(1)
data(i,:) = {'Minimum area of transversal bars', 'A_T_min', ...
    '[m^2]', max(0.25*M.A_s, 0.001*M.A_c), ...
    2*M.L/M.s_T * pi() * (M.o_T/2)^2}; i = i+1;

% (2)
data(i,:) = {'Maximum spacing between transversal bars', 's_T_max', ...
    '[mm]', 400, M.s_T *10^3}; i = i+1;

% {9.5.3 - Transverse reinforcement - p.162}
data(i,:) = {'Minimum geometric reinforcement ratio (RC)', 'rho_min', ...
    '[-]', 0.002, M.rho};

% {BK - 3.3 Armeringsafstande - Fig. 3.13}
data(i,:) = {['Minimum distance between diameter circumference and ' ...
    'wall edge'], 'c_1_min', '[mm]', ...
    max([max(M.o_L)+5*10^(-3), M.c]) *10^3, (M.c + M.o_T)*10^3}; i = i+1;

% {BK - 3.3 Armeringsafstande - Fig. 3.13}
% (Assuming that the tolerance addition is 5 mm)
data(i,:) = {'Minimum distance between diameter circumferences', ...
    'a_min', '[mm]', ...
    max([max(M.o_L)*10^3, 32+5, 20]), min(M.w/max(M.n_L)*10^3, ...
    M.s_T*10^3)};

end


% Ending function
end



%--------------------------------------------------------------------------
%% HTML text

function text = html_text(M)

% Font size and line height
[FS, LH] = deal('9', '250');

% Colors
[R, Y, G, B] = deal('#ff4747','#f9ff47', '#54ff71', '#000000');

% Possiable environmental class
EC = {'extra aggressive', 'aggressive', 'moderate', 'passive', 'ERROR'};

% Font size and line height (line spacing)
fslh = ['font-size:' FS 'px; line-height:' LH '%"'];

% Initializing and string counter for number of rows
i = 1; t = {};

% Header
t{i} = ['<p style="font-size:14px"><strong>' ...
        'Assessment of other parameters:' ...
        '</p></strong>']; i=i+1;

% Horizontal spacing and starting bulleted list
t{i} = '<p>&nbsp;</p><ul>'; i=i+1;

% Concrete class
j = M.f_ck*10^(-6);
if j<12 || j>90  c=R; elseif j<15 || j>80  c=Y; else c=B; end
t{i} = ['<li><p style="color:' c '; ' fslh '>' ...
        'The concrete class is: C' sprintf('%.f',M.f_ck*10^(-6)) ...
        '</p></li>']; i=i+1;

% Environment class
t{i} = ['<li><p style="color:' c '; ' fslh '>' ...
        'The environment class is according to:' ...
        '</p></li>']; i=i+1;

    % Starting integrated bulleted list
    t{i} = '<ul>'; i=i+1;

    % Cover layer
    j = M.c*10^(3) + 5;
    if j>=40 a=EC{1}; elseif j>=30 a=EC{2}; elseif j>=20 a=EC{3};
    elseif j>=10 a=EC{4}; else a=EC{5}; end
    if j<10 || j>60  c=R; elseif j==10 || j>50  c=Y; else c=B; end
    t{i} = ['<li><p style="color:' c '; ' fslh '>' ...
            'Minimum cover layer, ' a ' (c = ' ...
            sprintf('%.f',M.c*10^3+5) ' mm)'...
            '</p></li>']; i=i+1;

    % Concrete strength class
    j = M.f_ck*10^(-6);
    if j>=40 a=EC{1}; elseif j>=35 a=EC{2}; elseif j>=25 a=EC{3};
    elseif j>=12 a=EC{4}; elseif j<=12 a=EC{5}; else a=EC{5}; end
    if j<10 || j>60  c=R; elseif j==10 || j>50  c=Y; else c=B; end
    t{i} = ['<li><p style="color:' c '; ' fslh '>' ...
            'Minimum characteristic strength, ' a ' (f_ck = ' ...
            sprintf('%.f',M.f_ck*10^(-6)) ' MPa)'...
            '</p></li>']; i=i+1;

    % Ending integrated bulleted list
    t{i} = '</ul>'; i=i+1;

% Reinforcement steel (characteristic strain at maximum force)
j = M.varepsilon_uk*100;
if j<2.5 a='ERROR'; elseif j<5.0 a='A'; elseif j<7.5 a='B'; 
else a='C'; end
if j<2.5 || j>15  c=R; elseif j>10 c=Y; else c=B; end
t{i} = ['<li><p style="color:' c '; ' fslh '>' ...
        'The steel class of the reinforcement steel is: ' a ...
         ' (varepsilon_uk = ' sprintf('%.1f',j) ' %)'...
        '</p></li>']; i=i+1;

% Reinforcement steel (characteristic yield strength)
j = M.f_yk*10^(-6);
if j<=400 || j>600 a='WRONG'; else a='Correct'; end
if j<400 || j>600  c=R; elseif j<=400 || j>=600 c=Y; else c=B; end
t{i} = ['<li><p style="color:' c '; ' fslh '>' ...
        a ' char. yield strength of the reinforcement ' ...
        'steel (f_yk = ' sprintf('%.f',j) ' MPa)' ...
        '</p></li>']; i=i+1;

% Reinforcement steel (Young's modulus)
j = M.E_s*10^(-9);
if j~=200 a='WRONG'; else a='Correct'; end
if j~=200 c=R; else c=B; end
t{i} = ['<li><p style="color:' c '; ' fslh '>' ...
        a ' modolus of elesticity of reinforcement ' ...
        'steel (E = ' sprintf('%.f',j) ' GPa)' ...
        '</p></li>']; i=i+1;

% Control class
j = M.gamma_3;
if j==0.95 a='Intensified'; elseif j==1.00 a='Normal'; 
elseif j==1.10 a='Relaxed'; else a='ERROR'; end
if j<0.95 || j>1.10  c=R; else c=B; end
t{i} = ['<li><p style="color:' c '; ' fslh '>' ...
        a ' control class' ...
         ' (gamma_3 = ' sprintf('%.2f',j) ')' ...
        '</p></li>']; i=i+1;

% Concrete partial factor
j = M.gamma_c;
if j==1.45 a='Correct'; else a='Wrong'; end
if j~=1.45 c=R; else c=B; end
t{i} = ['<li><p style="color:' c '; ' fslh '>' ...
        a ' reinforced concrete partial factor, comp.' ...
         ' (gamma_c = ' sprintf('%.2f',j) ')' ...
        '</p></li>']; i=i+1;

% Steel partial factor
j = M.gamma_s;
if j==1.20 a='Correct'; else a='Wrong'; end
if j~=1.20 c=R; else c=B; end
t{i} = ['<li><p style="color:' c '; ' fslh '>' ...
        a ' steel partial factor' ...
         ' (gamma_s = ' sprintf('%.2f',j) ')' ...
        '</p></li>']; i=i+1;

% Ending bulleted list
t{i} = '</ul>';

% Joining the string array with no spacing
text = join(t,'');

end