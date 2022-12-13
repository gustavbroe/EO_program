function print_figures(path_folder, name_file, number_figure)

%--------------------------------------------------------------------------
% Prints the figures which are currently open, excluding any 'Input check'
% windows.
%
% Input:    The path to the folder (string array)
%           The desired file name (string array)
%           The number of the figure (number array)
% Output:   None (other than the files printed)
%--------------------------------------------------------------------------

% Storing path to current folder and changing to results folder
oldFolder = cd(path_folder);

for i = 1:size(number_figure, 2)

    % String indicating figure to print
    sif = append('-f', string(number_figure(i)));

    % Printing the current figure
        % as PDF
        print(name_file(i), sif, '-dpdf', '-fillpage', '-r600')

        % as EPS (image/openGL format)
        print(name_file(i), sif, '-depsc', '-image', '-r600')

        % as EPS (vector format)
        print(append(name_file(i), '_v'), sif, '-depsc', ...
            '-vector', '-r600')

end

% Returning to the old current folder
cd(oldFolder)