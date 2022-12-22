function MI = material_info(info_lev)

%--------------------------------------------------------------------------
% Outputs all relevant information about the products used in this thesis.
% This includes accurate value for climate information and density.
%
%
% Input:    The requested amount of printed information
% Output:   Material information about products (structure array)
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% CLIMATE INFORMATION
%
% UNITS: All units are convertet to [kg CO2-ep / m^3]
%
% Densities originate from the material data section of the EPD.
%
% PROCESS PHASE EXPLANATION
%
%   Parentheses around a process stage means that it is:
%       - MND = Modole Not Declared
%
%   Square parentheses mean that the process is:
%       - MNR = Module Not Relevant


    % Concrete products


    % GWP from 1 cubic meter of C16 with FUTURECEM cem. (Product EPD) + Density
    % (EPD Danmark x Unicon A/S)
    % {https://www.epddanmark.dk/media/jufkifzq/md-22092-da.pdf}
    MI.c16F.GWP = 1.60e2 + 3.28 + 0 + (-1.12e1) +    0   + 1.20e1 + 6.50 + 6.77 + 4.93; %(-4.56);
    %            A1-A3     A4   (A5)    B1        [B2-B7]    C1      C2     C3     C4       D
    MI.c16F.rho = 2.22e3;   %[kg/m^3]

    % GWP from 1 cubic meter of C16 with RAPID cem. (Product EPD) + Density
    % (EPD Danmark x Unicon A/S)
    % {https://www.epddanmark.dk/media/t1ge0gnr/md-22090-da.pdf}
    MI.c16R.GWP = 1.88e2 + 3.28 + 0 + (-1.38e1) +    0   + 1.20e1 + 6.50 + 6.67 + 4.93; %(-4.56);
    %            A1-A3     A4   (A5)    B1        [B2-B7]    C1      C2     C3     C4       D
    MI.c16R.rho = 2.22e3;   %[kg/m^3]


    % GWP from 1 cubic meter of C20 with FUTURECEM cem. (Product EPD) + Density
    % (EPD Danmark x Unicon A/S)
    % {https://www.epddanmark.dk/media/ksgclghy/md-22093-da.pdf}
    MI.c20F.GWP = 1.78e2 + 3.26 + 0 + (-1.28e1) +    0   + 1.19e1 + 6.45 + 6.62 + 4.89; %(-4.53);
    %            A1-A3     A4   (A5)    B1        [B2-B7]    C1      C2     C3     C4       D
    MI.c20F.rho = 2.22e3;   %[kg/m^3]

    % GWP from 1 cubic meter of C20 with RAPID cem. (Product EPD) + Density
    % (EPD Danmark x Unicon A/S)
    % {https://www.epddanmark.dk/media/unajvaoe/md-22091-da.pdf}
    MI.c20R.GWP = 2.10e2 + 3.29 + 0 + (-1.56e1) +    0   + 1.20e1 + 6.51 + 6.68 + 4.94; %(-4.57);
    %            A1-A3     A4   (A5)    B1        [B2-B7]    C1      C2     C3     C4       D
    MI.c20R.rho = 2.22e3;   %[kg/m^3]


    % GWP from 1 cubic meter of C25 with FUTURECEM cem. (Product EPD) + Density
    % (EPD Danmark x Unicon A/S)
    % {https://www.epddanmark.dk/media/dovcgbn4/md-21108-da.pdf}
    MI.c25F.GWP = 1.98e2 + 3.29 + 0 + (-9.51) +    0   + 1.20e1 + 6.51 + 6.68 + 4.93; %(-4.57);
    %            A1-A3     A4   (A5)    B1        [B2-B7]    C1      C2     C3     C4       D
    MI.c25F.rho = 2.22e3;   %[kg/m^3]

    % GWP from 1 cubic meter of C25 with RAPID cem. (Product EPD) + Density
    % (EPD Danmark x Unicon A/S)
    % {https://www.epddanmark.dk/media/bzln31ik/md-21021-da_unicon.pdf}
    MI.c25R.GWP = 2.31e2 + 3.30 + 0 + (-1.14e1) +    0   + 1.21e1 + 6.54 + 6.71 + 4.96; %(-4.59);
    %            A1-A3     A4   (A5)    B1        [B2-B7]    C1      C2     C3     C4       D
    MI.c25R.rho = 2.23e3;   %[kg/m^3]


    % GWP from 1 cubic meter of C30 with FUTURECEM cem. (Product EPD) + Density
    % (EPD Danmark x Unicon A/S)
    % {https://www.epddanmark.dk/media/dkclcn3n/md-21109-da.pdf}
    MI.c30F.GWP = 2.27e2 + 3.32 + 0 + (-1.10e1) +    0   + 1.22e1 + 6.58 + 6.75 + 4.99; %(-4.62);
    %            A1-A3     A4   (A5)    B1        [B2-B7]    C1      C2     C3     C4       D
    MI.c30F.rho = 2.25e3;   %[kg/m^3]

    % GWP from 1 cubic meter of C30 with RAPID cem. (Product EPD) + Density
    % (EPD Danmark x Unicon A/S)
    % {https://www.epddanmark.dk/media/kfceit4y/md-21026-da_unicon.pdf}
    MI.c30R.GWP = 2.89e2 + 3.33 + 0 + (-1.38e1) +    0   + 1.22e1 + 6.60 + 6.77 + 5.00; %(-4.63);
    %            A1-A3     A4   (A5)    B1        [B2-B7]    C1      C2     C3     C4       D
    MI.c30R.rho = 2.25e3;   %[kg/m^3]


    % GWP from 1 cubic meter of C35 with  LAVALKALI SULFATBESTANDIG cem. (Product EPD) + Density
    % (EPD Danmark x Unicon A/S)
    % {https://www.epddanmark.dk/media/dkclcn3n/md-21109-da.pdf}
    MI.c35L.GWP = 4.14e2 + 3.41 + 0 + (-5.43) +    0   + 1.25e1 + 6.76 + 6.94 + 5.14; %(-4.75);
    %            A1-A3     A4   (A5)    B1        [B2-B7]    C1      C2     C3     C4       D
    MI.c35L.rho = 2.31e3;   %[kg/m^3]

    % GWP from 1 cubic meter of C35 with RAPID cem. (Product EPD) + Density
    % (EPD Danmark x Unicon A/S)
    % {https://www.epddanmark.dk/media/kfceit4y/md-21026-da_unicon.pdf}
    MI.c35R.GWP = 3.66e2 + 3.37 + 0 + (-5.16) +    0   + 1.23e1 + 6.67 + 6.85 + 5.06; %(-4.68);
    %            A1-A3     A4   (A5)    B1        [B2-B7]    C1      C2     C3     C4       D
    MI.c35R.rho = 2.28e3;   %[kg/m^3]


    % GWP from 1 cubic meter of C40 with  LAVALKALI SULFATBESTANDIG cem. (Product EPD) + Density
    % (EPD Danmark x Unicon A/S)
    % {https://www.epddanmark.dk/media/dkclcn3n/md-21109-da.pdf}
    MI.c40L.GWP = 4.62e2 + 3.41 + 0 + (-4.09) +    0   + 1.25e1 + 6.76 + 6.93 + 5.12; %(-4.74);
    %            A1-A3     A4   (A5)    B1        [B2-B7]    C1      C2     C3     C4       D
    MI.c40L.rho = 2.30e3;   %[kg/m^3]

    % GWP from 1 cubic meter of C40 with RAPID cem. (Product EPD) + Density
    % (EPD Danmark x Unicon A/S)
    % {https://www.epddanmark.dk/media/kfceit4y/md-21026-da_unicon.pdf}
    MI.c40R.GWP = 4.32e2 + 3.39 + 0 + (-1.26e-1) +    0   + 1.24e1 + 6.71 + 6.88 + 5.09; %(-4.71);
    %            A1-A3     A4   (A5)    B1        [B2-B7]    C1      C2     C3     C4       D
    MI.c40R.rho = 2.29e3;   %[kg/m^3]



    % Steel products

    % GWP from 1 tonne of reinforcement products (Product EPD) + Density
    % (EPD International AB x Celsa Steel Service A/S)
    % {https://api.environdec.com/api/v1/EPDLibrary/Files/b8f8bd47-8c39-480b-3a09-08d98fadb225/Data}
    MI.rp.GWP = 422 + 7.21 + 0.341 + 6.85 + 0; %+ 120;
    %          A1-A3   A4      C1    C2  C3-C4  D
    MI.rp.rho = 7700;   %[kg/m^3]
        
        % Converting unit from [kg CO2-eq pr. ton] to [kg CO2-eq pr. m^3]
        MI.rp.GWP = MI.rp.GWP * 1/1000 * MI.rp.rho;
    
    
    % GWP from 1 kg of steel rebar (Product EPD) + Density
    % (The Norwegian EPD Foundation x E.A. Smith AS)
    % {https://www.epd-norge.no/getfile.php/1313430-1642584515/EPDer/Byggevarer/St%C3%A5lkonstruksjoner/NEPD-2193-988_Steel-rebar.pdf}
    MI.sr.GWP = 3.93e-1 + 3.55e-3 + 0 +   0   + 6.09e-4 + 7.30e-3 + 0 + 1.28e-3; %+ 2.17e-1;
    %            A1-A3      A4    (A5) (B1-B7)    C1        C2     C3     C4         D
    MI.sr.rho = 7850;   %[kg/m^3]   {assumed, but likely accurate}
    
        % Converting unit from [kg CO2-eq pr. ton] to [kg CO2-eq pr. m^3]
        MI.sr.GWP = MI.sr.GWP * MI.sr.rho;
    
    
    % GWP from 1 metric tonne of Steel reinforcement products [Danish sites] (Product EPD) + Density
    % (The International EPD® System x Lemvigh-Müller)
    % {https://www.lemu.dk/da/service/information-og-dokumentation/miljoevaredeklarationer-epd}
    MI.srp.GWP = 7.73e2 + 2.48e1;
    %            A1-A3      A4   
    MI.srp.rho = 7850;   %[kg/m^3]   {assumed, but likely accurate}
    
        % Converting unit from [kg CO2-eq pr. ton] to [kg CO2-eq pr. m^3]
        MI.srp.GWP = MI.srp.GWP * 1/1000 * MI.srp.rho;



%--------------------------------------------------------------------------
% Printing informataion if requested

% Climate information
if info_lev == 1
    
    % Header, including unit
    fprintf(['Global warming potential (GWP) from building ' ...
        'product EPDs: \n[kg CO2-ep pr. m^3] \n\n'])

    % Products
    fprintf('%.2f / C16 with FUTURECEM cement (c16F) \n', MI.c16F.GWP)
    fprintf('%.2f / C16 with RAPID cement (c16R) \n', MI.c16R.GWP)
    fprintf('%.2f / C20 with FUTURECEM cement (c20F) \n', MI.c20F.GWP)
    fprintf('%.2f / C20 with RAPID cement (c20R) \n', MI.c20R.GWP)
    fprintf('%.2f / C25 with FUTURECEM cement (c25F) \n', MI.c25F.GWP)
    fprintf('%.2f / C25 with RAPID cement (c25R) \n', MI.c25R.GWP)
    fprintf('%.2f / C30 with FUTURECEM cement (c30F) \n', MI.c30F.GWP)
    fprintf('%.2f / C30 with RAPID cement (c30R) \n', MI.c30R.GWP)
    fprintf(['%.2f / C35 with LAVALKALI SULFATBESTANDIG cement (c35L) ' ...
        '\n'], MI.c35L.GWP)
    fprintf('%.2f / C35 with RAPID cement (c35R) \n', MI.c35R.GWP)
    fprintf(['%.2f / C40 with LAVALKALI SULFATBESTANDIG cement (c40L) ' ...
        '\n'], MI.c40L.GWP)
    fprintf('%.2f / C40 with RAPID cement (c40R) \n', MI.c40R.GWP)
    fprintf('%.2f / Reinforcement products (rp) \n', MI.rp.GWP)
    fprintf('%.2f / Steel reinforcement products (srp) \n', MI.srp.GWP)
    fprintf('%.2f / Steel rebar (sr) \n', MI.sr.GWP)

    % Footer
    fprintf(['-----------------------------------------------' ...
        '-------------------------------\n'])

end


% Ending function
end