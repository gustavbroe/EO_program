function [price_tot, price_c, price_s, o_LF] = calc_Price(info_lev, M, concrete, steel)

%--------------------------------------------------------------------------
% This function calcualtes the price of materials used to construct a given
% wall or column. Valid with the following manufacturers and assumptions.
%
% Concrete from UNICON (in accordance with DS/EN 206 DK NA)
%  - Every truckload has to be pumped.
% 
% Steel from Celsa Steel Service
%  - Prices are calculated for small customers in November 2022.
%  - Different calculation procedure apply for hand-welded steel and nets.
%  - Total weight and longitudinal diameter determines price.
%  - Requirement about the maximum difference of 4 mm in diameter. (neglected)
%  - Transportation in Zone 1 (approx. 39 from Ølstykke to DTU) 15-20 tons.
%    (The total weight of steel in almost optimal walls and columns)
%
% Input:    Information about the member (structure array)  
% Output:   Price of products [DKK]
%--------------------------------------------------------------------------

% Import pricing information
PI = pricing_info(info_lev, steel);

% Checking if the member is a wall
if M.w > 4*M.h
    is_wall = logical(1);
else
    is_wall = logical(0);
end

% Rounding up to the nearest valid rebar diameter
% (Very inefficient, but it gets the job done)
if ~is_wall
    if max(M.o_L)*10^(3) <= 6
        asp = PI.s.hwp.o6; 
        o_LF = 6;
        M = def_c(M.f_ck*10^6, M.h*10^3, M.w*10^3, ...
            repmat(6,1,size(M.o_L,2)), M.n_L, M.d_L(2:end-1));
    elseif max(M.o_L)*10^(3) <= 8
        asp = PI.s.hwp.o8;
        o_LF = 8;
        M = def_c(M.f_ck*10^6, M.h*10^3, M.w*10^3, ...
            repmat(8,1,size(M.o_L,2)), M.n_L, M.d_L(2:end-1));
    elseif max(M.o_L)*10^(3) <= 10
        asp = PI.s.hwp.o10;
        o_LF = 10;
        M = def_c(M.f_ck*10^6, M.h*10^3, M.w*10^3, ...
            repmat(10,1,size(M.o_L,2)), M.n_L, M.d_L(2:end-1));
    elseif max(M.o_L)*10^(3) <= 12
        asp = PI.s.hwp.o12;
        o_LF = 12;
        M = def_c(M.f_ck*10^6, M.h*10^3, M.w*10^3, ...
            repmat(12,1,size(M.o_L,2)), M.n_L, M.d_L(2:end-1));
    elseif max(M.o_L)*10^(3) <= 14
        asp = PI.s.hwp.o14;
        o_LF = 14;
        M = def_c(M.f_ck*10^6, M.h*10^3, M.w*10^3, ...
            repmat(14,1,size(M.o_L,2)), M.n_L, M.d_L(2:end-1));
    elseif max(M.o_L)*10^(3) <= 16
        asp = PI.s.hwp.o16;
        o_LF = 16;
        M = def_c(M.f_ck*10^6, M.h*10^3, M.w*10^3, ...
            repmat(16,1,size(M.o_L,2)), M.n_L, M.d_L(2:end-1));
    elseif max(M.o_L)*10^(3) <= 20
        asp = PI.s.hwp.o20;
        o_LF = 20;
        M = def_c(M.f_ck*10^6, M.h*10^3, M.w*10^3, ...
            repmat(20,1,size(M.o_L,2)), M.n_L, M.d_L(2:end-1));
    else
        asp = NaN;
    end

elseif is_wall
    if max(M.o_L)*10^(3) <= 6
        asp = PI.s.cn.o6;
        o_LF = 6;
        M = def_w(M.f_ck*10^6, M.h*10^3, M.w*10^3, ...
            repmat(6,1,size(M.o_L,2)), M.n_L, M.d_L(2:end-1));
    elseif max(M.o_L)*10^(3) <= 8
        asp = PI.s.cn.o8;
        o_LF = 8;
        M = def_w(M.f_ck*10^6, M.h*10^3, M.w*10^3, ...
            repmat(8,1,size(M.o_L,2)), M.n_L, M.d_L(2:end-1));
    elseif max(M.o_L)*10^(3) <= 10
        asp = PI.s.cn.o10;
        o_LF = 10;
        M = def_w(M.f_ck*10^6, M.h*10^3, M.w*10^3, ...
            repmat(10,1,size(M.o_L,2)), M.n_L, M.d_L(2:end-1));
    elseif max(M.o_L)*10^(3) <= 12
        asp = PI.s.cn.o12;
        o_LF = 12;
        M = def_w(M.f_ck*10^6, M.h*10^3, M.w*10^3, ...
            repmat(12,1,size(M.o_L,2)), M.n_L, M.d_L(2:end-1));
    elseif max(M.o_L)*10^(3) <= 16
        asp = PI.s.cn.o16;
        o_LF = 16;
        M = def_w(M.f_ck*10^6, M.h*10^3, M.w*10^3, ...
            repmat(16,1,size(M.o_L,2)), M.n_L, M.d_L(2:end-1));
    else
        asp = NaN;
    end
end


% Volumens

    % Steel (both vertical and horizontal)
    V.s = (M.A_s*M.L) + M.V_sT;     %[m^3]

    % Concrete
    V.c = (M.A_c*M.L) - V.s;        %[m^3]

% Calculating prices

    % Concrete
    price_c = V.c * PI.c.(concrete(1:end-1));        %[DKK]

    % Steel
    price_s = V.s * asp;                    %[DKK]

    % Total
    price_tot = price_c + price_s;          %[DKK]

   
%--------------------------------------------------------------------------
% Printing informataion if requested

% Basic
if info_lev == 1

    fprintf('price = %.1f [DKK] | ', price_tot)

end

% Extended
if info_lev == 2

    fprintf('price(tot, c, s) = (%.1f, %.1f, %.1f) [DKK] | ' ...
        ,price_tot, price_c, price_s);

end



% Ending the function
end



%--------------------------------------------------------------------------
% Struct array with pricing information
function PI = pricing_info(info_lev, steel)

% Loading matertial 
MI = material_info(0);

% PRICING FOR CONCRETE
%
% Special thanks to Catrine from UNICON northern zealand sales team for 
% answering additional and clarifying questions.
% 
% SOURCES:
%   Manufacturer: UNICON
%   (https://www.unicon.dk/media/1852/pricelist-november-2022-web-lr.pdf)


    % Concrete products

    % Mandotory surcharge (Maybe, Extended Control Class: 25 [DKK/m^3])
    % (Price of transportation is included in the price if the order is
    % greater 7 cubic meters and no amount limit per truck is requested)
    
    %      CO2 Surcharge    Energy Supplement    Environmental Supplement
    extra =     47       +        90          +             70          ...
        +       32               +           114;       %[DKK/m^3]
    %  Environmental Fee (Pump)      Concrete pumping (31m)


    % Price for 1 cubic meter of UNI-GREEN / Passive Slump concrete (C16)
    PI.c.c16 = 1280 + extra;        %[DKK/m^3]

    % Price for 1 cubic meter of UNI-GREEN / Passive Slump concrete (C20)
    PI.c.c20 = 1315 + extra;        %[DKK/m^3]

    % Price for 1 cubic meter of UNI-GREEN / Passive Slump concrete (C25)
    PI.c.c25 = 1351 + extra;        %[DKK/m^3]

    % Price for 1 cubic meter of UNI-GREEN / Moderate Slump concrete (C30)
    PI.c.c30 = 1388 + extra;        %[DKK/m^3]

    % Price for 1 cubic meter of UNI-GREEN / Aggressiv Slump concrete (C35)
    PI.c.c35 = 1425 + extra;        %[DKK/m^3]

    % Price for 1 cubic meter of UNI-GREEN / Extra aggresive Slump concrete (C40)
    PI.c.c40 = 1891 + extra;        %[DKK/m^3]



% PRICING FOR STEEL
%
% Special thanks to Claus Hammer, Project konsultant from Celsa Steel 
% Service for taken his time to e-mail a pricing list.
% 
% SOURCES:
%   Manufacturer: Celsa Steel Service A/S
%   (Sales department)


    % Steel products

    % Per order paymet
    PI.s.pop = 2000         +        1500;
    %   Invoicing surcharge     Bending lists

    % Hand-welded products (Column, Beams, etc.)
    % (Price per ton changes with longitudinal diameter)
    % [550B quality, in accordance with EN/ISO 17660-2]

        % Transportation
        PI.s.hwp.trans = 225;

        % Different longitudinal diameters
        PI.s.hwp.o6 = (PI.s.hwp.trans + 27830) * MI.(steel).rho*10^(-3);
        PI.s.hwp.o8 = (PI.s.hwp.trans + 25830) * MI.(steel).rho*10^(-3);
        PI.s.hwp.o10 = (PI.s.hwp.trans + 22830) * MI.(steel).rho*10^(-3);
        PI.s.hwp.o12 = (PI.s.hwp.trans + 19830) * MI.(steel).rho*10^(-3);
        PI.s.hwp.o14 = (PI.s.hwp.trans + 17830) * MI.(steel).rho*10^(-3);
        PI.s.hwp.o16 = (PI.s.hwp.trans + 14830) * MI.(steel).rho*10^(-3);
        PI.s.hwp.o20 = (PI.s.hwp.trans + 13830) * MI.(steel).rho*10^(-3);


    % Mesh reinforcement (Celsanets, with bars per 150 mm)
    % (Price per ton changes with longitudinal diameter)
    % [550B quality, in accordance with DS/EN 10080 and DS/INF 165]

        % Transportation
        PI.s.cn.trans = 150;

        % Different longitudinal diameters
        PI.s.cn.o6 = (PI.s.cn.trans + 11325) * MI.(steel).rho*10^(-3);
        PI.s.cn.o8 = (PI.s.cn.trans + 10950) * MI.(steel).rho*10^(-3);
        PI.s.cn.o10 = (PI.s.cn.trans + 10750) * MI.(steel).rho*10^(-3);
        PI.s.cn.o12 = (PI.s.cn.trans + 10675) * MI.(steel).rho*10^(-3);
        PI.s.cn.o16 = (PI.s.cn.trans + 11475) * MI.(steel).rho*10^(-3);


% Printing information
if info_lev == 2

    fprintf(['The mean value of building material prices \n' ...
        '[DKK/m^3] \n\n'])

    fprintf('Slump concrete, total: %.1f \n', ...
        mean(cell2mat(struct2cell(PI.c))))
    fprintf('Hand-welded steel, column: %.1f \n', ...
        mean(cell2mat(struct2cell(PI.s.hwp))))
    fprintf('Mesh reinforcement, celsanets: %.1f \n', ...
        mean(cell2mat(struct2cell(PI.s.cn))))

    % Footer
    fprintf(['-----------------------------------------------' ...
        '-------------------------------\n\n'])

end

% Ending sub-function
end