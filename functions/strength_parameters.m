function P = strength_parameters(M)

%--------------------------------------------------------------------------
% This function determines the strength parameters for the rectangular 
% stress distribution, this is determined by the concrete stregnth class.
% 
% Input:    Information about the member (structure array)
% Output:   The four strength parameters in a struc array
%--------------------------------------------------------------------------


% {"Ultimate limit state calculation" - 4.3 - p.109}
if (M.f_ck >= 12*10^(6)) && (M.f_ck <= 50*10^(6))

    % Factor determining the effective height of the compression zone
    P.lambda = 0.8;                                                   %[-]

    % Factor defining the effective strength
    P.eta = 1.0;                                                      %[-]

    % Compressive strain in the concrete at the peak stress, f_c
    P.varepsilon_c3 = 0.175 /100;                                     %[%]

    % Ultimate compressive strain in the concrete
    P.varepsilon_cu3 = 0.35 /100;                                     %[%]

elseif (M.f_ck > 50*10^(6)) && (M.f_ck <= 90*10^(6))

    P.lambda = 0.8 - (M.f_ck*10^(-6) - 50)/400;                       %[-]
    P.eta = 1.0 - (M.f_ck*10^(-6) - 50)/200;                          %[-]
    P.varepsilon_c3 = 0.175 + 0.055 * (M.f_ck*10^(-6) - 50)/40;       %[%]
    P.varepsilon_cu3 = 0.26 + 3.5 * ((90 - M.f_ck*10^(-6))/100)^4;    %[%]

end


end