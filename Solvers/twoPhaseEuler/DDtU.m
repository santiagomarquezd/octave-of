% {
%     DDtUa =
%         fvc::ddt(Ua)
%       + fvc::div(phia, Ua)
%       - fvc::div(phia)*Ua;
       DDtUa = ...
          fvc_ddt(Ua, Ua0, dt)...
        + fvc_div_cell(Ua, w, xC, xF, Sf, V) ...
        - fvc_div_face(phia, V).*Ua.internal;
% 
%     DDtUb =
%         fvc::ddt(Ub)
%       + fvc::div(phib, Ub)
%       - fvc::div(phib)*Ub;
% }

       DDtUb = ...
          fvc_ddt(Ub, Ub0, dt)...
        + fvc_div_cell(Ub, w, xC, xF, Sf, V) ...
        - fvc_div_face(phib, V).*Ub.internal;