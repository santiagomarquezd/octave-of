function [Kc]=fvc_ddtPhiCoeff(U, phi, w, xC, xF, Sf)
    % Gives the coefficient for fvcDdtPhi correction calculation
    %
    % [Kc]=fvc_ddtPhiCoeff(U, phi, xC, xF, Sf)
    %
    % U: velocity field
    % phi: volumetric flux at faces
    % w: interpolation weights
    % xC: cell centres positions
    % xF: face centres positions
    % Sf: face areas

 
    Kc=1-min(abs(phi-(Sf.*fvc_interpolate(U,w,xC,xF)))./(abs(phi) + 1E-37),1);
 
    if (U.left.type=='V')
      Kc(1)=0;
    elseif (U.right.type=='V')
      Kc(end)=0;
    end
      
end