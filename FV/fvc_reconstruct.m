function [a]=fvc_reconstruct(af,Sf)
    % Reconstructs surface field to cell field
    % fvm_div_flux_cell(phi,Sf,method)
    %
    % [a]=fvc_reconstruct(af,Sf)
    %
    % a: cell based reconstructed field
    % af: surface field & Sf
    % Sf: face areas
    
    
    a=((-Sf(1:end-1)).*(-Sf(1:end-1))./abs(Sf(1:end-1))+Sf(2:end).*Sf(2:end)./abs(Sf(2:end))).^(-1).* ...
    (af(1:end-1)+af(2:end));

end
