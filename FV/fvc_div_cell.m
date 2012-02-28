function [div]=fvc_div_cell(U, w, xC, xF, Sf, V)
    % Gives the explicit the divergence of a cell field
    %
    % [grad]=fvc_grad(U, w, xC, xF, Sf, V)
    %
    % div: the field divergence
    % U: given field
    % w: interpolation weights
    % xC: cell centers
    % xF: face centers
    % Sf: face areas
    % V: cell volumes
    
    % Face interpolation of U
    Uf=fvc_interpolate(U, w, xC, xF);
    
    % Divergence calculation
    div=(Uf(2:end).*Sf(2:end)-Uf(1:end-1).*Sf(1:end-1))./V;
end