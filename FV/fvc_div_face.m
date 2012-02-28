function [div]=fvc_div_face(Uf, V)
    % Gives the explicit divergence of a face field
    %
    % [div]=fvc_div_face(U, w, xC, xF, Sf, V)
    %
    % div: the field divergence
    % Uf: given face field
    % V: cell volumes
    
    % Divergence calculation
    div=(Uf(2:end)-Uf(1:end-1))./V;
    V;
end