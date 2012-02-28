function [A] = fvm_Sp(SpField, V)
    % Gives the matrix for an implicit source operator
    %
    % [A] = fvm_Sp(Spfield, V)
    %
    % SpField: source field  
    % V: cell volumes
    %
    
    A=diag(Spfield.*V);

end
