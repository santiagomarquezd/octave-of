function A = fvm_Sp(Spfield,V)
    % Gives the matrix
    % fvm_Sp(Spfield,V)
    %
    % Spfield:  
    % V: Volume
    %
    % fvm.diag() += mesh.V()*sp.field();
    
    A=diag(Spfield.*V);

end
