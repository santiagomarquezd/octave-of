function [grad]=fvc_snGrad(a,xC,xF)
    % Finds the face gradient dotted by Sf vector
    % (gradient at faces for 1D)
    %
    % [grad]=fvc_snGrad(a,xC)
    %
    % grad: face gradient
    % a: cell based field
    % xC: cell centres
    % xF: face centres
    
    % Memory allocation
    grad=0*xF;
    
    % Calculation at internal faces 
    grad(2:end-1)=(a.internal(2:end)-a.internal(1:end-1))./(xC(2:end)-xC(1:end-1));
    
    % Calculation at boundaries
    if (a.left.type=='V')
      grad(1)=(a.internal(1)-a.left.setvalue)/(xC(1)-xF(1));
    elseif (a.left.type=='G')
      grad(1)=a.left.setvalue;
    end

    if (a.right.type=='V')
      grad(end)=(a.right.setvalue-a.internal(end))/(xF(end)-xC(end));
    elseif (a.right.type=='G')
      grad(end)=a.right.setvalue;
    end

    
    

end
