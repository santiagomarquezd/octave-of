function [A, RHS]=fvm_ddt(factor,factor0,xn,V,dt,method)
    % Gives the temporal derivative matrix and corresponding RHS
    % d/dt x^(n+1)-x^(n)*V
    % 
    % fvm_ddt(V,method)
    %
    % factor: slacar factor in term at time n+1
    % factor0: slacar factor in term at time n
    % xn: column vector with x n state
    % V: column vector with volumes of each cell
    % dt: time-step
    % method: 1. Backward Euler, 2. Crank-Nicholson

    if (method==1)
      % Allocate matrix and RHS
      A=diag(ones(size(V)).*V/dt.*factor.internal);
      RHS=xn.internal.*V/dt.*factor0.internal;
    elseif (method==2)
    
    end
end
