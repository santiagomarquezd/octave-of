function [flux]=alphaEqnNoVmFlux(u,V0,rhol,rhog,a)
    % Gives the flux for alphaEqn in mixture model with no Vm term
    %
    % [flux]=alphaEqnNoVmFlux(u,V0,rhol,rhog,a)
    %
    % flux: obtained flux
    % u: independent unknown
    % V0: constant for relative velocity in flux calculation
    % rho: density of continue phase
    % rhog: density of dispersed phase
    % a: exponent for relative velocity in flux calculation

    % Takes the data size
    N=size(u.internal,1);

    rhom=assign(assign(constField(rhog,N),u,'*'),assign(assign(constField(1,N),u,'-'),constField(rhol,N),'*'),'+');
    cp=assign(assign(constField(rhog,N),u,'*'),rhom,'/');
    flux=assign(assign(assign(constField(V0,N),assign(constField(1,N),u,'-'),'*'),u,'*'),assign(constField(1,N),cp,'-'),'*');
 
end