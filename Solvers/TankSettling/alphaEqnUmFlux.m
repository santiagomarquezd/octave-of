function [flux]=alphaEqnUmFlux(u,V0,rhol,rhog,a)
    % Gives the flux for alphaEqn in mixture model with Um term
    %
    % [flux]=alphaEqnUmFlux(u,V0,rhol,rhog,a)
    %
    % flux: obtained flux
    % u: independent unknown
    % V0: constant for relative velocity in flux calculation
    % rho: density of continue phase
    % rhog: density of dispersed phase
    % a: exponent for relative velocity in flux calculation

    % Takes the data size
    N=size(u.internal,1);


    Urlg=assign(constField(-V0,N),assign(u,constField(a,N),'^'),'*');
    % Um is calculated from an analytical expression since the complete solver is not available
    % Um is zero in this problem
    Um=constField(0,N);
    flux=assign(assign(Urlg,assign(u,assign(constField(1,N),u,'-'),'*'),'*'),assign(Um,u,'*'),'+');
end