function [fluxFlat]=alphaEqnNoUmFluxFlat(uFlat,V0,rhol,rhog,a)
    % Gives the flux for alphaEqn in mixture model with no Um term
    % u flatten version (u has plain data)
    %
    % [fluxFlat]=alphaEqnNoUmFlux(uFlat,V0,rhol,rhog,a)
    %
    % flux: flatten obtained flux
    % uFlat: independent unknown in flat form
    % V0: constant for relative velocity in flux calculation
    % rho: density of continue phase
    % rhog: density of dispersed phase
    % a: exponent for relative velocity in flux calculation

    % Deflatten
    u.internal=uFlat;
    u.left.setvalue=uFlat(1);
    u.right.setvalue=uFlat(end);

    % Takes the data size
    N=size(u.internal,1);

    Urlg=assign(constField(-V0,N),assign(u,constField(a,N),'^'),'*');
    % Um is calculated from an analytical expression since the complete solver is not available
    % Um is zero in this problem
    Um=constField(0,N);
    flux=assign(Urlg,assign(u,assign(constField(1,N),u,'-'),'*'),'*');
 
    % Flatten
    fluxFlat=flux.internal;


end