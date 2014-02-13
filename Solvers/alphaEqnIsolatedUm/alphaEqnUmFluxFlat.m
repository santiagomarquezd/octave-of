function [fluxFlat]=alphaEqnUmFluxFlat(uFlat,Um,V0,rhol,rhog,a)
    % Gives the flux for alphaEqn in mixture model with Um term
    % u flatten version (u has plain data)
    %
    % [fluxFlat]=alphaEqnUmFluxFlat(u,V0,Um,rhol,rhog,a)
    %
    % fluxFlat: obtained flux
    % u: independent unknown as a vector (no BC's)
    % V0: constant for relative velocity in flux calculation
    % Um: center of volume velocity
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
    % DEPRECATED, full Um value is passed
    % Um is calculated from an analytical expression since the complete solver is not available
    % Um is zero in this problem
    %Um=constField(0,N);
    flux=assign(assign(Urlg,assign(u,assign(constField(1,N),u,'-'),'*'),'*'),assign(Um,u,'*'),'+');

    % Flatten
    fluxFlat=flux.internal;

end