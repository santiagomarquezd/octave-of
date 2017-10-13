function [fluxFlat] = alphaEqnUmFluxFlat(uFlat, Um, V0, rhol, rhog, a)
    % Gives the flux for alphaEqn in mixture model with Um term
    % u flatten version (u has plain data)
    %
    % [fluxFlat] = alphaEqnUmFluxFlat(u, V0, Um, rhol, rhog, a)
    %
    % fluxFlat: obtained flux
    % u: independent unknown as a vector (no BC's)
    % V0: constant for relative velocity in flux calculation
    % Um: center of volume velocity
    % rho: density of continue phase
    % rhog: density of dispersed phase
    % a: exponent for relative velocity in flux calculation

    % Deflatten
    u.internal= uFlat;
    u.left.setvalue = uFlat(1);
    u.right.setvalue = uFlat(end);

    % Takes the data size
    N = size(u.internal, 1);

    % Urlg = -V0*u;
    Urlg = assign(constField(-V0, N), assign(u, constField(a, N), '^'), '*');
    % flux = (1 - u)*u*Urlg + Um*u
    flux = assign(assign(Urlg, assign(u, assign(constField(1, N), u, '-'),...
           '*'), '*'), assign(Um, u, '*'), '+');

    % Flatten
    fluxFlat = flux.internal;

end
