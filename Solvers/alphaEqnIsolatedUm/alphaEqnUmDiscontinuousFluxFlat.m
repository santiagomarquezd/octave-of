function [fluxFlat]=alphaEqnUmDiscontinuousFluxFlat(uFlat,Um,V0,rhol,rhog,a)
    % Gives the flux for alphaEqn in mixture model with Um term
    % u flatten version (u has plain data). Flux is discontinuous in u
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
    
    % Selection value
    sel=0.25;
    
    % Allocation of selector
    selection=u;
    % Values
    selection.internal=u.internal>sel;
    selection.left.setvalue=u.left.setvalue>sel;
    selection.right.setvalue=u.right.setvalue>sel;
   
    % Relative velocity calculation
    Urlg=assign(assign(assign(constField(-V0,N),assign(u,constField(a,N),'^'),'*'),selection,'*'), 
    assign(assign(constField(-V0*0.5,N),assign(u,constField(a,N),'^'),'*'),assign(constField(1,N),selection,'-'),'*'),'+');   
    
    %keyboard; pause;
    
    flux=assign(assign(Urlg,assign(u,assign(constField(1,N),u,'-'),'*'),'*'),assign(Um,u,'*'),'+');

    % Flatten
    fluxFlat=flux.internal;

end