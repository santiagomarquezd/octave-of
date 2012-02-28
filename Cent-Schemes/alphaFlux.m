function [flux]=alphaFlux(Vm,alphag,Vr,rhol,rhog,model)
    % Gives the flux for alphag for several mixtureModels
    %
    % [flux]=alphaFlux(alphag,Vm,cp,Vr,model)
    %
    % flux: obtained flux
    % Vm: mixture velocity
    % alphag: alphag field
    % rhol: liquid density
    % rhog: gas density
    % model: 1. No Vr; 2: with Vr
    
    if (model==1)
	  flux=Vm.*alphag;    
    elseif (model==2)
	  disp('Not implemented yet')
    end
end