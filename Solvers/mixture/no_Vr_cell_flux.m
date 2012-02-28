function [flux1,flux2]=no_Vr_cell_flux(U,alphag,rhom)
    % Gives the fluxes by cells for no Vr mixture model
    % (also calculates flux at BC's)
    %
    % [flux1,flux2]=no_Vr_cell_flux(U,alphag,rhom)
    %
    % flux1,flux2: cell fluxes
    % U: cell velocity
    % alphag: cell void fraction
    % rhom: cell density
    
    % Allocation
    flux1=U;
    flux2=U;
    
    % Internal field calculation
    flux1.internal=rhom.internal.*U.internal.^2;
    flux2.internal=alphag.internal.*U.internal;
    
    % BC's
    flux1.left.type='V';
    flux1.left.value=flux1.left.setvalue=rhom.left.setvalue*U.left.setvalue^2;
    
    flux1.right.type='V';
    flux1.right.value=flux1.right.setvalue=rhom.right.setvalue*U.right.setvalue^2;
    
    flux2.left.type='V';
    flux2.left.value=flux2.left.setvalue=alphag.left.setvalue*U.left.setvalue;
    
    flux2.right.type='V';
    flux2.right.value=flux2.right.setvalue=alphag.right.setvalue*U.right.setvalue;




end