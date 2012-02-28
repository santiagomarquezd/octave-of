function [eigenvectors]=no_Vr_eigenvectors(U,alphag,rhol,rhog)
  % Gives the eigenvectors as matrix columns for no 
  % Vr mixture model
  % 
  % [eigenvectors]=no_Vr_eigenvectors(U,alphag)
  % 
  % eigenvectors: eigevector matrix with eigevectos as columns
  % U: mixture velocity field
  % alphag: void fraction for gas
  % rhol: liquid density
  % rhog: gas density

  eigenvectors=[1,0;-(alphag*rhom)/(rhol*Vm),1];
  
end
