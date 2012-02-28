% volVectorField Ur(Ua - Ub);
Ur = assign(Ua,Ub,'-');
% volScalarField magUr(mag(Ur));
magUr = assign(assign(Ur,Ur,'*'),constField(1/2,N),'^'); 
% 
% volScalarField Ka(draga->K(magUr));
% volScalarField K(Ka);
K = dragmodel_K(magUr,InterfacialProperties,transportProperties);
% 
% if (dragPhase == "b")
% {
%     volScalarField Kb(dragb->K(magUr));
%     K = Kb;
% }
% else if (dragPhase == "blended")
% {
%     volScalarField Kb(dragb->K(magUr));
%     K = (beta*Ka + alpha*Kb);
% }
if (InterfacialProperties.dragPhase == 'b')
elseif (InterfacialProperties.dragPhase == 'blended')
end
% 
% volVectorField liftCoeff(Cl*(beta*rhob + alpha*rhoa)*(Ur ^ fvc::curl(U)));
liftCoef = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LiftCoeficient:
% La fuerza de Lift o lubricación es Nula en Casos 1D debido a que este
% termino actua en la direccion transversal: curl(U)=0
% volVectorField liftCoeff = Cl*(beta*rhob + alpha*rhoa)*(Ur ^fvc::curl(U));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
