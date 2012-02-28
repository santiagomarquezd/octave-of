% calcVdrpSimple

% Implements gas velocity calculation with simple equation (ever positive)
Vg=U;
Vg.internal = V0.*(1-alphag0.internal);

% To adjust BC's to dictionary
Vg=setBC(Vg,rhom,xC,xF,g);




