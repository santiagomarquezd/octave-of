for i=1:1

% Inverse of density at faces
rRhomF=fvc_interpolate(assign(constField(1,N),rhom,'/'),w,xC,xF);

% Velocity at faces for staggered grid scheme
rhomPhi = fvc_interpolate(rhom,w,xC,xF).*(fvc_interpolate(U,w,xC,xF).*Sf);
% For debugging purposes
rhomPhiP=rhomPhi;

% Assembling (implicit terms)
if (fullVerbose==1)
  printf('-----')
  disp('Assembling p_rghEqn')
end

% Diffusion coefficient
diffusion=ones(N+1,1);
[pM, pRHS]=fvm_laplacian(p,diffusion,xC,xF,Sf);

% For debugging purposes
AA=fvc_ddt(rhom, rhom0, dt);
BB=fvc_div_face(rhomPhi,V);

% Calculating explicit terms (to RHS)
ERHS=(fvc_ddt(rhom, rhom0, dt)+fvc_div_face(rhomPhi,V)).*V./dt;

% Solve
if (fullVerbose==1)
  disp('Solving for p_rgh')
end
p.internal=pM\(pRHS+ERHS);  

% U flux at faces correction
relax=0.5;
rhomPhi -= flux(pM,pRHS,p)*dt*relax;

% Reconstruction of velocities at cell centers
U.internal=fvc_reconstruct(rhomPhi.*rRhomF,Sf);
U=setBC(U,constField(0,N),xC,xF,g);

% Correct alphag
%  calcVdrp
%  BAlphaRFChorin

end
