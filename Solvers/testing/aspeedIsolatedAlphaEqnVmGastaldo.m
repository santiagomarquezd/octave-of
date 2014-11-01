function [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnVmGastaldo(alphag,rhol,rhog,V0)
    % Gives the max of spectral radius at interfaces for 
    % Kurganov and Tadmor's kind of schemes
    % in simplified mixture problem for isolated alpha equation
    % Gastaldo's version of problem (constant Vpq)
    %
    % d alphag/dt+d(alphag*U+alphag*(1-cp)*Vr)/dx=0
    %
    % U is given by alphag*(rhog/rhom-1)*V0
    %
    % The eigenvalues are calculated without having rhom and Vr
    % wich can be calculated from theoretical formulas
    % rhom=rhog*alphag+(1-alphag)*rhol
    % Vr=V0*(1-alphag)
    %
    % [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnVmGastaldo(alphag,rhol,rhog,V0)
    %
    % a_j_minus_half: max of spectral radius at left interfaces
    % a_j_plus_half: max of spectral radius at right interfaces
    % alphag: gas volume fraction
    % rhol: liquid density
    % rhog: gas density
    % V0: velocity constant for relative velocity law
   
    
    % Data is flatten for an easier calculus (BC's and internal at the same time)
    alphag=[alphag.left.setvalue; alphag.internal; alphag.right.setvalue];

    %eig1=(((3.*alphag.^4-6.*alphag.^3+2.*alphag.^2+2.*alphag-1).*rhol.^2+(-6.*alphag.^4+8.*alphag.^3-2.*alphag.^2).*rhog.*rhol+(3.*alphag.^4-2.*alphag.^3).*rhog.^2).*V0)./((alphag.^2-2.*alphag+1).*rhol.^2+(2.*alphag-2.*alphag.^2).*rhog.*rhol+alphag.^2.*rhog.^2);
    eig1=-(((3.*alphag.^3-3.*alphag.^2-alphag+1).*rhol+(-3.*alphag.^3+3.*alphag.^2-alphag).*rhog).*V0)./((alphag-1).*rhol-alphag.*rhog);
    %eig1=(alphag>0.5)*(-1)+(alphag<0.5);
    %eig1=alphag.*0;

    disp('paso');

    eig2=eig1;

    % Local velocities at internal field
    if 0
      a=max(abs(eig1),abs(eig2));
    else
      a=max(abs(eig2),abs(eig2));
    end
      
    aLeft=a(1);
    aRight=a(end);
    
    % a is turned back to internal values only
    a=a(2:end-1);

    a_j_minus_half=[aLeft;max(a(1:end-1),a(2:end))];
    a_j_plus_half=[max(a(1:end-1),a(2:end));aRight];
	  
end	  