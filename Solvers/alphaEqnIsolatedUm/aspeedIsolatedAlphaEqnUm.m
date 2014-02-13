function [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnUm(alphal,Um,rhol,rhog,V0,a)
    % Gives the max of spectral radius at interfaces for 
    % Kurganov and Tadmor's kind of schemes
    % in simplified mixture problem for isolated alpha equation
    %
    % d alphal/dt+d(alphal*Um-alphal^a*(1-alphal)*Urlg)/dx=0
    %
    % [a_j_minus_half,a_j_plus_half]=aaspeedIsolatedAlphaEqn(alphal,Um,rhol,rhog,V0,a)
    %
    % a_j_minus_half: max of spectral radius at left interfaces
    % a_j_plus_half: max of spectral radius at right interfaces
    % alphal: liquid volume fraction
    % Um: velocity of center of volume
    % rhol: liquid density
    % rhog: gas density
    % V0: velocity constant for relative velocity law
    % a: exponent constant for relative velocity law
   
    
    % Data is flatten for an easier calculus (BC's and internal at the same time)
    alphal=[alphal.left.setvalue; alphal.internal; alphal.right.setvalue];
    Um=[Um.left.setvalue; Um.internal; Um.right.setvalue];

    % Eigenvalues calculation
    %(-V0*((a+1)*Foam::pow(alpha1_pos,a)-(a+2)*(alpha1_pos,a+1))+ui)+phiV
    eig1=-V0.*((a+1).*(alphal).^a-(a+2).*alphal.^(a+1))+Um;

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