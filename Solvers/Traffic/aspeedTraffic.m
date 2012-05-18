function [a_j_minus_half,a_j_plus_half]=aspeedTraffic(u,V0)
    % Gives the max of spectral radius at interfaces for 
    % Kurganov and Tadmor's kind of schemes
    % in Traffic Flow problem
    %
    %  [a_j_minus_half,a_j_plus_half]=aspeedTraffic(u,V0)
    %
    % a_j_minus_half: max of spectral radius at left interfaces
    % a_j_plus_half: max of spectral radius at right interfaces
    % u: solved unknown
    % V0: constant for velocity definition

 
    % u is flatten for an easier calculus (BC's and internal at the same time)
    u=[u.left.setvalue; u.internal; u.right.setvalue];

    eig1=V0-2*V0*u;
    eig2=V0-2*V0*u;

    % Local velocities at internal field
    a=max(abs(eig1),abs(eig2));
	  
    aLeft=a(1);
    aRight=a(end);
	
    % a is turned back to internal values only
    a=a(2:end-1);
    
    a_j_minus_half=[aLeft;max(a(1:end-1),a(2:end))];
    a_j_plus_half=[max(a(1:end-1),a(2:end));aRight];
	  
end	  