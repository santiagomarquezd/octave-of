function [a_j_minus_half,a_j_plus_half]=aspeedShallowWater(h,hu,g)
    % Gives the max of spectral radius at interfaces for 
    % Kurganov and Tadmor's kind of schemes
    % in 1D Shallow Water problem
    %
    % [a_j_minus_half,a_j_plus_half]=aspeedShallowWater(h,hu,g)
    %
    % a_j_minus_half: max of spectral radius at left interfaces
    % a_j_plus_half: max of spectral radius at right interfaces
    % h: height
    % hu: momentum
    % g: gravitational acceleration

    % u calculation from conservative values
    u=assign(hu,h,'/');
   
    % h and u are flatten for an easier calculus (BC's and internal at the same time)
    h=[h.left.setvalue; h.internal; h.right.setvalue];
    u=[u.left.setvalue; u.internal; u.right.setvalue];

    eig1=u-sqrt(g*h);
    eig2=u+sqrt(g*h);

    % Local velocities at internal field
    a=max(abs(eig1),abs(eig2));
	  
    aLeft=a(1);
    aRight=a(end);
	
    % a is turned back to internal values only
    a=a(2:end-1);
    
    a_j_minus_half=[aLeft;max(a(1:end-1),a(2:end))];
    a_j_plus_half=[max(a(1:end-1),a(2:end));aRight];
	  
end	  