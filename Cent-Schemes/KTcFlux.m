function [u]=KTcFlux(u0,uLimited,ufluxfunction,a,cFlux,uInt,dx,dt,rhom,rhom0,rhom0Int,V0,rhol,rhog,aexp)
    % Advances u field in time by means of Kurganov & Tadmor scheme
    % adding a centered flux
    %
    % [u]=KTcFlux(u0,uLimited,ufluxfunction,a,cFlux,uInt,dx,dt,Vpq,rhom,rhom0,rhom0Int,cp)
    %
    % u0: actual value of u
    % u: first field advanced in time
    % uLimited: limited values for u at actual tiem
    % ufluxfunction: flux function for u
    % a: local velocities at interfaces (eigenvalues)
    % cFlux: centered flux at faces
    % uInt: u variable interpolated at faces
    % dx: spatial step
    % dt: time step
    % Vpq: relative velocity
    % rhom: mixture density
    % rhom0: previous time-step mixture density
    % rhom0Int: previous time-step mixture density at faces, needed for stabilization
    % cp: mass fraction of secondary phase
    
    % Allocation
    u=u0;
    
    % Pre-processing
    lambda=dt/dx;

    % Stabilization tunning
    stab=1; %2.5;

    % Warning message about stabilization
    if (stab>1)
	  disp('***************************************');
	  disp('Stabilization factor is used en KTcFlux');
	  disp('***************************************');
    end
    
    % Since KTcFLux uses a flux given at faces the value of Um is zeroed
    N=size(u0.internal,1);
    Um=constField(0,N);
       
    % Star fluxes calculation (see http://en.wikipedia.org/wiki/MUSCL_scheme)
    if 0
      % Last formulation is Sao Carlos
      Fu_star_minus_half=1/2*(ufluxfunction(uLimited.u_i_m_h_r,V0,rhol,rhog,aexp)+...
						  ufluxfunction(uLimited.u_i_m_h_l,Um,V0,rhol,rhog,aexp)-...
						  rhom0Int(1:end-1).*a(1:end-1).*(uLimited.u_i_m_h_r-uLimited.u_i_m_h_l));
      Fu_star_plus_half=1/2*(ufluxfunction(uLimited.u_i_p_h_r,V0,rhol,rhog,aexp)+...
						  ufluxfunction(uLimited.u_i_p_h_l,Um,V0,rhol,rhog,aexp)-...
						  rhom0Int(2:end).*a(2:end).*(uLimited.u_i_p_h_r-uLimited.u_i_p_h_l));
    else
      % Formulation in B-45
      % Generation of mixture densities at faces with constant values by cells
      % First or last value is repeated like a ghost cell
      rhom0_i_p_h_r=[rhom0.internal(2:end);rhom0.internal(end)];
      rhom0_i_p_h_l=[rhom0.internal(1:end)];
      rhom0_i_m_h_r=[rhom0.internal(1:end)];
      rhom0_i_m_h_l=[rhom0.internal(1);rhom0.internal(1:end-1)];

      % F star flux generation
      Fu_star_plus_half=1/2*(ufluxfunction(uLimited.u_i_p_h_r,Um,V0,rhol,rhog,aexp)+...
						  ufluxfunction(uLimited.u_i_p_h_l,Um,V0,rhol,rhog,aexp)-...
						  stab.*a(2:end).*rhog.*(uLimited.u_i_p_h_r-uLimited.u_i_p_h_l));
      Fu_star_minus_half=1/2*(ufluxfunction(uLimited.u_i_m_h_r,Um,V0,rhol,rhog,aexp)+...
						  ufluxfunction(uLimited.u_i_m_h_l,Um,V0,rhol,rhog,aexp)-...
						  stab.*a(1:end-1).*rhog.*(uLimited.u_i_m_h_r-uLimited.u_i_m_h_l));
    end

    % Impermeable walls BC's (don't work at all if cflux if not zero at boundaries)
    if 1
      Fu_star_minus_half(1)=0;
      Fu_star_plus_half(end)=0;
    end

    %keyboard; pause;

    %Time advancement
    if 1
      u.internal=u0.internal.*rhom0.internal./rhom.internal-...
		lambda*(Fu_star_plus_half-Fu_star_minus_half)./rhom.internal-...
		lambda*(cFlux(2:end).*uInt(2:end)-cFlux(1:end-1).*uInt(1:end-1))./rhom.internal;
      %keyboard; pause;
    else
      u.internal=u0.internal-lambda*(Fu_star_plus_half-Fu_star_minus_half);
    end

end
