% Solves three waves of non-linear advection
% using various methods for non-linear hyperbolic
% equations
% du/dt+d/dx[F(u)]=0; 
% u=[u1;u2;u3], F(u)=[v1*u1;v2*u2;v3*u3]
% Where:
%       vi=vri*(1-ui)-sum_{j!=i}uj*vrj
%       vr=V0i*(1-beta)^ai
%       beta=sum_1^N ui 
%
% Solution of v_m, u and p are obtained by
% Finite Difference Equations obtained
% from the original Finite Volume formulation

% Variables clearance
clear all;
%close all;
page_screen_output(0);
warning ("off", "Octave:broadcast");

% Physical paramaters
% Gravity
g = -10;
% Dispersed phase density
rhop = 1; %1000; %CHANGED
% Continuous phase density
rhoq = 1000; %800/800;  %CHANGED
% Dispersed phase viscosity
mup = 1E-3;
% Continuous phase viscosity
muq = 0.3;


if (1)
    % Test case from thesis
    % Volume fractions 2 and 3 are set with
    % small values to represent only one 
    % dispersed phase case
    method='KT'; %'LxF';
    % Constants for advective velocities
    V0=[1;0;0]; %[-1;0;0]; %CHANGED
    % Exponents for advective velocities
    a=[1;1;1]; %a=[1;1;1];
    % Volume fraction of dense packed-layer
    alphaDPL=1;
    % Initial values
    % Two section iniatilization
    layers=1;
    layerL1=1;
    ULeft1=0.3;
    layerL2=1;
    ULeft2=1E-3;
    layerL3=1;
    ULeft3=1E-3;
    % Time-step
    dt= 0.001; %0.04; %0.001;
    % Number of timesteps
    timesteps=500;
    % Number of cells
    N=400;%10; %400;
    % Number of PISO correctors
    nCorrectors = 3;
else
    % Base case
    % Solver ENIEF 2014 Polydisperse paper
    % Method of integration (LxF, Rusanov, Roe, KT)
    method='KT'; %'LxF';
    % Constants for advective velocities
    V0=[-0.5;-1;-2];
    % Exponents for advective velocities
    a=[1;1;1]; %a=[1;1;1];
    % Volume fraction of dense packed-layer
    alphaDPL=0.74; %1;
    % Initial values
    % Two section iniatilization
    layers=1;
    layerL1=1;
    ULeft1=0.1; %0.4;
    layerL2=1;
    ULeft2=0.15; %0.2;
    layerL3=1;
    ULeft3=0.2; %0.1;
    % For LxF
    % dt = 0.001/25;
    % N = 400*25;
    % timesteps = 25000;
    % Time-step
    dt= 0.001;
    % Number of timesteps
    timesteps=1000;
    % Number of cells
    N=400;
end

% Domain extension (fixed to 1)
xleft=0;
xright=1;

% Dummy cross section
S = 1;

% BC's
%  vLeft1=1;
%  vLeft2=2;

% Numerical Pre-processing

% Cell centers and face centers calculation
% Equi-spaced cells
dx=(xright-xleft)/N;
xC=((xleft+dx/2):dx:(xright-dx/2))';
xF=(xleft:dx:xright)';
lambda=dt/dx;

% Face areas
Sf=ones(size(xF))*S;
% Cell volumes
Vol=ones(size(xC))*S*dx;
% Interpolation weights calculation
w=weights(xC, xF);

% Fields initialization
% u1
u1.internal=zeros(N,1);
u1.left.type='G';
u1.left.gradient=0;
u1.right.type='G';
u1.right.gradient=0;
if (layers)
  % Detect number of cells in left zone
  nCellsLL=sum(xC<layerL1);
  u1.internal(1:nCellsLL)=ULeft1;
end
u1=setBC(u1,constField(0,N),xC,xF,0);

% u2
u2.internal=zeros(N,1);
u2.left.type='G';
u2.left.gradient=0;
u2.right.type='G';
u2.right.gradient=0;
if (layers)
  % Detect number of cells in left zone
  nCellsLL=sum(xC<layerL2);
  u2.internal(1:nCellsLL)=ULeft2;
end
u2=setBC(u2,constField(0,N),xC,xF,0);

% u3
u3.internal=zeros(N,1);
u3.left.type='G';
u3.left.gradient=0;
u3.right.type='G';
u3.right.gradient=0;
if (layers)
  % Detect number of cells in left zone
  nCellsLL=sum(xC<layerL3);
  u3.internal(1:nCellsLL)=ULeft3;
end
u3=setBC(u3,constField(0,N),xC,xF,0);

% rho
rho.internal=u1.internal*rhop+(1-u1.internal)*rhoq;
rho.left.type='G';
rho.left.gradient=0;
rho.right.type='G';
rho.right.gradient=0;
rho=setBC(rho,constField(0,N),xC,xF,g);

% mu
mu.internal=u1.internal*mup+(1-u1.internal)*muq;
mu.left.type='G';
mu.left.gradient=0;
mu.right.type='G';
mu.right.gradient=0;
mu=setBC(mu,constField(0,N),xC,xF,g);

% p (fixed value at top, fixed gradient at bottom)
p.internal=zeros(N,1);
p.left.type='G';
p.left.gradient = 0; % To initialize with zero
p.right.type='V';
p.right.value=0;
p = setBC(p, constField(0,N), xC, xF, 0);

% V 
Vm.internal=zeros(N,1);
Vm.left.type='V';
Vm.left.value=0;
Vm.right.type='V';
Vm.right.value=0;
Vm=setBC(Vm,constField(0,N),xC,xF,g);

% Vqp field from V field
Vqp=Vm;
Vqp.left.type='V';
Vqp.left.value=0;
Vqp.right.type='V';
Vqp.right.value=0;
Vqp=setBC(Vqp,constField(0,N),xC,xF,g);

% Dummy phi initialization
phi = zeros(N + 1, 1);

% One Jacobian per inter-cell
A=zeros(3,3,N-1);

% Temporal loop
for i=1:timesteps

    % Prints present timestep
    i

    if (strcmp(method,'LxF'))

        % Temporal data
        u1tmp=u1.internal;
        u1left=u1.left.setvalue;
        u1right=u1.right.setvalue;

        u2tmp=u2.internal;
        u2left=u2.left.setvalue;
        u2right=u2.right.setvalue;

        u3tmp=u3.internal;
        u3left=u3.left.setvalue;
        u3right=u3.right.setvalue;

        usum=u1.internal+u2.internal+u3.internal;

        if (0)
            % Standard implementation

            % Flux calculations
            F1=(V0(1,1).*(1-u1tmp).*(alphaDPL-usum).^a(1,1)-1*(V0(2,1).*(alphaDPL-usum).^a(2,1).*u2tmp+V0(3,1).*(alphaDPL-usum).^a(3,1).*u3tmp)).*u1tmp;
            F2=(V0(2,1).*(1-u2tmp).*(alphaDPL-usum).^a(2,1)-1*(V0(1,1).*(alphaDPL-usum).^a(1,1).*u1tmp+V0(3,1).*(alphaDPL-usum).^a(3,1).*u3tmp)).*u2tmp;
            F3=(V0(3,1).*(1-u3tmp).*(alphaDPL-usum).^a(3,1)-1*(V0(1,1).*(alphaDPL-usum).^a(1,1).*u1tmp+V0(2,1).*(alphaDPL-usum).^a(2,1).*u2tmp)).*u3tmp;
            
            % u1 temporal advancement
        
            % Non boundary cells
            u1.internal(2:N-1)=1/2*(u1tmp(3:N)+u1tmp(1:N-2))-dt/(dx/2)*(F1(3:N)-F1(1:N-2));

            % Boundary cells           
            % First cell
            u1.internal(1)=1/2*(u1tmp(1)+u1tmp(2))-dt/(dx/2)*(F1(1)+F1(2)); % Impermeable wall
            % Last cell
            u1.internal(N)=1/2*(u1tmp(N-1)+u1tmp(N))+dt/(dx/2)*(F1(N-1)+F1(N)); % Impermeable wall


            % u2 temporal advancement

            % Non boundary cells
            u2.internal(2:N-1)=1/2*(u2tmp(3:N)+u2tmp(1:N-2))-dt/(dx/2)*(F2(3:N)-F2(1:N-2));

            % Boundary cells           
            % First cell
            u2.internal(1)=1/2*(u2tmp(1)+u2tmp(2))-dt/(dx/2)*(F2(1)+F2(2)); % Impermeable wall
            % Last cell
            u2.internal(N)=1/2*(u2tmp(N-1)+u2tmp(N))+dt/(dx/2)*(F2(N-1)+F2(N)); % Impermeable wall

            % u3 temporal advancement

            % Non boundary cells
            u3.internal(2:N-1)=1/2*(u3tmp(3:N)+u3tmp(1:N-2))-dt/(dx/2)*(F3(3:N)-F3(1:N-2));

            % Boundary cells           
            % First cell
            u3.internal(1)=1/2*(u3tmp(1)+u3tmp(2))-dt/(dx/2)*(F3(1)+F3(2)); % Impermeable wall
            % Last cell
            u3.internal(N)=1/2*(u3tmp(N-1)+u3tmp(N))+dt/(dx/2)*(F3(N-1)+F3(N)); % Impermeable wall
        else
            % Implementation using functions
            %[u1.internal,u2.internal,u3.internal]=LxFPoly(u1.internal,u2.internal,u3.internal,V0,a,alphaDPL,dx,dt);
        
            if (1)

            % Runke-Kutta 4
            % Data allocation
            k1_1=u1.internal;
            k2_1=u1.internal;
            k3_1=u1.internal;
            k4_1=u1.internal;

            k1_2=u2.internal;
            k2_2=u2.internal;
            k3_2=u2.internal;
            k4_2=u2.internal;

            k1_3=u3.internal;
            k2_3=u3.internal;
            k3_3=u3.internal;
            k4_3=u3.internal;

            [k1_1,k1_2,k1_3]=RKLxFPoly(u1.internal,u2.internal,u3.internal,V0,a,alphaDPL,dx,dt);

            [k2_1,k2_2,k2_3]=RKLxFPoly(u1.internal+dt/2*k1_1,u2.internal+dt/2*k1_2,u3.internal+dt/2*k1_3,V0,a,alphaDPL,dx,dt);

            [k3_1,k3_2,k3_3]=RKLxFPoly(u1.internal+dt/2*k2_1,u2.internal+dt/2*k2_2,u3.internal+dt/2*k2_3,V0,a,alphaDPL,dx,dt);

            [k4_1,k4_2,k4_3]=RKLxFPoly(u1.internal+dt*k3_1,u2.internal+dt*k3_2,u3.internal+dt*k3_3,V0,a,alphaDPL,dx,dt);

            u1.internal=u1.internal+dt/6*(k1_1+2*k2_1+2*k3_1+k4_1);
            u2.internal=u2.internal+dt/6*(k1_2+2*k2_2+2*k3_2+k4_2);
            u3.internal=u3.internal+dt/6*(k1_3+2*k2_3+2*k3_3+k4_3);
            
            end    

        end

    elseif (strcmp(method,'Rusanov'))

        % Rusanov solution

        % Fluxes
        % Cell centered fluxes
        % One flux vector per cell
        u=[(u1.internal)' 
            (u2.internal)' 
            (u3.internal)'];

        F=arrayPFlux(u,V0,a,alphaDPL);
        % Face fluxes (impermeable walls)
        F=[[0; 0; 0] (F(:,1:(end-1))+F(:,2:(end)))/2 [0; 0; 0]];

        % Arrays for all inter-cells (one advection matriz per inter-cell)
        % Vectorized version
        u=[((u1.internal(1:end-1)+u1.internal(2:end))/2)' 
            ((u2.internal(1:end-1)+u2.internal(2:end))/2)' 
            ((u3.internal(1:end-1)+u3.internal(2:end))/2)'];   
        
        A=arrayPFluxJacobian(u,V0,a,alphaDPL);

        lambda=arrayMaxAbsEig(A);
        LAMBDA=[lambda';lambda';lambda'];

        % Rusanov fluxes (the boundary fluxes are left zero)
        u=[u1.internal u2.internal u3.internal]'; 
        F(:,2:end-1)=F(:,2:end-1)-1/2.*LAMBDA.*(u(:,2:end)-u(:,1:end-1));

        % Rusanov integration
        u=u-dt/dx*(F(:,2:end)-F(:,1:end-1));

        u1.internal=u(1,:)';
        u2.internal=u(2,:)';
        u3.internal=u(3,:)';

    elseif (strcmp(method,'Roe'))

        % Roe solution

        % Fluxes
        % Cell centered fluxes
        % One flux vector per cell
        u=[(u1.internal)' 
            (u2.internal)' 
            (u3.internal)'];

        F=arrayPFlux(u,V0,a,alphaDPL);
        % Face fluxes (impermeable walls)
        F=[[0; 0; 0] (F(:,1:(end-1))+F(:,2:(end)))/2 [0; 0; 0]];

        % Arrays for all inter-cells (one advection matriz per inter-cell)
        % Vectorized version
        u=[((u1.internal(1:end-1)+u1.internal(2:end))/2)' 
            ((u2.internal(1:end-1)+u2.internal(2:end))/2)' 
            ((u3.internal(1:end-1)+u3.internal(2:end))/2)'];   
        

        A=arrayPFluxJacobian(u,V0,a,alphaDPL);

        [V,VT,LAMBDA]=arrayEig(A);

        % Quartepelle p. 81
        % V^(-1)*A*V=LAMBDA
        % V*LAMBDA*V^(-1)=A
        % V is right eigenvector
        % Left eigenvector L=R^(-1)

        % Some complex eigenvalues and eigenvector appear
        %ARoe=matrixArrayProd(matrixArrayProd(real(V),abs(LAMBDA)),matrixArrayInverse(real(V)));
        ARoe=matrixArrayProd(matrixArrayProd(V,abs(LAMBDA)),matrixArrayInverse(V));        

        % Roe fluxes (the boundary fluxes are left zero)
        for j=1:N-1
            uL=[u1.internal(j,1); u2.internal(j,1); u3.internal(j,1)]; 
            uR=[u1.internal(j+1,1); u2.internal(j+1,1); u3.internal(j+1,1)]; 
            F(:,j+1)=F(:,j+1)-1/2*ARoe(:,:,j)*(uR(:,1)-uL(:,1));
        end

        % Data reshaping
        u=[u1.internal u2.internal u3.internal]'; 

        % Roe integration
        u=u-dt/dx*(F(:,2:end)-F(:,1:end-1));

        % Data reshaping
        u1.internal=u(1,:)';
        u2.internal=u(2,:)';
        u3.internal=u(3,:)';

    elseif (strcmp(method,'KT'))

        % Kurganov and Tadmor solution

        % Limiting
        % Sweby's fuction calculation

        % phiAlphag=superbee(rvalue(u,1E-9));
        % phiAlphag=vanLeer(rvalue(u,1E-9));
    
        if (1)    
            phiU1=vanLeer(rvalue(u1,1E-9)); 
            phiU2=vanLeer(rvalue(u2,1E-9));
            phiU3=vanLeer(rvalue(u3,1E-9));
        else
            phiU1=minmod(rvalue(u1,1E-12))*0; % Constant values by cells (mimiking Rusanov?)
            phiU2=minmod(rvalue(u2,1E-12))*0;
            phiU3=minmod(rvalue(u3,1E-12))*0;
        end    


        % Limited values calculation
        [u1Limited]=limitedValues(u1,phiU1,dx,dt);
        [u2Limited]=limitedValues(u2,phiU2,dx,dt);
        [u3Limited]=limitedValues(u3,phiU3,dx,dt);

        % m_h: from left boundary to last intercell
        % p_h: from first intercell to right boundary

        % Left and right Jacobians at intercells (N-1 intercells)
        uMinus=[(u1Limited.u_i_p_h_l(1:N-1))' 
                (u2Limited.u_i_p_h_l(1:N-1))' 
                (u3Limited.u_i_p_h_l(1:N-1))'];

        uPlus=[(u1Limited.u_i_p_h_r(1:N-1))' 
                (u2Limited.u_i_p_h_r(1:N-1))' 
                (u3Limited.u_i_p_h_r(1:N-1))'];

        AMinus=arrayPFluxJacobian(uMinus,V0,a,alphaDPL);
        APlus=arrayPFluxJacobian(uPlus,V0,a,alphaDPL);

        % Filter Jacobians (use only diagonal part of jacobian)
        if (0)
            for j=1:size(AMinus,3)
                AMinus(:,:,j)=diag(diag(AMinus(:,:,j)));
                APlus(:,:,j)=diag(diag(APlus(:,:,j)));
            end
        end

        % Left and right spectral radius
        if (1)
            % Exact
            lambdaMinus=arrayMaxAbsEig(AMinus);
            lambdaPlus=arrayMaxAbsEig(APlus);
        else
            % Estimate by matrix norm
            lambdaMinus=arrayMatrixNorm(AMinus,1);
            lambdaPlus=arrayMatrixNorm(APlus,1);
        end

        % K&T local speeds at intercells
        aspeeds=max([lambdaMinus';lambdaPlus']);
        
%              uMinusHL=[(u1Limited.u_i_m_h_l)' 
%                        (u2Limited.u_i_m_h_l)' 
%                        (u3Limited.u_i_m_h_l)'];
%      
%              uMinusHR=[(u1Limited.u_i_m_h_r)' 
%                        (u2Limited.u_i_m_h_r)' 
%                        (u3Limited.u_i_m_h_r)'];

%              uPlusHL=[(u1Limited.u_i_p_h_l)' 
%                       (u2Limited.u_i_p_h_l)' 
%                       (u3Limited.u_i_p_h_l)'];
%      
%              uPlusHR=[(u1Limited.u_i_p_h_r)' 
%                       (u2Limited.u_i_p_h_r)' 
%                       (u3Limited.u_i_p_h_r)'];
                    

%              FMinusHL=arrayPFlux(uMinusHL,V0,a,alphaDPL);
%              FMinusHR=arrayPFlux(uMinusHR,V0,a,alphaDPL);
%              FPlusHL=arrayPFlux(uPlusHL,V0,a,alphaDPL);
%              FPlusHR=arrayPFlux(uPlusHR,V0,a,alphaDPL);

        % Fluxes at intercells
        FMinus=arrayPFlux(uMinus,V0,a,alphaDPL);
        FPlus=arrayPFlux(uPlus,V0,a,alphaDPL);

        ASPEEDS=[aspeeds;aspeeds;aspeeds];

        F=1/2*(FMinus+FPlus)-1/2*ASPEEDS.*(uPlus-uMinus);

        % Final face fluxes (impermeable walls)
        F=[[0; 0; 0] F [0; 0; 0]];

        % Data reshaping        
        u=[u1.internal u2.internal u3.internal]'; 

        % K&T integration
        u=u-dt/dx*(F(:,2:end)-F(:,1:end-1));  %CHANGED

        % Data reshaping        
        u1.internal=u(1,:)';
        u2.internal=u(2,:)';
        u3.internal=u(3,:)';

    end % Methods selection

    % Apply BC's
    u1=setBC(u1,constField(0,N),xC,xF,0);
    u2=setBC(u2,constField(0,N),xC,xF,0);
    u3=setBC(u3,constField(0,N),xC,xF,0);

    % ************************************************************
    % Momentum predictor 
    % ************************************************************
    
    if 0
        % Override u1 values from FOAM
        u1.internal = [0.2198892; 0.2961108; 0.3; 0.3; 0.3; 0.3; 0.3; 0.3; 0.3001764; 0.3838236];
    end
    
    % Update rho storing last rho
    rho0 = rho;
    rho.internal=u1.internal*rhop+(1-u1.internal)*rhoq;
    rho=setBC(rho,constField(0,N),xC,xF,g);
    
    % Update mu
    mu.internal=u1.internal*mup+(1-u1.internal)*muq;
    mu=setBC(mu,constField(0,N),xC,xF,g);
    
    % In order to assemble U equation the center-of-mass velocity face flux
    % is required. It is calculated only for one dispersed phase (u1)
    % rhoPhi = phiAlphaQ*(rhoq - rhop) + phi*rhop
    % 
    % where:
    %
    % rhoPhi: is the center-of-mass velocity face flux by rho
    % phiAlphaQ: is the total face flux of the alpha equation (F)
    % rhoq, rhop: are the fluid densities
    % phi: is the center-of-volume velocity face flux
    % This formula is for primary phase advection equation
    % For secondary phase flux is in the opposite direction
    rhoPhi = -(F(1, :)'*(rhoq - rhop) + phi*rhop);
    
    if 0
        % Override rhoPhi values from FOAM
        rhoPhi = [0, -200.076723, -209.79, -209.79, -209.79, -209.79, -209.79, -209.79, -209.79, -209.349441, 0]';
    end
    
    % Precomputed quantities
    % Mixture density at faces
    rhof = fvc_interpolate(rho, w, xC, xF);
    % Mixture viscosities at faces
    muf = fvc_interpolate(mu, w, xC, xF);
    
    vd = zeros(N, 1);
    vl = vu = zeros(N - 1, 1);
    rowd = (1:N)';
    columnd = (1:N)';
    % Assemble internal matrix' equations
    % using sparse matrices
    % Diagonal values (a_P)
    vd(2:end-1) = rho.internal(2:end-1).*Vol(2:end-1)./dt + ...
                    abs(rhoPhi(3:end-1)).*sign(rhoPhi(3:end-1) > 0) + abs(rhoPhi(2:end-2)).*sign(rhoPhi(2:end-2) < 0); % + ...                      %CHANGED
                    %(muf(3:end-1) + muf(2:end-2))/dx;          %CHANGED
    % Low triangular off diagonal values
    vl = -abs(rhoPhi(2:end-1)).*sign(rhoPhi(2:end-1) > 0); % - muf(2:end-1)/dx; %CHANGED
    rowl = (2:N)';
    columnl = (1:N-1)';
    % Upper triangular off diagonal values
    vu = -abs(rhoPhi(2:end-1)).*sign(rhoPhi(2:end-1) < 0); %-muf(2:end-1)/dx; %CHANGED
    rowu = (1:N-1)';
    columnu = (2:N)';
    % Assembling
    M = sparse([rowd; rowl; rowu], [columnd; columnl; columnu], [vd; vl; vu]);
    
    % Assemble base RHS
    % Gravitational term
    b = g.*Vol.*rho.internal;
    % Temporal term
    % Isolated assembling for the use with H operator
    bt = rho0.internal.*Vm.internal.*Vol(1:end)./dt;
    b += bt;
    % Presure gradient term (internal cells)
    b(2:end-1) -= (p.internal(3:end) - p.internal(1:end-2))/2;
    
    % Boundary conditions effect
    % In the matrix
    % Considering impermeable walls (V=0)
    M(1, 1) = rho.internal(1).*Vol(1)./dt + ...
                    abs(rhoPhi(2))*sign(rhoPhi(2) > 0); % + ... %CHANGED
                    %muf(2)/dx + muf(1)/(dx/2); %CHANGED
    M(1, 2) -= muf(2)/dx*0; %CHANGED
    M(N, N) = rho.internal(N).*Vol(N)./dt + ...
                    abs(rhoPhi(end - 1))*sign(rhoPhi(end - 1) < 0); % + ... %CHANGED
                    %muf(end)/(dx/2) + muf(end - 1)/dx; %CHANGED
    M(N, N - 1) += -muf(end - 1)/dx*0; %CHANGED
    % In the RHS
    % From viscous term
    % (Should be null)
    %b(1) += muf(1)*V.left.setvalue/(dx/2); %CHANGED
    %b(end) += muf(end)*V.right.setvalue/(dx/2); %CHANGED
    % From pressure
    b(1) -= (p.internal(2) + p.internal(1))/2 - p.left.setvalue;
    b(end) -= p.right.setvalue - (p.internal(end - 1) + p.internal(end))/2;

    % FOAM like assembling
    [ddtM, ddtRHS] = fvm_ddt(rho, rho0, Vm, Vol, dt, 1);
    [convM, convRHS] = fvm_div_flux_cell(rhoPhi, Vm, xC, xF, w, 1);
    
    % Solve for V
    Vm.internal = M\b;
    if 0
        % Override velocity with exact solution
        Vm.internal = -u1.internal*(rhop/rho.internal - 1)*Vqp.internal;
    end
    % Set boundary conditions
    Vm = setBC(Vm,constField(0,N),xC,xF,g);
    
    % Drift flux needed for principal velocities fluxes
    % transformation formulas
    u1tmp=u1.internal;
    u2tmp=u2.internal;
    u3tmp=u3.internal;
    usum=u1.internal+u2.internal+u3.internal;
    % Valid only for one dispersed phase
    Vqp.internal = -V0(1,1).*(alphaDPL-usum).^a(1,1);
    phir = fvc_interpolate(Vqp, w, xC, xF);
    % Multiplying factor
    mFactor = assign(assign(assign(assign(constField(1,N), u1, '-'), u1, '*'), constField(rhoq - rhop, N), '*'), rho, '/');  
    driftPhi = fvc_interpolate(mFactor, w, xC, xF).*phir;
    if 1
        % By direct interpolation of the product
        driftPhi = fvc_interpolate(assign(mFactor,Vqp,'*'), w, xC, xF);
        % Extracting it from the volume fraction integrator flux
        %driftPhi = (-F(1, :)' - phi)*(rhoq - rhop)./fvc_interpolate(rho, w, xC, xF);
    end
    
    % ************************************************************
    % Solve pressure by PISO loop
    % ************************************************************
    for k = 1:nCorrectors
        % Calculate rUA
        rUA = arrayToField(1./Aop(M, Vol), N);
        % Calculate HbyA
        HopArray = Hop(M, bt, Vm, Vol);
        HbyA = assign(arrayToField(HopArray, N), rUA, '*');
        % Correct boundary conditions to mimic V
        HbyA.left.type='V';
        HbyA.left.value=0;
        HbyA.right.type='V';
        HbyA.right.value=0;
        HbyA = setBC(HbyA, constField(0,N), xC, xF, g);
        % Diffusivities at faces
        rUAf = fvc_interpolate(rUA, w, xC, xF); 
        % First approximation of center-of-mass velocity flux
        phiHbyA = fvc_interpolate(HbyA, w, xC, xF);
        if 1
            % Override phiHbyA
            directionFlux = sign(HbyA.internal(2:end) - HbyA.internal(1:end - 1))*(-1);
            %directionFlux = -1*ones(N - 1, 1);
            phiHbyA = fvc_general_interpolate(HbyA, xC, xF, -1, directionFlux);
        end
%          % Override rho ************************************
%          rho.internal(:, 1)=rhop;
%          rho=setBC(rho,constField(0,N),xC,xF,g);
        
        % Gravity effect added
        if 0
            phig = fvc_interpolate(rho, w, xC, xF).*g.*rUAf;
        else
            phig = fvc_interpolate(assign(rho, rUA, '*'), w, xC, xF).*g;
        end
        phiHbyA += phig;

        if 0
            % Override data from FOAM
            %rUAf = [1.4249741513629238e-06, 1.4249727153552731e-06, 1.4249712793476223e-06, 1.4249712793476221e-06, 1.4249712793476221e-06, 1.4249712793476221e-06, 1.4249712793476221e-06, 1.4249712793476221e-06, 1.4249724530670622e-06, 1.4264665363691187e-06, 1.427959445951735e-06]';
            
            %phiHbyA = [0.01, -0.010011510015425865, -0.010000011546346783, -0.010000000002468871, -0.0099999999999999985, -0.0099999999999997244, -0.0099999999998687111, -0.0099999999372612017, -0.0099999690546100101, -0.0099884418019184473, -0.0099790464952536247]';
            
            driftPhi = [0, -0.25764885449591207, -0.29763095546053775, -0.29957161216621453, -0.29957161216621453, -0.29957161216621453, -0.29957161216621453, -0.29957161216621453, -0.29965963227513331, -0.34147436414064131, 0]';
        end
        
        % Assemble pressure equation
        % Diagonal contribution
        % Internal cells
        vd(2:end-1) = -rUAf(3:end - 1)/dx - rUAf(2:end - 2)/dx;
        % Upper triangular off-diagonal contribution
        vu = rUAf(2:end - 1)/dx;
        % Lower triangular off-diagonal contribution
        vl = rUAf(2:end - 1)/dx;
        % Matrix assembling
        P = sparse([rowd; rowl; rowu], [columnd; columnl; columnu], [vd; vl; vu]);
        
        % Assemble base RHS
        % Divergence of phi tilde
        bp = phiHbyA(2:end) - phiHbyA(1:end - 1) - (driftPhi(2:end) - driftPhi(1:end -1));
        % To assure flux consistency in gradient BC's for p
        bp(1) += phiHbyA(1);
        
        % Boundary conditions effects
        % In the matrix
        P(1, 1) = - rUAf(2)/dx;
        P(N, N) = - rUAf(end - 1)/dx - rUAf(end)/(dx/2);
        % In the RHS, assembled in a different vector in order to 
        % evaluate the p matrix flux correctly
        bbp = bp*0;
        %bbp(1) += rUAf(1)*phig(1);      % Este termino en FOAM por alguna razon no está
        bbp(end) -= rUAf(end)*p.right.setvalue/(dx/2);
        bp += bbp;
        
        
        % pEqn initial residual
        pRes = norm(P*p.internal - bp);
        printf("Residual for p %g:\n", pRes);
        % Solve pEqn
        p.internal = P\bp;
        p.left.type='G';
        p.left.gradient = rho.left.setvalue*g;   % ES MAS O MENOS ACA, POR QUE DEPENDE DE COMO SEA EL GRADIENTE ADENTRO
        p=setBC(p, constField(0,N), xC, xF, 0);
        
        % Flux correction
        pFlux = flux(P, bbp, p);
        % Recover flux at gradient BC boundary
        pFlux(1) += phiHbyA(1);
        phiV = phiHbyA - pFlux;
        
        % Velocity of center-of-mass correction
        gradp = zeros(N, 1);
        gradp(2:end - 1, 1) = (p.internal(3:end) - p.internal(1:end - 2))/(2*dx);
        gradp(1) = ((p.internal(2) + p.internal(1))/2 - p.left.setvalue)/dx;
        gradp(end) = (p.right.setvalue - (p.internal(end) + p.internal(end - 1))/2)/dx;
        Vm.internal = HbyA.internal - rUA.internal.*gradp; 
        if 1
            % Cell centered correction for gravity
            Vm.internal = Vm.internal + rho.internal.*g.*rUA.internal;
        else    
            % Extended stencil correction for gravity
            rhoCorr = gradp*0;
            rhoCorr(2:end - 1) = rho.internal(1:end - 2)/4 + rho.internal(2:end -1 )/2 + ...
                                rho.internal(3:end)/4;
            rhoCorr(1) = rho.left.setvalue/2 + (rho.internal(1) + rho.internal(2))/4;
            rhoCorr(end) = (rho.internal(N - 1) + rho.internal(N))/4 + rho.right.setvalue/2;
            Vm.internal = Vm.internal + rhoCorr.*g.*rUA.internal;
        end
        % Correct boundary fields
        Vm = setBC(Vm, constField(0,N), xC, xF, 0);

        % Velocity of center-of-volume flux correction
        phi = phiV - driftPhi;
        
        % Mass imbalance checking 
        massI = phi(2:end) - phi(1:end - 1);
        
        printf("Sum local cont. errors %g:\n", mean(abs(massI))); 
        printf("Global cont. errors %g:\n", mean(massI)); 
        
    end

end %Temporal lop

%plot(xC,u1.internal,'m'); hold on; plot(xC,u2.internal, 'r'); plot(xC,u3.internal, 'c');plot(xC,u1.internal+u2.internal+u3.internal, 'k'); plot(xC,1-u1.internal-u2.internal-u3.internal, 'b')
figure(2);
plot(xC,1-u1.internal-u2.internal-u3.internal,";petroleum;",xC,u1.internal,";water1;",xC,u2.internal,";water2;",xC,u3.internal,";water3;",xC,u1.internal+u2.internal+u3.internal,'color','black',";water;",'linewidth',3);
%close all; figure(1); plot(xC,u1.internal,'m'); hold on; plot(xC,u2.internal, 'r'); plot(xC,u3.internal, 'c');plot(xC,u1.internal+u2.internal+u3.internal, 'k'); plot(xC,1-u1.internal-u2.internal-u3.internal, 'b')
%plot(xC,u1.internal,'m');
axis([0 1 0 1])

%u=[xC u1.internal u2.internal u3.internal]; save -ascii alphaDPL_0_74_LxF_t_60_ref.dat u

% eval(["print('-dtex', 'profile_a_" num2str(a) "_" num2str(b) "_d.tex')"]);
% print('-djpg', 'Rusanov_2_5.jpg')
% print('-djpg', 'LxF_60.jpg')

