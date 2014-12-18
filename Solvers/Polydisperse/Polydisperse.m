% Solves three waves of non-linear advection
% using various methods for non-linear hyperbolic
% equations
% du/dt+d/dx[F(u)]=0; 
% u=[u1;u2;u3], F(u)=[v1*u1;v2*u2;v3*u3]
% Where:
%       vi=vri*(1-ui)-sum_{j!=i}uj*vrj
%       vr=V0i*(1-beta)^ai
%       beta=sum_1^N ui 

% Variables clearance
clear all;
%close all;
page_screen_output(0);
warning ("off", "Octave:broadcast");

% Physical paramaters
% Dummy value, only for methods compatibility
g=0;

if (0)
    % Mono-polydispersity report case
    % Method of integration (LxF, Rusanov, Roe, KT)
    method='LxF';

    % Constants for advective velocities
    V0=[-0.5;-1;-2];

    % Exponents for advective velocities
    a=[1;1;1]; %a=[1;1;1];

    % Volume fraction of dense packed-layer
    alphaDPL=1; %0.7;

    % Initial values
    % Two section iniatilization
    layers=1;
    layerL1=1;
    ULeft1=0.4;
    layerL2=1;
    ULeft2=0.2;
    layerL3=1;
    ULeft3=0.1;

    % Time-step
    dt=0.001/5; %0.0001

    % Number of timesteps
    timesteps=2500*5; %2500;

    % Number of cells
    N=400*5;
else
    % Test case
    % Method of integration (LxF, Rusanov, Roe, KT)
    method='LxF'; %'LxF';
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
    % Time-step
    dt=0.001/25; %0.001;
    % Number of timesteps
    timesteps=8000/8*25; %8000;
    % Number of cells
    N=400*25; %400;
end

% Domain extension
xleft=0;
xright=1;

% BC's
vLeft1=1;
vLeft2=2;

% Numerical Pre-processing

% Cell centers and face centers calculation
% Equi-spaced cells
dx=(xright-xleft)/N;
xC=((xleft+dx/2):dx:(xright-dx/2))';
xF=(xleft:dx:xright)';
lambda=dt/dx;

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
            u1.internal(2:N-1)=1/2*(u1tmp(3:N)+u1tmp(1:N-2))-dt/dx/2*(F1(3:N)-F1(1:N-2));

            % Boundary cells           
            % First cell
            u1.internal(1)=1/2*(u1tmp(1)+u1tmp(2))-dt/dx/2*(F1(1)+F1(2)); % Impermeable wall
            % Last cell
            u1.internal(N)=1/2*(u1tmp(N-1)+u1tmp(N))+dt/dx/2*(F1(N-1)+F1(N)); % Impermeable wall


            % u2 temporal advancement

            % Non boundary cells
            u2.internal(2:N-1)=1/2*(u2tmp(3:N)+u2tmp(1:N-2))-dt/dx/2*(F2(3:N)-F2(1:N-2));

            % Boundary cells           
            % First cell
            u2.internal(1)=1/2*(u2tmp(1)+u2tmp(2))-dt/dx/2*(F2(1)+F2(2)); % Impermeable wall
            % Last cell
            u2.internal(N)=1/2*(u2tmp(N-1)+u2tmp(N))+dt/dx/2*(F2(N-1)+F2(N)); % Impermeable wall

            % u3 temporal advancement

            % Non boundary cells
            u3.internal(2:N-1)=1/2*(u3tmp(3:N)+u3tmp(1:N-2))-dt/dx/2*(F3(3:N)-F3(1:N-2));

            % Boundary cells           
            % First cell
            u3.internal(1)=1/2*(u3tmp(1)+u3tmp(2))-dt/dx/2*(F3(1)+F3(2)); % Impermeable wall
            % Last cell
            u3.internal(N)=1/2*(u3tmp(N-1)+u3tmp(N))+dt/dx/2*(F3(N-1)+F3(N)); % Impermeable wall
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
    
        if (0)    
            phiU1=vanLeer(rvalue(u1,1E-9)); % Constant values by cells (mimiking Rusanov?)
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

        % Filter Jacobians
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
        u=u-dt/dx*(F(:,2:end)-F(:,1:end-1));

        % Data reshaping        
        u1.internal=u(1,:)';
        u2.internal=u(2,:)';
        u3.internal=u(3,:)';

    end % Methods selection

    % Apply BC's
    u1=setBC(u1,constField(0,N),xC,xF,0);
    u2=setBC(u2,constField(0,N),xC,xF,0);
    u3=setBC(u3,constField(0,N),xC,xF,0);

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

