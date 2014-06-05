% Solves two waves of non-linear advection
% using Flux Vector Splitting
% du/dt+d/dx[F(u)]=0; 
% u=[u1;u2], F(u)=[v1*u1;v2*u2]
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

if (1)
    % Mono-polydispersity report case

    % Constants for advective velocities
    V0=[-0.25;-0.5;-1];

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
    dt=0.0001; %0.0001

    % Number of timesteps
    timesteps=600000;

    % Number of cells
    N=10000;
else
    % Test case
    % Constants for advective velocities
    V0=[-0.5;-1;-2];
    % Exponents for advective velocities
    a=[1;1;1]; %a=[1;1;1];
    % Volume fraction of dense packed-layer
    alphaDPL=0.7/0.7;
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
    dt=0.001; %0.0001
    % Number of timesteps
    timesteps=60000; %2500;
    % Number of cells
    N=400;
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
%      usumleft=u1.left.setvalue+u2.left.setvalue;
%      usumright=u1.right.setvalue+u2.right.setvalue;

    if (0)
    
        % Decoupled

        % Flux calculations
        F1=V0(1,1).*(1-u1tmp).*(1-u1tmp).^a(1,1).*u1tmp;
        F2=V0(2,1).*(1-u2tmp).*(1-u2tmp).^a(2,1).*u2tmp;
        F3=V0(3,1).*(1-u3tmp).*(1-u3tmp).^a(3,1).*u3tmp;

    
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

        if (1)

            % Coupled

            % Decoupled
        
            % Flux calculations
            if (0)
                F1=(V0(1,1).*(1-u1tmp).*(1-usum).^a(1,1)-1*(V0(2,1).*(1-usum).^a(2,1).*u2tmp+V0(3,1).*(1-usum).^a(3,1).*u3tmp)).*u1tmp;
                F2=(V0(2,1).*(1-u2tmp).*(1-usum).^a(2,1)-1*(V0(1,1).*(1-usum).^a(1,1).*u1tmp+V0(3,1).*(1-usum).^a(3,1).*u3tmp)).*u2tmp;
                F3=(V0(3,1).*(1-u3tmp).*(1-usum).^a(3,1)-1*(V0(1,1).*(1-usum).^a(1,1).*u1tmp+V0(2,1).*(1-usum).^a(2,1).*u2tmp)).*u3tmp;

                u=[u1.internal u2.internal u3.internal]'; 
                F=arrayPFlux(u,V0,a);
                F1=F(1,:)';
                F2=F(2,:)';
                F2=F(2,:)';    
            else
                F1=(V0(1,1).*(1-u1tmp).*(alphaDPL-usum).^a(1,1)-1*(V0(2,1).*(alphaDPL-usum).^a(2,1).*u2tmp+V0(3,1).*(alphaDPL-usum).^a(3,1).*u3tmp)).*u1tmp;
                F2=(V0(2,1).*(1-u2tmp).*(alphaDPL-usum).^a(2,1)-1*(V0(1,1).*(alphaDPL-usum).^a(1,1).*u1tmp+V0(3,1).*(alphaDPL-usum).^a(3,1).*u3tmp)).*u2tmp;
                F3=(V0(3,1).*(1-u3tmp).*(alphaDPL-usum).^a(3,1)-1*(V0(1,1).*(alphaDPL-usum).^a(1,1).*u1tmp+V0(2,1).*(alphaDPL-usum).^a(2,1).*u2tmp)).*u3tmp;
            end
            
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

        end
    end

    % Apply BC's
    u1=setBC(u1,constField(0,N),xC,xF,0);
    u2=setBC(u2,constField(0,N),xC,xF,0);
    u3=setBC(u3,constField(0,N),xC,xF,0);

end

close all; plot(xC,u1.internal,'m'); hold on; plot(xC,u2.internal, 'r'); plot(xC,u3.internal, 'c');plot(xC,u1.internal+u2.internal+u3.internal, 'k'); plot(xC,1-u1.internal-u2.internal-u3.internal, 'b')
%figure(2); plot(xC,u1.internal,'m'); hold on; plot(xC,u2.internal, 'r'); plot(xC,u3.internal, 'c');plot(xC,u1.internal+u2.internal+u3.internal, 'k'); plot(xC,1-u1.internal-u2.internal-u3.internal, 'b')
axis([0 1 0 1])

% eval(["print('-dtex', 'profile_a_" num2str(a) "_" num2str(b) "_d.tex')"]);
% print('-djpg', 'Rusanov_2_5.jpg')
print('-djpg', 'LxF_60.jpg')
