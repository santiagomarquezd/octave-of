% Solves an 1D one phase flow problem like OpenFOAM

clear all;
page_screen_output(0);

% Case initialization
test1

% **************************** MAIN PROGRAM ***************************

% Gravity treatment terms
ghf = g*xF;
gh = g*xC;

% Set fields as 'old' states
rho0 = rho;
U0 = U;

% phi field initialization from U
% Creation of flux direction
directionFlux = fvc_interpolate(U0, w, xC, xF);
directionFlux = sign(directionFlux(2:end-1,1)+1E-9);
phi = fvc_interpolate(U0, w, xC, xF).*Sf;

% Set fields as 'old' states
phi0 = phi;

if (initia != 1)
  load("dump.dat");
end

if 1 % Enables temporal loop

% Temporal loop
for step=1:timesteps
    fprintf('---------------------------------\n')
    fprintf('Timestep: %d. Time: %g\n',step,step*dt)

    % UEqn
    UEqn

    %PISO loop
    if 1  
        for corr = 1:nCorr
            pEqn
        end
    end

    % Set fields as 'old' states
    rho0 = rho;
    U0 = U;
    phi0 = phi;


    % Printing and saving
    if 0
        hold on;
            if (rem(step,100)==0 || step==1)
                %plot(xC,u.internal,'r*-')
                eval(['save mixtureSolver-UADE-' num2str(step) ...
                '.dat alphag rhom U p    xC        xF'])
        end
    end

end

end % Ends condition for temporal loop

