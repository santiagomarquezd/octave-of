function [Uddt]=fvc_ddt(U, U0, dt)
    % Gives the temporal explicit derivative
    % for a given field
    % 
    % [Uddt]=fvc_ddt(U, U0)
    %
    % U: actual field value
    % U0: previous timestep field value
    % dt: time-step

    Uddt=(U.internal-U0.internal)/dt;
    dt;
end
