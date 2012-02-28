function [w]=weights(xC, xF)
    % Gives the linear interpolation weights for internal faces
    % 
    % [w]=weights(xC, xF)
    %
    % w: linear interpolation weights
    % xC: cell centers
    % xF: face centers

	% Face interpolated field allocation
	w=abs(xC(2:end)-xF(2:end-1))./abs(xC(2:end)-xC(1:end-1));
end
