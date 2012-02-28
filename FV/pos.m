function [out]=pos(vec)
    % Returns a vector with one values in positive entries
    % 
    % [out]=pos(vec)
    % 
    % vec: column vector of entering values
    % out: column vector with ones in positive values places 
    % and zeros otherwise

	out=vec>0;
end
