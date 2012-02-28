function [out]=neg(vec)
    % Returns a vector with one values in negative entries
    % 
    % [out]=neg(vec)
    % 
    % vec: column vector of entering values
    % out: column vector with ones in negative values places 
    % and zeros otherwise

	out=vec<0;
end
