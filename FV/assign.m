function [sol]=assign(first, second, op)
  % Applies given operation to first and second fields including BC's
  %
  % [sol]=assign(first, second, op)
  %
  % sol: solution field
  % first: first field to be added
  % second: second field to be added
  % op: operation ('+', '-', '*', '/','^')
  
  % Check for constant fields
%    if (not(isfield(second,'left')))
%  	second=constField(size(first,1),second);
%    end

  % Construct BC's
  sol.left.type='V';
  sol.right.type='V';
  
  % Do operations

  if (op=='+')
	% Internal field
	sol.internal=first.internal+second.internal;

	% Boundary conditions
	sol.left.value=first.left.setvalue+second.left.setvalue;
	sol.right.value=first.right.setvalue+second.right.setvalue;
	
  elseif (op=='-')
	% Internal field
	sol.internal=first.internal-second.internal;

	% Boundary conditions
	sol.left.value=first.left.setvalue-second.left.setvalue;
	sol.right.value=first.right.setvalue-second.right.setvalue;
   
  elseif (op=='*')
	% Internal field
	sol.internal=first.internal.*second.internal;

	% Boundary conditions
	sol.left.value=first.left.setvalue*second.left.setvalue;
	sol.right.value=first.right.setvalue*second.right.setvalue;
  
  elseif (op=='/')
	% Internal field
	sol.internal=first.internal./second.internal;

	% Boundary conditions
	sol.left.value=first.left.value/second.left.setvalue;
	sol.right.value=first.right.value/second.right.setvalue;
  
  elseif (op=='^')
	% Internal field
	sol.internal=first.internal.^second.internal;

	% Boundary conditions
	sol.left.value=first.left.setvalue^second.left.setvalue;
	sol.right.value=first.right.setvalue^second.right.setvalue;
	
  end
  
  % Evaluate BC's
  sol=setBC(sol,0,0,0,0);

end