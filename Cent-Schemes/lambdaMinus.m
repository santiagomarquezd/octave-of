function [lambdaM]=lambdaMinus(lambda)
  % Gives the matrix of negative eigenvalues from full eigenvalues matrix
  %
  % [lambdaM]=lambdaMinus(lambda)
  %
  % lambdaP: array of matrices of negative eigenvalues
  % lambda: array of full eigenvalues matrices

  % Gets the number of matrices in array
  N=size(lambda,3);

  % Memory allocation
  lambdaM=lambda; 
  
  for i=1:N 
    lambdaM(:,:,i)=diag(diag(lambda(:,:,i)).*neg(diag(lambda(:,:,i))));
  end
end