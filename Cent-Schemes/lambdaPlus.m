function [lambdaP]=lambdaPlus(lambda)
  % Gives the matrix of positive eigenvalues from full eigenvalues matrix
  %
  % [lambdaP]=lambdaPlus(lambda)
  %
  % lambdaP: array of matrices of positive eigenvalues
  % lambda: array of full eigenvalues matrices

  % Gets the number of matrices in array
  N=size(lambda,3);

  % Memory allocation
  lambdaP=lambda; 
  
  for i=1:N 
    lambdaP(:,:,i)=diag(diag(lambda(:,:,i)).*pos(diag(lambda(:,:,i))));
  end
end