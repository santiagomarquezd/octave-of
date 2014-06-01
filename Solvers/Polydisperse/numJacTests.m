% Test for numerical jacobian calculation
% Simple coupled flux functions for polydisperse 
% cases
% F1:(V01*(1-alpha1)*(1-sumAlpha))*alpha1;
% F2:(V02*(1-alpha2)*(1-sumAlpha))*alpha2;
% sumAlpha=alpha1+alpha2
%
% See simpleCoupledVr

V01=-1;
V02=0.5;

alpha1=0.25;
alpha2=0.8;

% Exact jacobian
Je=[((2*alpha1-1)*alpha2+3*alpha1^2-4*alpha1+1)*V01 (alpha1^2-alpha1)*V01;
(alpha2^2-alpha2)*V02 (3*alpha2^2+(2*alpha1-4)*alpha2-alpha1+1)*V02]

% Numerical jacobian via forward differencing
alpha=[alpha1;alpha2];
sumAlpha=sum(alpha);
V=[V01;V02];
N=size(alpha,1);   
    
% Memory allocation
Jn=zeros(N,N);

% u variable lives in [0, 1]
dAlpha=1/1000; 

% Jacobian calculation loop
for i=1:N
    % Derivative of fluxes 1:N respect to 
    % variable i
    Falpha=(V.*(1-alpha).*(1-sumAlpha)).*alpha;
    alphaPlusdAlpha=alpha;
    alphaPlusdAlpha(i,1)=alpha(i,1)+dAlpha;
    sumAlphaPlusdAlpha=sum(alphaPlusdAlpha);
    FPlusDalpha=(V.*(1-alphaPlusdAlpha).*(1-sumAlphaPlusdAlpha)).*alphaPlusdAlpha;
    Jn(:,i)=(FPlusDalpha-Falpha)/dAlpha;
end

Jn

eig(Je)
eig(Jn)

