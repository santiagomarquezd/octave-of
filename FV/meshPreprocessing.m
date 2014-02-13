% Cell centers and face centers calculation
% Equi-spaced cells
dx=(xright-xleft)/N;
xC=((xleft+dx/2):dx:(xright-dx/2))';
xF=(xleft:dx:xright)';
% Interpolation weights calculation
w=weights(xC, xF);
% Face areas
S=1;
Sf=ones(size(xF))*S;