% Run located in octave-of root
[a,b]=system("pwd");
eval(['addpath ' b(1:end-1) '/Solvers/twoPhaseEuler'])
eval(['addpath ' b(1:end-1) '/Solvers/mixture'])
eval(['addpath ' b(1:end-1) '/FV'])
eval(['addpath ' b(1:end-1) '/Cent-Schemes'])
eval(['addpath ' b(1:end-1) '/Solvers/testing'])
eval(['addpath ' b(1:end-1) '/Solvers/Traffic'])
eval(['addpath ' b(1:end-1) '/Solvers/alphaEqnIsolatedVm'])
eval(['addpath ' b(1:end-1) '/Solvers/alphaEqnIsolatedUm'])
eval(['addpath ' b(1:end-1) '/Solvers/alphaEqnIsolated'])
eval(['addpath ' b(1:end-1) '/Solvers/PISOSolver'])
eval(['addpath ' b(1:end-1) '/Solvers/SIMPLESolver'])


