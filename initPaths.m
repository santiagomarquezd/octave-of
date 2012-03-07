% Run located in octave-of root
[a,b]=system("pwd");
eval(['addpath ' b(1:end-1) '/Solvers/twoPhaseEuler'])
eval(['addpath ' b(1:end-1) '/Solvers/mixture'])
eval(['addpath ' b(1:end-1) '/FV'])
eval(['addpath ' b(1:end-1) '/Cent-Schemes'])




