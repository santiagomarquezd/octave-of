%
%  add to path my directories needed for Matlab computations
%

% aqui poner en HOME el directorio del cual cuelgan estas cosas 
% en mi caso, mi maquina yo tengo esto
%HOME = '/home/nnigro/Soft';
% con la distro esta en mi maquina estara aqui
%HOME = '/home/nnigro/Soft/MY_TOOLS/femcode/version4.1/version4.1-distro';
%HOME = '';
HOME = 'C:\Users\Sebastian Corso\Documents\octave-of';
% en la maquina de uds deben poner en HOME lo que equivalga a
%     /home/nnigro/Soft/MY_TOOLS/femcode/version4.1
% en su directorio

eval(['addpath ''' HOME '/Solvers\twoPhaseEuler'''])
eval(['addpath ''' HOME '/FV'''])
eval(['addpath ''' HOME '/Cent-Schemes'''])


