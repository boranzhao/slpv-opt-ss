The data and codes are for the following paper:

P. Zhao and R. Nagamune, Switching linear-parameter-varying control with improved local performance and
optimized switching surfaces," International Journal of Robust and Nonlinear Control, accepted in Feb 2018.


Note that the code for Algorithm~1 in the paper is based on the particle swarm algorithm (PSO) code developed by Sam, which is available at https://www.mathworks.com/matlabcentral/fileexchange/25986-constrained-particle-swarm-optimization. The original PSO code is slightly modified for optimizing the switching surfaces in switching LPV controller design. The major revisions were done in pso.m (modified to be pso_ssopt.m), psoboundspenalize.m (modified to be psoboundspenalize_ssopt.m) and \private\psocreationuniform.m (modified to be \private\psocreationuniform_ssopt.m). 
