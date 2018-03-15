
function state = psocreationuniform(options,nvars,ss_dist_min)
% Generates uniformly distributed swarm based on options.PopInitRange.
% customized for switching surface optimization by Pan Zhao

n = options.PopulationSize;
itr = options.Generations;

[state,nbrtocreate] = psogetinitialpopulation(options,n,nvars) ;

% Initialize particle positions
state.Population(n-nbrtocreate+1:n,1) = ...
    repmat(options.PopInitRange(1,1),nbrtocreate,1) + ...
    repmat((options.PopInitRange(2,1) - options.PopInitRange(1,1)),...
    nbrtocreate,1).*rand(nbrtocreate,1) ;
% make sure the distance between switching surfaces are smaller than 
for i =2:nvars 
    state.Population(n-nbrtocreate+1:n,i) = ...        
    state.Population(n-nbrtocreate+1:n,i-1) + ss_dist_min +... 
    (repmat(options.PopInitRange(2,i),nbrtocreate,1)- state.Population(n-nbrtocreate+1:n,i-1)-ss_dist_min)...
        .*rand(nbrtocreate,1);
end

% Initial particle velocities are zero by default (should be already set in
% PSOGETINTIALPOPULATION).

% Initialize the global and local fitness to the worst possible
state.fGlobalBest = ones(itr,1)*inf; % Global best fitness score
state.fLocalBests = ones(n,1)*inf ; % Individual best fitness score

% Initialize global and local best positions
state.xGlobalBest = ones(1,nvars)*inf ;
state.xLocalBests = ones(n,nvars)*inf ;