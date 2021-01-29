% Squashed 2D grid indexing to 1D indexing
% Rest of simulation uses this convention, which is kinda standard
function idx = indexf(i,j,syst)
idx = i+syst.n*(j-1);
end