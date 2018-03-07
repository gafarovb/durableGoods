 
cal = struct('theta',0.5,...
    'beta',0.99,...
    'sigma',5,...
    'a',1,...
    'gamma',2);
cSS = initialGuessCSS( cal )
sol = deterministicSolution( cal ,cSS)
 

%%
  sol = sol.findFixPoint;
disp('Saving the results')
 
%%
close all
sol.plotConsumption;

sol.plotMarkup;
sol.steadyState
 
