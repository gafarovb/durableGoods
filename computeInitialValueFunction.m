function policy_0 = computeInitialValueFunction(policy_0,cal)
%%
MarkovChain= cal.mc ;

nDraws = 4000;
tPeriods = 200;
cGrid = policy_0.cGrid;
nGrid = size(cGrid,1);
V_L = NaN(nGrid,1);
V_H = NaN(nGrid,1);
h = waitbar(0,'Computing myopic value function');

for iC = 1:nGrid
    
    waitbar(iC / nGrid,h )
    lowState = 1; % 1 is L , 2 is high
    statesMatrix_L  = simulateMC( MarkovChain,  lowState, nDraws, tPeriods );
    
    highState = 2;
    statesMatrix_H  = simulateMC( MarkovChain, highState, nDraws, tPeriods );
    
    C_0 = policy_0.cGrid(iC);
    
    [  C_L, mu_L ] = simulateCandMu( statesMatrix_L, C_0, policy_0 );
    V_L(iC) = computeValueFunction( C_L, mu_L ,lowState, cal );
    
    
    [  C_H, mu_H ] = simulateCandMu( statesMatrix_H, C_0, policy_0 );
    V_H(iC) = computeValueFunction( C_H, mu_H ,highState, cal );
    
end
 delete(h)

 policy_0.('V_L')=V_L;
 policy_0.('V_H')=V_H;

end



 
%plot(cGrid,V_L,cGrid,V_H)
