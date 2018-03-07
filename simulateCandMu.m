function [ consumptionTS, markupTS ] = simulateCandMu( statesMatrix, C_1, policy )
%SIMULATECANDMU Summary of this function goes here
%   Detailed explanation goes here
[nDraws, tPeriods] = size(statesMatrix);

lowState = 1; % 1 is L , 2 is high
highState = 2;


consumptionTS = C_1 * ones(nDraws, tPeriods);


for t = 2:tPeriods
    
    c_L  = policy.c_L(consumptionTS(:,t-1));
    c_H  = policy.c_H(consumptionTS(:,t-1));
    state = statesMatrix(:,t);
    
    
    consumptionTS(:,t) = c_L.*( state == lowState) + c_H.*( state == highState);
    
end

markupTS = NaN(nDraws, tPeriods);

stateMat = statesMatrix(:, 2:tPeriods) ;
mu_L =  policy.mu_L(  consumptionTS(:, 1:(tPeriods-1))  );
mu_H =  policy.mu_H(  consumptionTS(:, 1:(tPeriods-1))  );

markupTS(:, 2:tPeriods) =   mu_L.*( stateMat == lowState) + mu_H.*( stateMat == highState);


end