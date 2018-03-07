function [ X ] = simulateMC( MarkovChain, initialState, nDraws, tPeriods )
%SIMULATEMC Summary of this function goes here
%   Detailed explanation goes here


X = ones(nDraws,tPeriods);
X(:,1) =   NaN(nDraws,1);

X(:,2) = initialState * ones(nDraws,1);
 
eps = rand(nDraws,tPeriods);
 
eps( 1:(nDraws/2),:) = 1-eps((nDraws/2+1):end,:);


for tx = 3:tPeriods
    
for ix =1:nDraws
    prevState = X(ix,tx-1);
    altState = 1 *(1~=prevState) + 2 *(2~=prevState);
    pStay = MarkovChain.transitionP(prevState,prevState);
    stayed = eps(ix,tx)<pStay;
    
    X(ix,tx) =  prevState *( stayed) + altState *( ~stayed);
    
end

end


 

 