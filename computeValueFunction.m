function [ value ] = computeValueFunction( C_T, mu_T ,state_2, cal )
%COMPUTEVALUEFUNCTION Summary of this function goes here
%   Detailed explanation goes here
a_2 = cal.mc.values(state_2);
[nDraws, tPeriods] = size(C_T);

lagC_T = 0 *C_T;
lagC_T(:,2:(tPeriods)) = C_T(:,1:(tPeriods-1));
payoffs = (mu_T-1) /a_2 .* (C_T  - cal.theta * lagC_T);
payoffShort =payoffs(:,2:end );



NPV = cal.beta.^  (0:(tPeriods-2))';
NPV(end) = NPV(end)/(1-cal.beta);
discountedPayoffs = payoffShort*NPV;

value = mean(discountedPayoffs,1);
errorEstimate =3*std(discountedPayoffs)/sqrt(nDraws);

 end

