function [ policy_out ] = switchToC(policy )
%BLENDPOLICY Summary of this function goes here
%   Detailed explanation goes here
policy_out = policy;

c_x=policy.c_x;
policy_out.('x_c')=@(c)-1 + 2*(c -c_x(-1));
funcList ={'c_H','c_L','V_H','V_L','mu_H','mu_L'};

for ix=1:6
    policy_out.(funcList{ix})=@(c)policy.(funcList{ix})(policy_out.x_c(c) );
end

 policy_out.cGrid=c_x(policy_out.xGrid);
 
end

