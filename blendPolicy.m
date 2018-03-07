function [ pblend ] = blendPolicy(lambda, p1,p2 )
%BLENDPOLICY Summary of this function goes here
%   Detailed explanation goes here
pblend = p1;

funcList ={'c_H','c_L','V_H','V_L','mu_H','mu_L'};

for ix=1:6
    pblend.(funcList{ix})=@(x)(lambda*(p1.(funcList{ix})(x))+(1-lambda)*(p2.(funcList{ix})(x)));
end

 
[ pblend ] = convertPolynomialAfterBlend( pblend );

if lambda==1
    pblend = p1;
end
end

