function [ cSS , cSScommitement , cSSub ] = initialGuessCSS( cal )
%INITIALGUESSCSS Summary of this function goes here
%   Detailed explanation goes here
shift = 1;

cSScommitement = ( cal.a / (1- cal.beta* cal.theta)*(cal.sigma - 1 )/( cal.sigma))^(1/cal.gamma) ;
cSSub =( cal.a / (1- cal.beta* cal.theta))^(1/cal.gamma);  
cSS = cSScommitement + shift* cal.theta*(  cSSub - cSScommitement);

end

