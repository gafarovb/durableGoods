function [ polynomPolicy ] = convertPolynomial_step2( policy )
%CONVERTPOLYNOMIAL Summary of this function goes here
%   Detailed explanation goes here

xGrid = policy.xGrid;
nConstraints = size(xGrid,1);
step = 2 /( nConstraints-1 );
x = xGrid';
deg = nConstraints-1;
 

convex = zeros(nConstraints, deg+1);
decreasing = convex;
for iDeg = 0:deg
   decreasing(:,iDeg+1) = ( iDeg * x.^max([iDeg-1,0]))';
   convex(:,iDeg+1) = -  (max([iDeg-1,0])* max([iDeg,0])* x.^max([iDeg-2,0]))';
end
concave = -convex;
increasing = - decreasing;

polynomPolicy = policy;
xGridTran = xGrid';
% fit V
A_V =[decreasing ;convex]; %[decreasing;convex];
[~,polynomPolicy.V_L] = fitConvexPolynomial( policy.V_L, A_V, deg,xGridTran );
[~,polynomPolicy.V_H] = fitConvexPolynomial( policy.V_H, A_V, deg,xGridTran );
% fit Mu
 
A_mu =[decreasing;convex]; %decreasing;% [decreasing;convex];

[~,polynomPolicy.mu_L] = fitConvexPolynomial(policy.mu_L,A_mu, deg,xGridTran ); 
[~,polynomPolicy.mu_H ]= fitConvexPolynomial(policy.mu_H,A_mu, deg,xGridTran );


% fit C

A_c =  [increasing;concave ];%concave
 
[~,polynomPolicy.c_L] = fitConvexPolynomial(policy.c_L,A_c, deg,xGridTran );
[~,polynomPolicy.c_H ]= fitConvexPolynomial(policy.c_H,A_c, deg,xGridTran );

 
end

