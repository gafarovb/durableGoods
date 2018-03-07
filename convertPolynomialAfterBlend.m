function [ polynomPolicy ] = convertPolynomialAfterBlend( policy )
%CONVERTPOLYNOMIAL Summary of this function goes here
%   Detailed explanation goes here

xGrid = policy.xGrid;
nConstraints = size(xGrid,1);
 
x = xGrid';
deg = nConstraints-1;


cons = policy.c_x(xGrid');
restC = @(c)min(max(c,cons(1)),cons(end));

convex = zeros(nConstraints, deg+1);
decreasing = convex;
for iDeg = 0:deg
   decreasing(:,iDeg+1) = ( iDeg * x.^max([iDeg-1,0]))';
   convex(:,iDeg+1) = -  (max([iDeg-1,0])* max([iDeg,0])* x.^max([iDeg-2,0]))';
end
concave = -convex;
increasing = - decreasing;

polynomPolicy = policy;

% fit V
A_mu = [ ];
Y_vL = policy.V_L(xGrid);
Y_vH = policy.V_H(xGrid);



[~,polynomPolicy.V_L] = fitConvexPolynomial( Y_vL, A_mu, deg );
[~,polynomPolicy.V_H] = fitConvexPolynomial( Y_vH, A_mu, deg );
% fit Mu
Y_muL = policy.mu_L(xGrid);

[~,polynomPolicy.mu_L] = fitConvexPolynomial(Y_muL,A_mu, deg );

Y_muH = policy.mu_H(xGrid);

[~,polynomPolicy.mu_H ]= fitConvexPolynomial(Y_muH,A_mu, deg );


% fit C

A_c =  [];
Y_CL =restC( policy.c_L(xGrid));

[~,polynomPolicy.c_L] = fitConvexPolynomial(Y_CL,A_c, deg );

Y_CH = restC(policy.c_H(xGrid));
[~,polynomPolicy.c_H ]= fitConvexPolynomial(Y_CH,A_c, deg );
 
end

