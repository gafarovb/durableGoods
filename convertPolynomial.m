function [ polynomPolicy ] = convertPolynomial( policy,deg )
%CONVERTPOLYNOMIAL Summary of this function goes here
%   Detailed explanation goes here
cGrid = policy.cGrid;
nConstraints = size(cGrid,1);
step = 2 /( nConstraints-1 );
x = -1:step:1;

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
A_mu = [decreasing;convex];
[~,polynomPolicy.V_L] = fitConvexPolynomial( policy.V_L, A_mu, deg );
[~,polynomPolicy.V_H] = fitConvexPolynomial( policy.V_H, A_mu, deg );
% fit Mu
Y_muL = policy.mu_L(cGrid);

[~,polynomPolicy.mu_L] = fitConvexPolynomial(Y_muL,A_mu, deg );

Y_muH = policy.mu_H(cGrid);

[~,polynomPolicy.mu_H ]= fitConvexPolynomial(Y_muH,A_mu, deg );


% fit C

A_c =  [increasing;concave];
Y_CL = policy.c_L(cGrid);

[~,polynomPolicy.c_L] = fitConvexPolynomial(Y_CL,A_c, deg );

Y_CH = policy.c_H(cGrid);
[~,polynomPolicy.c_H ]= fitConvexPolynomial(Y_CH,A_c, deg );

polynomPolicy.('c_x') = @(x) (cGrid(1) + (cGrid(end)-cGrid(1))* (x+1));
polynomPolicy.('xGrid') = - cos((0:deg)*pi/deg);

end

