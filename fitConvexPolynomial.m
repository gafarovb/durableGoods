function [ a,p_a ] = fitConvexPolynomial( Y,A, deg,x )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nConstraints = size(Y,1);
step = 2 /( nConstraints-1 );



if nargin<4
    x = -1:step:1;
end
nConstraints = size(A,1);
X =(x'.^(0:1:deg))';
b =  zeros(nConstraints,1);
H= 2 * X*(X');
H=(H+H')/2;
f = - (2  * X *Y);
options = optimoptions('quadprog',...
    'Algorithm','interior-point-convex','Display','off');
a =  quadprog(H,f,A,b,[],[],[],[],[],options);
p_a = @(x) (x.^(0:1:deg))*a;

end

