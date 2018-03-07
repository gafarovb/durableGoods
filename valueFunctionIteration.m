

maxIterations = 200;

mu_max = cal.sigma/(cal.sigma-1);
xGrid = pPolicy_0.xGrid;
fineGrid = -1:0.01:1;
deg = size(xGrid,2)-1;
mu_L_new = zeros(deg+1,1);
V_L_new = mu_L_new;
x_L_new = mu_L_new;
mu_H_new = zeros(deg+1,1);
V_H_new = mu_H_new;
x_H_new = mu_H_new;
p = pPolicy_0;
a_L = cal.mc.values(1);
a_H = cal.mc.values(2);
transP =  cal.mc.transitionP;
cons = p.c_x(xGrid');
g_L =sol_L.cSS*ones(deg+1,1); % p.c_L(xGrid');
g_H =sol_H.cSS*ones(deg+1,1);
%%
maxV = cons(end)*(1-cal.theta)*(mu_max-1)/(1-cal.beta);
minV =  0;
restV = @(V)  min(max(V,minV),maxV);
restrict = @(mu) min(max(mu,1.0),mu_max);
restC = @(c)min(max(c,cons(1)),cons(end));
%%
h = waitbar(0,'Fixed point interation');


 theta=cal.theta;
x_c = @(c) 2 *(c - cons(1) )/(cons(end) - cons(1) )-1;

for iter = 1:maxIterations
    
    
    waitbar(iter / maxIterations,h,sprintf('Fixed point interation'))
    
    
    expMu_L = @(x) restrict (p.mu_L(x) * transP(1,1) /a_L + p.mu_H(x) * transP(1,2) /a_H);
    expMu_H = @(x) restrict (p.mu_L(x) * transP(2,1) /a_L + p.mu_H(x) * transP(2,2) /a_H);
    
    mu_L =@(x,iNode)   ( a_L *((p.c_x(x).^(-cal.gamma-(cal.sigma^-1))).*g_L(iNode).^(cal.sigma^-1) ...
        + cal.beta* cal.theta * expMu_L(x)));
    mu_H =@(x,iNode) ( a_H *((p.c_x(x).^(-cal.gamma-(cal.sigma^-1))).*g_H(iNode).^(cal.sigma^-1) ...
        + cal.beta* cal.theta * expMu_H(x)));
    expV_L = @(x) restV(p.V_L(x) * transP(1,1)  + p.V_H(x) * transP(1,2));
    expV_H = @(x) restV(p.V_L(x) * transP(2,1)  + p.V_H(x) * transP(2,2));
    
    obj_L = @(x,iNode) (mu_L(x,iNode)-1).*(p.c_x(x) - cal.theta * cons(iNode))/a_L + cal.beta * expV_L(x);
    obj_H = @(x,iNode) (mu_H(x,iNode)-1).*(p.c_x(x) - cal.theta * cons(iNode))/a_H + cal.beta * expV_H(x);
    
    xV_L_new = @(iNode) fminbnd( @(x)-obj_L(x,iNode),max(x_c(theta*cons(iNode)),-1) ,1);
    xV_H_new = @(iNode) fminbnd( @(x)-obj_H(x,iNode),max(x_c(theta*cons(iNode)),-1),1);
    
    for iNode=1:(deg+1)
        [x_L_new(iNode),V_L_new(iNode)] = xV_L_new(iNode);
        mu_L_new(iNode) =mu_L( x_L_new(iNode),iNode);
        
        V_L_new(iNode) = restV(-V_L_new(iNode));
        [x_H_new(iNode),V_H_new(iNode)] = xV_H_new(iNode);
        mu_H_new(iNode) = mu_L( x_H_new(iNode),iNode);
        V_H_new(iNode) = restV(-V_H_new(iNode));
        
    end
    C_L_new = p.c_x(x_L_new);
    C_H_new = p.c_x(x_H_new);
    
    
    policy_new = struct('mu_L',mu_L_new,'c_L',C_L_new,'V_L',V_L_new,...
        'mu_H',mu_H_new,'c_H',C_H_new,'V_H',V_H_new,'xGrid',xGrid','c_x',p.c_x);
    
    [ p_New ] = convertPolynomial_step2( policy_new );
    lambda= 1;
    policy = policy_new;
    
    [ p ] = blendPolicy(lambda, p_New,p );
    g_L =  C_L_new ;
    g_H =  C_H_new ;
    
    
    
    %    plot(p.c_x(fineGrid'),p.mu_L(fineGrid'),p.c_x(fineGrid'),p.mu_H(fineGrid'));
  %  plot(p.c_x(fineGrid'),p.c_L(fineGrid'),p.c_x(fineGrid'),p.c_H(fineGrid'));
%
    
 % plot(p.c_x(fineGrid'),p_New.c_L(fineGrid'),p.c_x(fineGrid'),p_New.c_H(fineGrid'));
  plot(p.c_x(fineGrid'),p.V_L(fineGrid'),p.c_x(fineGrid'),p.V_H(fineGrid'));  
end
delete(h)

%%

plot(p.c_x(fineGrid'),p_New.c_L(fineGrid'),p.c_x(fineGrid'),p_New.c_H(fineGrid'));
%%

  plot(p.c_x(fineGrid'),p.mu_L(fineGrid'),p.c_x(fineGrid'),p.mu_H(fineGrid'));  
%%

  plot(p.c_x(fineGrid'),p.V_L(fineGrid'),p.c_x(fineGrid'),p.V_H(fineGrid'));  

%%
figure(1)
plot(p.c_x(fineGrid'),pPolicy_0.c_L(fineGrid'),p.c_x(fineGrid'),pPolicy_0.c_H(fineGrid'));
figure(2)

plot(p.c_x(fineGrid'),pPolicy_0.mu_L(fineGrid'),p.c_x(fineGrid'),pPolicy_0.mu_H(fineGrid'));
figure(3)

plot(p.c_x(fineGrid'),pPolicy_0.V_L(fineGrid'),p.c_x(fineGrid'),pPolicy_0.V_H(fineGrid'));


%%

figure(1)
plot(p.c_x(fineGrid'),p_New.c_L(fineGrid'),p.c_x(fineGrid'),p_New.c_H(fineGrid'));
figure(2)
plot(p.c_x(fineGrid'),p_New.mu_L(fineGrid'),p.c_x(fineGrid'),p_New.mu_H(fineGrid'));
figure(3)
plot(p.c_x(fineGrid'),p_New.V_L(fineGrid'),p.c_x(fineGrid'),p_New.V_H(fineGrid'));
%%

policy_c = switchToC(p);

%policy_c_sim = computeInitialValueFunction(policy_c,cal);
