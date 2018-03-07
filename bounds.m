 
cal_L = struct('theta',0.9,...
    'beta',0.95,...
    'sigma',2,...
    'a',1,...
    'gamma',1);
cal_U = cal_L;
cal_U.a = 1.1;


MarkovChain = struct('transitionP',[0.6 0.4;0.4 0.6],...
    'values', [cal_L.a, cal_U.a]);


cal  = cal_L;
cal.('mc') = MarkovChain;
 %%
[cSS_L,lb,~] = initialGuessCSS( cal_L );
[cSS_U,~,ub] = initialGuessCSS( cal_U );


boundsX = struct( 'low', lb, ...
            'up', ub);

sol_L = deterministicSolution( cal_L ,cSS_L,boundsX);
sol_L = sol_L.findFixPoint;

sol_H = deterministicSolution( cal_U  ,sol_L.cSS, boundsX);
sol_H = sol_H.findFixPoint;

 policy_0 = struct  ('cGrid',sol_L.fullGrid.x,...
    'c_L',sol_L.consumptionPolicy,...
    'mu_L',sol_L.markupPolicy,...
    'c_H',sol_H.consumptionPolicy,...
    'mu_H',sol_H.markupPolicy);
%%

policy_0 = computeInitialValueFunction(policy_0,cal);

%%

degree = 30;
pPolicy_0 = convertPolynomial(policy_0, degree);