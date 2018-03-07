
function out =  equationInputsMex1(f,g,f1,C_ss,sigma,theta,beta)
out =  (g.^(-sigma^-1)) * (C_ss^((sigma^-1) -1)) + beta * theta .* fg - f ;
end

 
 
