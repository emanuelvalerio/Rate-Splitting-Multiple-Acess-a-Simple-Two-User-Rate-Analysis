function [rate] =  generalRateFunction(h1,h2,P1,P2,Pc,fc,rho)
     rate = log2(1+((norm(h1)^2)*rho*P1))+log2(1+((norm(h2)^2)*rho*P2)+(abs(h2'*sqrt(Pc)*fc)^2));
end