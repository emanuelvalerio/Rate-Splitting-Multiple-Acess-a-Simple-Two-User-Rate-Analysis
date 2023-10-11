function  [a,b,c,d] = rateParameters(h2,rho,Gamma,P,fc)
      b = rho*P*0.5;
      a = 1 + ((Gamma * b) / P);
      d = ((norm(h2)^2)*rho*P*0.5) - ((abs(h2'*fc)^2) * P);
      c = 1 - (Gamma*d / P) + ((abs(h2'*fc)^2) * (P - Gamma));
end