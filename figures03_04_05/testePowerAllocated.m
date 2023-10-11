function [t_optimal,rate] = testePowerAllocated(Nt,P,h1,h2,rho,Gamma,gama,t0)
     
     t = linspace(0,1,100);
     Pk = t0*P;
     Pc = (1-t0)*P;
     g1 = (norm(h1)^2)*rho;
     g2 = (norm(h2)^2)*rho;
     g = [g1;g2];
     [WF] = waterPouring(Pk,g,Nt);
     P1 = WF(1);
     P2 = WF(2);
  
     mu = (1/2)*(t0*P + (1/g1) + 1/g2);
     bestRate = zeros(1,length(t));

     if mu <=  1/g2

            for y = 1:length(t)
                 P1 = t(y)*P;
                 Pc = (1-t(y))*P;
                 fc = calculateFc(h1,h2,P1,P2,rho,gama);
                 [rate] =  generalRateFunction(h1,h2,P1,P2,Pc,fc,rho);
                 bestRate(y) = rate;               
            end
           [rate,maxIndex] = max(bestRate);
           t_optimal = t(maxIndex);
 
     end

     if mu > 1/g2
           fc = calculateFc(h1,h2,P1,P2,rho,gama);
           if (t0 > 0 && t0<=1 && rho > 0 && P1==P2 && P2==t0*P*0.5)  % High SNR
              e = (((norm(h2)^2)*rho)/4) - ((abs(h2'*fc)^2)/2);
              f = ((abs(h2'*fc)^2)/2);
              t_optimal = min((-f/(2*e)),1);
              rate = (log2(rho) + 2*log2(P)+log2(e*t_optimal^2 + f*t_optimal));
           else
             [a,b,c,d] = rateParameters(h2,rho,Gamma,P,fc);
             t_optimal = min(((-a/(2*b))-c/(2*d)),1);         
             [rate] =  rateCalculate(a,b,c,d,t_optimal); 
           end
     end
end
 
