function [powerAlocated , regime,t_optimal,idx,rate] = newPowerAllocated(Nt,P,h1,h2,rho,Gamma,gama,t)
     g1 = (norm(h1)^2)*rho;
     g2 = (norm(h2)^2)*rho;
     g = [g1;g2];
     Pk = t*P;
     Pc = (1-t)*P;
     [WF] = waterPouring(Pk,g,Nt);
     P1 = WF(1);
     P2 = WF(2);
     tolerance = 10^-5; %convergence tolerance;
     mu = (1/2)*(t*P + (1/g1) + 1/g2);
     if mu <=  1/g2

                if abs(P1-0)<tolerance
                    regime = "Multicast";
                    idx = 1;
                elseif abs(Pc-0)<tolerance
                     regime = "OMA";
                     idx = 2;
                else
                    regime = "NOMA";
                    idx = 3;
                end

     
     else
         if abs(Pc-0)<tolerance
             regime = "SDMA";
             idx=4;
         else
             regime = "RSMA";
             idx=5;
         end

     end
     powerAlocated(1) = P1;
     powerAlocated(2) = P2;
     powerAlocated(3) = Pc;
     t_optimal = 0;
     fc = calculateFc(h1,h2,P1,P2,rho,gama);
     rate = generalRateFunction(h1,h2,P1,P2,Pc,fc,rho);
end