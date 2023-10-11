function [rateMulticast,rateOMA,rateNOMA,rateSDMA,rateRSMA] = rateSchemes(Nt,gamadb,theta,matrixT,P)
        rateMulticast = zeros(length(gamadb),length(theta));
        rateOMA = zeros(length(gamadb),length(theta));
        rateNOMA = zeros(length(gamadb),length(theta));
        rateSDMA = zeros(length(gamadb),length(theta));
        rateRSMA = zeros(length(gamadb),length(theta));

        for i = 1:length(gamadb)
            for j = 1:length(theta)
                gama = (10.^(gamadb(i)/20));
                [rho,h1,h2] = calculateRhoBasedOnGammaAndTheta(gama,theta(j)); % Calculate ρ
                g1 = (norm(h1)^2)*rho;
                g2 = (norm(h2)^2)*rho;
                g = [g1;g2];
                P1 = 0;
                P2 = 0;
                Pc = P; 
                fc = calculateFc(h1,h2,P1,P2,rho,gama);
                rateMulticast(i,j) =  generalRateFunction(h1,h2,P1,P2,Pc,fc,rho);

                P1 =P;
                P2 = 0;
                Pc = 0;
                fc = calculateFc(h1,h2,P1,P2,rho,gama);
                rateOMA(i,j) =  generalRateFunction(h1,h2,P1,P2,Pc,fc,rho);

                P1 = matrixT(i,j)*P;
                P2 = 0;
                Pc = (1-matrixT(i,j))*P;
                fc = calculateFc(h1,h2,P1,P2,rho,gama);
                rateNOMA(i,j) =  generalRateFunction(h1,h2,P1,P2,Pc,fc,rho);

                Pc = 0;
                [power] = waterPouring(matrixT(i,j)*P,g,Nt);
                P1 = power(1);
                P2 = power(2);
                rateSDMA(i,j) =  generalRateFunction(h1,h2,P1,P2,Pc,fc,rho);

                Pc = (1 - matrixT(i,j))*P;
                [power] = waterPouring(matrixT(i,j)*P,g,Nt);
                P1 = power(1);
                P2 = power(2);
                rateRSMA(i,j) =  generalRateFunction(h1,h2,P1,P2,Pc,fc,rho);
            end
        end
end