function fc = calculateFc(h1,h2,P1,P2,rho,gama)

        %   Detailed explanation goes here

        y1 = sqrt(1+(norm(h1)^2)*rho*P1);
        y2 = sqrt(1+(norm(h2)^2)*rho*P2);
        H1 = h1/y1;
        H2 = h2/y2;
        
        alphaMatrix = [H1'; H2'] * [H1 H2] ;
        a = (alphaMatrix(2,2) - abs(alphaMatrix(1,2))) ;
        b = (alphaMatrix(1,1) - abs(alphaMatrix(1,2)));
        muMatrix = (1 / (alphaMatrix(1,1) + alphaMatrix(2,2) - 2*abs(alphaMatrix(1,2))))*[ a;b ];
        lambda = (alphaMatrix(1,1)*alphaMatrix(2,2) - abs(alphaMatrix(1,2))^2) / ((alphaMatrix(1,1) + alphaMatrix(2,2) - 2*abs(alphaMatrix(1,2))));
        
        fc = (1/sqrt(lambda))*(muMatrix(1)*H1 + muMatrix(2)*H2*exp(-j*angle(alphaMatrix(1,2))));

end