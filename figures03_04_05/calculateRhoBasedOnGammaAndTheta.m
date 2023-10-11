% auxiliar functions
function [rho,h1,h2] = calculateRhoBasedOnGammaAndTheta(gamma, theta)
    h1 = (1/sqrt(2)).*([1, 1])';
    h2 = (gamma/sqrt(2)).*(([1, exp(1j * theta)]))';
    h_1 = h1./norm(h1);
    h_2 = h2./norm(h2);
    rho = 1 - (abs(h_1' * h_2)^2);
end