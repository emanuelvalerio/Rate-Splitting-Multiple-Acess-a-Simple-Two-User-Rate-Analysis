function [rate] =  rateCalculate(a,b,c,d,t)
    rate = log2(a*c + (a*d + b*c)*t + b*d*t^2);
end