% Test script for scaling factor

PFA = 10^(-3);
RefWindow = 32;
k = floor((3*RefWindow)/4);
PFA_range_upper_limit = PFA*1.001; 
PFA_range_lower_limit = PFA;

for a = 1:0.001:20
    numerator = factorial(RefWindow)*gamma(a+RefWindow-k+1);
    denominator = factorial(RefWindow-k)*gamma(a+RefWindow+1);
    eqn = numerator/denominator;
    if ((eqn <= PFA_range_upper_limit) && (eqn >= PFA_range_lower_limit))
        a
        break       
    end
end
