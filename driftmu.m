function [recdrift,mutdrift] = driftmu(x,varargin) % x = [x11 x12; x21 x22]
% Function to compute the drift mu_ij on p13 of the paper. Input: a 2x2
% matrix x and (optional) arguments rho, theta_A, theta_B, P_A, P_B.
% Returns the recombination and mutation components of mu_ij separately.
persistent tA tB PA PB rho;

if (nargin > 1)
    tA = varargin{1}
    tB = varargin{2}
    PA = varargin{3}
    PB = varargin{4}
    rho = varargin{5}
end

x11 = x(1,1);
x21 = x(2,1);
x12 = x(1,2);
x22 = x(2,2);

recdrift = rho*[(x11+x12)*(x11+x21)-x11 (x11+x12)*(x12+x22)-x12;
                 (x21+x22)*(x11+x21)-x21 (x21+x22)*(x12+x22)-x22];

mutdrift = tA/2*[(x11+x21)*PA(1) - x11 (x12+x22)*PA(1) - x12;
            (x11+x21)*PA(2) - x21 (x12+x22)*PA(2) - x22] ...
        + tB/2*[(x11+x12)*PB(1) - x11 (x11+x12)*PB(2) - x12;
                (x21+x22)*PB(1) - x21 (x21+x22)*PB(2) - x22];

assert(all(size(recdrift) == [2,2]));
assert(all(size(mutdrift) == [2,2]));
end