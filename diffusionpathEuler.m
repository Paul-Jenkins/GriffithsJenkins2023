% Function to implement the Euler-Maruyama method. Input: parameters
% theta_A, theta_B, P_A, P_B, rho, R (number of replicates), N (number of
% timesteps), dt (stepsize), x_0. Returns a table containing, for each
% replicate: rhohat,rhoerror,minx (minimum of X_{ij}(t) across
% (i,j,t),rhohatN (numerator of rho-hat),rhoerrorN (numerator of
% rho-error),I,Iinf (indicator for I_T = infinity).
function results = diffusionpathEuler(tA,tB,PA,PB,rho,R,N,dt,x0)

tol = 1e-10; % Error tolerance in some 'assert' checks
%T = dt*(N-1);
K = length(PA);
L = length(PB);
H = K*L;
% x = [x11 x12; x21 x22]

driftmu(x0,tA,tB,PA,PB,rho); % Populate static parameters inside function

rhohatN = zeros(R,1); % Numerator of rho_MLE
I = rhohatN;          % Information  (= rhohatD)
I2 = I;               % Information computed a different way
Iinf = I;             % Indicator for whether I(r) = Inf
rhohat = rhohatN;
rhoerrorN = rhohatN;  % rhoerror*I
rhoerror = rhohatN;   % rhohat - rho, as measured by a martingale term
minx = rhohatN;       % min x_ij(t) encountered

for r=1:R
    r
    x = zeros(K,L,N);
    x(:,:,1) = x0;
    dW = sqrt(dt)*randn(H,N-1);

    for n=2:N
        [recdrift,mutdrift] = driftmu(x(:,:,n-1));
        diffusionterm = reshape(diffusionsigma(x(:,:,n-1))*dW(:,n-1),[K,L]);
        x(:,:,n) = x(:,:,n-1) + (recdrift+mutdrift)*dt + diffusionterm;
        
        assert(all(size(recdrift) == [K,L]));
        assert(all(size(mutdrift) == [K,L]));
        assert(all(size(x(:,:,n-1)) == [K,L]));
        assert(all(size(x(:,:,n)) == [K,L]));
        if (any(x(:,:,n) < 0,'all') || any(x(:,:,n) > 1,'all'))
            M = x(:,:,n);
            M(M<0) = 0;
            M(M>1) = 1;
            x(:,:,n) = M/sum(M,'all');
        end
        assert(all(size(x(:,:,n)) == [K,L]));
    
        for k=1:K
            for l=1:L
                xkdot = sum(x(k,:,n-1));
                xdotl = sum(x(:,l,n-1));
                rhohatN(r) = rhohatN(r) + (xkdot*xdotl/x(k,l,n-1))*(x(k,l,n)-x(k,l,n-1) - mutdrift(k,l)*dt);
                if (all(x(:,:,n-1) > 0,'all'))
                    I(r) = I(r) + xkdot^2*xdotl^2*dt/x(k,l,n-1);
                    I2(r) = I2(r) + (x(k,l,n-1) - xkdot*xdotl)^2*dt/x(k,l,n-1);
                    rhoerrorN(r) = rhoerrorN(r) + xkdot*xdotl/x(k,l,n-1)*diffusionterm(k,l);
                else
                    I(r) = Inf;
                    I2(r) = Inf;
                    Iinf(r) = 1;
                end
            end
        end
    end
    minx(r) = min(x,[],'all');
    I(r) = I(r) - (N-1)*dt;
    rhohat(r) = rhohatN(r)/I(r);
    rhoerror(r) = rhoerrorN(r)/I(r);
end

% Icheck = abs(I-I2);
% assert(all(Icheck(not(Iinf)) < tol));

results = table(rhohat,rhoerror,minx,rhohatN,rhoerrorN,I,Iinf);

end