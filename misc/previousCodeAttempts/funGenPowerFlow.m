function [delta,Pg] = funGenPowerFlow(Y,Pg,V,delta0,H,PslackTol)
% Iteratively solve the power flow until the Pg balances with the
% intertias.  
% Inputs:
%   Y = Y bus matrix for generators
%   Pg = vector of injected real power from each gen.
%        Pg(1) is the swing gen.
%   V = vector of bus terminal voltage magnitudes
%   delta0 = initial angles for V
%   H = vector of generator inertia constants
%   PslackTol = tolerance error for Pslack
% Outputs:
%   delta = vector of solved angles
%   Pg = vector of solved injected powers

%% Settings, error checks, and initialize
tol = 1e-5; %power mismatch tolerance
MaxItPF = 100; %Max iterations for PF
MaxItPg = 20; %Max iterations for Pg
M = length(Pg);
if M~=length(delta0); error(' '); end
if max(abs(size(Y)-[M M]))>eps; error(' '); end
if length(H)~=M; error(' '); end

% Initialize
er = ones(M-1,1);
delta = delta0;
J = zeros(M-1,M-1);
dPslack = 1;
HT = sum(H);

%% Solve PF iteratively
kItPg = 0; %counter
while dPslack>PslackTol
    kItPg = kItPg + 1;
    
    %Solve PF
    kItPF = 0;
    while max(abs(er))>tol
        kItPF = kItPF + 1;
        %Build er and J
        for k=2:M
            er(k-1) = Pg(k);
            for n=1:M
                er(k-1) = er(k-1) - V(k)*abs(Y(k,n))*V(n)*cos(delta(k)-delta(n)-angle(Y(k,n)));
            end
            for m=2:M
                if m==k
                    J(k-1,k-1) = 0;
                    for n=1:M
                        if k~=n
                            J(k-1,k-1) = J(k-1,k-1) - V(k)*abs(Y(k,n))*sin(delta(k)-delta(n)-angle(Y(k,n)));
                        end
                    end
                else
                    J(k-1,m-1) = V(k)*abs(Y(k,m))*V(m)*sin(delta(k)-delta(m)-angle(Y(k,n)));
                end
            end
        end
        clear k n m

        %solve
        Ddelta = J\er;
        delta(2:end) = delta(2:end) + Ddelta;

        %check for convergence
        if kItPF>MaxItPF; error('No PF convergence'); end
    end

    %Solve for slack power error and adjust Pg
    dPslack = Pg(1);
    for n=1:M
        dPslack = dPslack - V(1)*abs(Y(1,n))*V(n)*cos(delta(1)-delta(n)-angle(Y(1,n)));
    end
    for k=1:M
        Pg(k) = Pg(k) - H(k)*dPslack/HT;
    end
    clear n k
    
    %converging?
    if kItPg>MaxItPg; error('No Pg convergence'); end
end

end

