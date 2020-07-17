function mexc_sig(k)
% Syntax: f = mexc_sig(k)
% 1:20 PM 15/08/97
% defines modulation signal for exciter control
% global exc_sig n_exc
% f=0; %dummy variable
global g

if g.exc.n_exc~=0
    % exc_sig(:,k)=zeros(n_exc,1);
    % exc_sig(1,k)=0.1;
    %end
    if g.sys.t(k)<=0.5
        g.exc.exc_sig(1,k) = 0;
    else
        %exc_sig(:,k) = zeros(n_exc,1);
        g.exc.exc_sig(1,k) = 0.01;
    end
end
return