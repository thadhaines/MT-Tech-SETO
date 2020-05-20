function [Gov] = funGov_TGOV1(Gov,Gen,ii,dt,Flag)

for jj = 1:length(Gov)
    if Flag == 0 % Initialize.

        % Parameters.
        if 1==1
            Gov(jj).k(1) = 1-Gov(jj).T2/Gov(jj).T3;
            Gov(jj).k(2) = Gov(jj).k(1)/Gov(jj).T3;
            Gov(jj).k(3) = Gov(jj).T2/Gov(jj).T3;
        end

        Gov(jj).x1(ii,1) = Gov(jj).Pref(ii,1);
        Gov(jj).x2(ii,1) = Gov(jj).k(1) * Gov(jj).Pref(ii,1);

    elseif Flag == 1; % Pmech.

        Gov(jj).Pm(ii,1) = Gov(jj).x2(ii,1) + ...
                           Gov(jj).k(3) * Gov(jj).x1(ii,1) - ...
                           Gov(jj).Dt * (Gen(jj).Wr(ii,1)-1);

    elseif Flag == 2 % Update states.

        dw = Gen(jj).Wr(ii,1)-1;
        d = Gov(jj).Pref(ii,1)-dw/Gov(jj).R;

%         if d >= Gov(jj).Vmax
%             d = Gov(jj).Vmax;
%         elseif d<=Gov(jj).Vmin
%             d = Gov(jj).Vmin;
%         end

        Gov(jj).x1dot(ii,1) = (-Gov(jj).x1(ii,1) + d)/Gov(jj).T1;
        Gov(jj).x2dot(ii,1) = -Gov(jj).x2(ii,1)/Gov(jj).T3 + Gov(jj).k(2)*Gov(jj).x1(ii,1);

        Gov(jj).x1(ii+1,1) = Gov(jj).x1(ii,1) + dt*(1.5*Gov(jj).x1dot(ii,1) - 0.5*Gov(jj).x1dot(ii-1,1));
        Gov(jj).x2(ii+1,1) = Gov(jj).x2(ii,1) + dt*(1.5*Gov(jj).x2dot(ii,1) - 0.5*Gov(jj).x2dot(ii-1,1));
        
        clear d dw
    else
        error('Invalid Flag')
    end
end

end