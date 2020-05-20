function [Gov] = funInitialize4AdamsBashTGOV1(Gov,ii)

for jj = 1:length(Gov)
    
    Gov(jj).Pref(ii+1,1) = Gov(jj).Pref(ii,1);
    
    Gov(jj).x1dot(ii,1) = 0;
    Gov(jj).x2dot(ii,1) = 0;

    Gov(jj).x1(ii+1,1) = Gov(jj).x1(ii,1);
    Gov(jj).x2(ii+1,1) = Gov(jj).x2(ii,1);
end

end