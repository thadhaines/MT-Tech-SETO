function pwrmod_indx(bus)
% syntax: f = pwrmod_indx
% determines the relationship between pwrmod and nc loads
% checks for pwrmod
% determines number of modulated injections
% f is a dummy variable
% Trudnowski, Feb. 2015

% global pwrmod_con n_pwrmod  pwrmod_idx 
% global bus_int

global load_con  

global g

g.pwr.n_pwrmod = 0;
g.pwr.pwrmod_idx = [];
if ~isempty(g.pwr.pwrmod_con)
    if isempty(load_con)
        error('you must have the power modulation bus declared as a consant-power non-conforming load using load_con'); 
    end
    g.pwr.n_pwrmod = length(g.pwr.pwrmod_con(:,1));
    g.pwr.pwrmod_idx = zeros(g.pwr.n_pwrmod,1);
    for j = 1:g.pwr.n_pwrmod
       index = find(g.pwr.pwrmod_con(j,1)==load_con(:,1));
       if ~isempty(index)
           if abs(sum(load_con(index,2:end)')-2)>1e-8; 
               error('pwrmod buses must be defined as 100% constant power or 100% constant current in load_con'); 
           end
           if (load_con(index,2)==1 && load_con(index,3)==1) || (load_con(index,4)==1 && load_con(index,5)==1)
               g.pwr.pwrmod_idx(j) = index;
           else
               error('pwrmod buses must be defined as 100% constant power or 100% constant current in load_con'); 
           end
           kk = g.sys.bus_int(g.pwr.pwrmod_con(:,1));
           if any(abs(bus(kk,10)-2))
               error('power modulation buses must be declared type 2 (PV) buses'); 
           end
           if max(abs(bus(kk,6)))>1e-10 || max(abs(bus(kk,7)))>1e-10
               error('power modulation buses must have zero load');
           end
           clear kk
       else
          error('you must have the power modulation bus declared as a 100% consant-power or 100% constant current non-conforming load')
       end
    end
end
       
    