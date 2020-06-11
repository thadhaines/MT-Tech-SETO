function f = pwrm_indx(bus)
% syntax: f = pwrm_indx
% determines the relationship between pwrmod and nc loads
% checks for pwrmod
% determines number of modulated buses
% f is a dummy variable
% Trudnowski, Feb. 2015
f = 0;
global pwrmod_con load_con  n_pwrmod  pwrmod_idx bus_int
n_pwrmod = 0;
pwrmod_idx = [];
if ~isempty(pwrmod_con)
    if isempty(load_con); error('you must have the power modulation bus declared as a consant-power non-conforming load using load_con'); end
    n_pwrmod = length(pwrmod_con(:,1));
    pwrmod_idx = zeros(n_pwrmod,1);
    for j = 1:n_pwrmod
       index = find(pwrmod_con(j,1)==load_con(:,1));
       if ~isempty(index)
           if abs(sum(load_con(index,2:end)')-2)>1e-8; error('pwrmod buses must be defined as 100% constant power or 100% constant current in load_con'); end
           if (load_con(index,2)==1 && load_con(index,3)==1) || (load_con(index,4)==1 && load_con(index,5)==1); 
               pwrmod_idx(j) = index;
           else
               error('pwrmod buses must be defined as 100% constant power or 100% constant current in load_con'); 
           end
           kk = bus_int(pwrmod_con(:,1));
           if any(abs(bus(kk,10)-2)); error('power modulation buses must be declared type 2 (PV) buses'); end
           if max(abs(bus(kk,6)))>1e-10 || max(abs(bus(kk,7)))>1e-10; error('power modulation buses must have a zero load'); end
           clear kk
       else
          error('you must have the power modulation bus declared as a 100% consant-power or 100% constant current non-conforming load')
       end
    end
end
       
    