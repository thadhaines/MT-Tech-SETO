%perturb the turbine governor states 
% part of svm_mgen
% 5:08 PM 15/08/97 

k_tg = find(mac_tg==k);
if ~isempty(k_tg)
   disp('disturb turbine governor')
   gh = g.tg.tg_idx(k_tg);
   j=j+1;
   pert = 0.0001*abs(g.tg.tg1(gh,1));   
   pert = max(pert,0.0001);
   g.tg.tg1(gh,2) = g.tg.tg1(gh,1) + pert;
   p_file   % m file of perturbations
   st_name(k,j) = 21;
   
   j = j + 1;
   pert = 0.0001*abs(g.tg.tg2(gh,1));   
   pert = max(pert,0.0001);
   g.tg.tg2(gh,2) = g.tg.tg2(gh,1) + pert;
   p_file   % m file of perturbations
   st_name(k,j)=22;
   
   
   j=j+1;
   pert = 0.0001*abs(g.tg.tg3(gh,1));   
   pert = max(pert,0.0001);
   g.tg.tg3(gh,2) = g.tg.tg3(gh,1) + pert;
   p_file   % m file of perturbations
   st_name(k,j) = 23;
end
k_tgh = find(mac_tgh==k);
if ~isempty(k_tgh)
   disp('disturb hydraulic turbine governor')
   j=j+1;
   gh = g.tg.tgh_idx(k_tgh)
   pert = 0.0001*abs(g.tg.tg1(gh,1));   
   pert = max(pert,0.0001);
   g.tg.tg1(gh,2) = g.tg.tg1(gh,1) + pert;
   p_file   % m file of perturbations
   st_name(k,j) = 21;
   
   j = j + 1;
   pert = 0.0001*abs(g.tg.tg2(gh,1));   
   pert = max(pert,0.0001);
   g.tg.tg2(gh,2) = g.tg.tg2(gh,1) + pert;
   p_file   % m file of perturbations
   st_name(k,j)=22;
   
   
   j=j+1;
   pert = 0.0001*abs(g.tg.tg3(gh,1));   
   pert = max(pert,0.0001);
   g.tg.tg3(gh,2) = g.tg.tg3(gh,1) + pert;
   p_file   % m file of perturbations
   st_name(k,j) = 23;
   
   j=j+1;
   pert = 0.0001*abs(g.tg.tg4(gh,1));   
   pert = max(pert,0.0001);
   g.tg.tg4(gh,2) = g.tg.tg4(gh,1) + pert;
   p_file   % m file of perturbations
   st_name(k,j) = 24;
   
   j=j+1;
   pert = 0.0001*abs(g.tg.tg5(gh,1));   
   pert = max(pert,0.0001);
   
   g.tg.tg5(gh,2) = g.tg.tg5(gh,1) + pert;
   p_file   % m file of perturbations
   st_name(k,j) = 25;
   
end