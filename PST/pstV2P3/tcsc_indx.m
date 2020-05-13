function f = tcsc_indx(tcsc_dc)
% syntax: f = tcsc_indx
% 10:29 am December 29, 1998
% determines the relationship between tcsc and nc loads
% checks for tcsc
% determines number of tcscs
% checks for user defined damping controls
% f is a dummy variable
f = 0;
global tcsc_con load_con  n_tcsc  tcscf_idx tcsct_idx % tcsc 
global n_tcscud dtcscud_idx  %user defined damping controls
n_tcsc = 0;
tcsc_idx = [];
dtcscud_idx = [];
n_tcscud = 0;
if ~isempty(tcsc_con)
   [n_tcsc npar] = size(tcsc_con);
   tcsc_idx = zeros(n_tcsc,1);
   for j = 1:n_tcsc
      index = find(tcsc_con(j,2)==load_con(:,1));
      if ~isempty(index)
         tcscf_idx(j) = index;
      else
         error('you must have the tcsc buses declared as a non-conforming load')
      end
      index = find(tcsc_con(j,3)==load_con(:,1));
      if ~isempty(index)
         tcsct_idx(j) = index;
      else
         error('you must have the tcsc buses declared as a non-conforming load')
      end
   end
   % check for user defined controls
   if ~isempty(tcsc_dc)
      [n_tcscud,dummy] = size(tcsc_dc);
      for j = 1:n_tcscud
         dtcscud_idx(j) = tcsc_dc{j,2};
      end
   end      
end

