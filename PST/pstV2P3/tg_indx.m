%governor indexes
% fix of function name 5/14/20
function f=tg_indx
%set global variables
global  tg_con  
global  tg_idx n_tg tgh_idx n_tgh
f=0;
n_tg=0; 
n_tgh=0;
%if length(tg_con)~=0 5/14/20 edit
if ~isempty(tg_con)
   tg_idx = find(tg_con(:,1)==1);
   if ~isempty(tg_idx)
      n_tg = length(tg_idx);
   else
      n_tg = 0;
   end 
   tgh_idx = find(tg_con(:,1)==2);
   if ~isempty(tgh_idx)
      n_tgh = length(tgh_idx);
   else
      n_tgh=0;
   end
end