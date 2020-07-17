function pss_indx()
%PSS_INDX forms indexes for the pss
% PSS_INDX Forms indexes for the pss and determines indexs for
% the generators and exciters to which pss are connected.
%
% Syntax: pss_indx
%
%   Input: 
%   VOID
%
%   Output: 
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   08/xx/97    xx:xx   Graham Rogers   Version 1.2
%   (c) Copyright Joe Chow/Cherry Tree Scientific Software 1997 All Rights Reserved
%   07/06/20    14:02   Thad Haines     Revised format of globals and internal function documentation

global g

if ~isempty(g.pss.pss_con)
  g.pss.pss_idx = find(g.pss.pss_con(:,1)==1|g.pss.pss_con(:,1)==2);
  g.pss.n_pss = length(g.pss.pss_idx);
  g.pss.pss_mb_idx = g.mac.mac_int(round(g.pss.pss_con(:,2)));
  for jpss = 1:g.pss.n_pss
     g.pss.pss_exc_idx(jpss) = find(round(g.pss.pss_con(jpss,2))==round(g.exc.exc_con(:,2))); 
     if isempty(g.pss.pss_exc_idx(jpss))
         error('you must have an exciter at the same generator as a pss');
     end  
  end
  if g.pss.n_pss~=0
     g.pss.pss_T = g.pss.pss_con(g.pss.pss_idx,4);
     g.pss.pss_T2 = g.pss.pss_con(g.pss.pss_idx,6);
     g.pss.pss_T4 = g.pss.pss_con(g.pss.pss_idx,8);
     g.pss.pss_T4_idx = find(g.pss.pss_T4>0.001);
     g.pss.pss_noT4_idx = find(g.pss.pss_T4<0.001);
     g.pss.pss_sp_idx = find(g.pss.pss_con(g.pss.pss_idx,1)==1);
     g.pss.pss_p_idx = find(g.pss.pss_con(g.pss.pss_idx,1)==2);
  end
else
  g.pss.n_pss = 0;
end