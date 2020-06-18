function exc_indx()
% Syntax: exc_indx()
% 9.02am 3/7/98
% forms indexes from exc_con to allow different
% exciter models to be used while retaining the vector option

% Version: 2.0
% Date:    July 1998
% Author:  Graham Rogers
% Modification: eliminate length checks in favour of isempty checks

% Version: 1.0
% Date:    June 1996
% Author:  Graham Rogers
% (c) copyright Joe Chow 1996, All rights reserved

%
% global exc_pot exc_con n_exc
% % simple exciter
% global smp_idx n_smp smppi_idx n_smppi dc_idx n_dc  dc2_idx n_dc2 st3_idx n_st3
% global smp_TA smp_TA_idx smp_noTA_idx smp_TB smp_TB_idx smp_noTB_idx
% global smp_TR smp_TR_idx smp_noTR_idx smppi_TR smppi_TR_idx smppi_noTR_idx
% % dc exciter
% global dc_TA dc_TA_idx dc_noTA_idx dc_TB dc_TB_idx dc_noTB_idx;
% global dc_TE dc_TE_idx dc_noTE_idx;
% global dc_TF dc_TF_idx  dc_TR dc_TR_idx dc_noTR_idx;
% % st3 exciter
% global st3_TA st3_TA_idx st3_noTA_idx st3_TB st3_TB_idx st3_noTB_idx;
% global st3_TR st3_TR_idx st3_noTR_idx;

% new global
global g

g.exc.n_smp = 0;
g.exc.n_smppi =0;
g.exc.n_dc = 0;
g.exc.n_dc1 = 0.; % not originally global -thad 06/17/20
g.exc.n_dc2 = 0;
g.exc.n_st3 = 0;
g.exc.n_exc = 0;

if ~isempty(g.exc.exc_con)
    %check for simple exciters
    g.exc.smp_idx = find(g.exc.exc_con(:,1)== 0);
    if ~isempty(g.exc.smp_idx)
        g.exc.n_smp = length(g.exc.smp_idx);
    end
    g.exc.smppi_idx = find(g.exc.exc_con(:,1)== 4);
    if ~isempty(g.exc.smppi_idx)
        g.exc.n_smppi = length(g.exc.smppi_idx);
    end
    
    %check for dc exciters
    g.exc.dc_idx = find( (g.exc.exc_con(:,1) == 1) || (g.exc.exc_con(:,1)==2) );
    if ~isempty(g.exc.dc_idx)
        g.exc.n_dc = length(g.exc.dc_idx);
    end
    
    g.exc.dc1_idx = find(g.exc.exc_con(:,1)==1);
    if ~isempty(g.exc.dc1_idx)
        g.exc.n_dc1 = length(g.exc.dc1_idx);
    end
    
    g.exc.dc2_idx = find(g.exc.exc_con(:,1)== 2);
    if ~isempty(g.exc.dc2_idx)
        g.exc.n_dc2 = length(g.exc.dc2_idx);
    end
    
    %check for type 3 exciter
    g.exc.st3_idx = find(g.exc.exc_con(:,1) == 3);
    if ~isempty(g.exc.st3_idx)
        g.exc.n_st3 = length(g.exc.st3_idx);
    end
    %
    %form  vectors for  time constants
    % Note: these do not account for situations where smp_TA == 0.001...
    % -thad
    if g.exc.n_smp ~= 0
        % TA
        g.exc.smp_TA = g.exc.exc_con(g.exc.smp_idx,5);
        g.exc.smp_TA_idx = find(g.exc.smp_TA>0.001);
        g.exc.smp_noTA_idx = find(g.exc.smp_TA<0.001);
        % TB & TC
        g.exc.smp_TB = g.exc.exc_con(g.exc.smp_idx,6);
        g.exc.smp_TB_idx = find(g.exc.smp_TB>0.001);
        g.exc.smp_noTB_idx = find(g.exc.smp_TB<0.001);
        % TR
        g.exc.smp_TR = g.exc.exc_con(g.exc.smp_idx,3);
        g.exc.smp_TR_idx = find(g.exc.smp_TR>0.001);
        g.exc.smp_noTR_idx = find(g.exc.smp_TR<0.001);
    end
    if g.exc.n_smppi~=0
        % TR
        g.exc.smppi_TR = g.exc.exc_con(g.exc.smppi_idx,3);
        g.exc.smppi_TR_idx = find(g.exc.smppi_TR>0.001);
        g.exc.smppi_noTR_idx = find(g.exc.smppi_TR<0.001);
    end
    if g.exc.n_dc ~= 0
        % TA
        g.exc.dc_TA = g.exc.exc_con(g.exc.dc_idx,5);
        g.exc.dc_TA_idx = find(g.exc.dc_TA >0.001);
        g.exc.dc_noTA_idx = find(g.exc.dc_TA<0.001);
        % TB & TC
        g.exc.dc_TB = g.exc.exc_con(g.exc.dc_idx,6);
        g.exc.dc_TB_idx = find(g.exc.dc_TB >0.001);
        g.exc.dc_noTB_idx = find(g.exc.dc_TB<0.001);
        % TE
        g.exc.dc_TE = g.exc.exc_con(g.exc.dc_idx,11);
        g.exc.dc_TE_idx = find(g.exc.dc_TE>0.001);
        g.exc.dc_noTE_idx = find(g.exc.dc_TE<0.001);
        % TF
        g.exc.dc_TF = g.exc.exc_con(g.exc.dc_idx,17);
        g.exc.dc_TF_idx = find(g.exc.dc_TF>0.001);
        % TR
        g.exc.dc_TR = g.exc.exc_con(g.exc.dc_idx,3);
        g.exc.dc_TR_idx = find(g.exc.dc_TR>0.001);
        g.exc.dc_noTR_idx = find(g.exc.dc_TR<0.001);
    end
    
    if g.exc.n_st3 ~=0
        % TA
        g.exc.st3_TA = g.exc.exc_con(g.exc.st3_idx,5);
        g.exc.st3_TA_idx = find(g.exc.st3_TA>0.001);
        g.exc.st3_noTA_idx = find(g.exc.st3_TA<0.001);
        % TB & TC
        g.exc.st3_TB = g.exc.exc_con(g.exc.st3_idx,6);
        g.exc.st3_TB_idx = find(g.exc.st3_TB>0.001);
        g.exc.st3_noTB_idx = find(g.exc.st3_TB<0.001);
        % TR
        g.exc.st3_TR = g.exc.exc_con(g.exc.st3_idx,3);
        g.exc.st3_TR_idx = find(g.exc.st3_TR>0.001);
        g.exc.st3_noTR_idx = find(g.exc.st3_TR<0.001);
    end
    %set size of exc_pot
    g.exc.n_exc = g.exc.n_smp + g.exc.n_smppi + g.exc.n_dc + g.exc.n_st3;
    g.exc.exc_pot = zeros(g.exc.n_exc,5);
end
