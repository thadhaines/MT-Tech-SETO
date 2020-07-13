%% Script to detect 'legacy' data input and adjust to new global g approach
% Usage of script to more easily detect legacy global definitons
% global g is already defined. Assumes data file structure does not change.
%
%   History:
%   Date        Time    Engineer        Description
%   06/02/20    08:51   Thad Haines     init - lmod_con
%   06/05/20    09:53   Thad Haines     addition of tg_con
%   06/12/20    11:53   Thad Haines     addition of livePlotFlag
%   06/15/20    14:21   Thad Haines     addition of rlmod_con
%   06/17/20    10:21   Thad Haines     addition of exc_con
%   06/18/20    14:21   Thad Haines     addition of mac_con and ibus_con
%   06/30/20    10:18   Thad Haines     addition of pwrmod_con
%   07/02/20    10:10   Thad Haines     addition of lmon_con
%   07/02/20    13:15   Thad Haines     addition of load_con for ncl
%   07/06/20    09:30   Thad Haines     addition of sw_con
%   07/06/20    11:15   Thad Haines     addition of pss_con and gain fix
%   07/08/20    15:32   Thad Haines     addition of svc_con
%   07/09/20    11:11   Thad Haines     addition of tcsc_con
%   07/13/20    09:44   Thad Haines     addition of igen_con

global g

%% lmod
if exist('lmod_con','var')
    g.lmod.lmod_con = lmod_con;
    clear lmod_con 
else
    g.lmod.lmod_con = [];
    g.lmod.n_lmod = 0;
end

%% rlmod
% indicies created in rlm_indx.m
if exist('rlmod_con','var')
    g.rlmod.rlmod_con = rlmod_con;
    clear rlmod_con 
else
    g.rlmod.rlmod_con = [];
end

%% tg
if exist('tg_con','var')
    g.tg.tg_con = tg_con;
    clear tg_con
else
    g.tg.tg_con = [];
end

%% exciter
if exist('exc_con','var')
    g.exc.exc_con = exc_con;
    clear exc_con
else
    g.exc.exc_con = [];
end
%% machines
if exist('mac_con','var')
    g.mac.mac_con = mac_con;
    clear mac_con
else
    g.mac.mac_con = [];
end
%% inifite bus (not really used?)
if exist('ibus_con','var')
    g.mac.ibus_con = ibus_con;
    clear ibus_con
else
    g.mac.ibus_con = [];
end
%% pwrmod
if exist('pwrmod_con','var')
    g.pwr.pwrmod_con = pwrmod_con;
    clear pwrmod_con;
else
    g.pwr.pwrmod_con = [];
end
%% lmon_con
if exist('lmon_con','var')
    g.sys.lmon_con = lmon_con;
    clear lmon_con;
else
    g.sys.lmon_con = [];
end
%% load_con
if exist('load_con','var')
    g.ncl.load_con = load_con;
    clear load_con;
else
    g.ncl.load_con = [];
end
%% sw_con
if exist('sw_con','var')
    g.sys.sw_con = sw_con;
    clear sw_con;
else
    g.sys.sw_con = [];
end

%% pss
if exist('pss_con','var')
    % account for washout gain alteration between version 2 and 3
    if exist('pssGainFix','var')
        if pssGainFix
            pss_con(:,3) = pss_con(:,3)./pss_con(:,4);
        end
    else
        pssGainFix = 0;
    end
    g.pss.pss_con = pss_con;
    g.pss.pssGainFix = pssGainFix;
    clear pss_con pssGainFix
else
    g.pss.pss_con = [];
    g.pss.pssGainFix = nan;
end

%% svc_con
if exist('svc_con','var')
    g.svc.svc_con = svc_con;
    clear svc_con;
else
    g.svc.svc_con = [];
end
%% tcsc_con
if exist('tcsc_con','var')
    g.tcsc.tcsc_con = tcsc_con;
    clear tcsc_con;
else
    g.tcsc.tcsc_con = [];
    g.tcsc.n_tcsc = 0;
end
%% igen_con
if exist('igen_con','var')
    g.igen.igen_con = igen_con;
    clear igen_con;
else
    g.igen.igen_con = [];
end
%% Global for plot flag
if exist('livePlotFlag','var')
    g.sys.livePlotFlag = livePlotFlag;
    clear livePlotFlag
else
    % default behavior
    g.sys.livePlotFlag = 1;
end