%% Script to detect 'legacy' data input and adjust to new global g approach
% Usage of script to more easily detect legacy global definitons
% global g is already defined. Assumes data file strucutre does not change.
%
%   History:
%   Date        Time    Engineer        Description
%   06/02/20    08:51   Thad Haines     init - lmod_con
%   06/05/20    09:53   Thad Haines     addition of tg_con
%   06/12/20    11:53   Thad Haines     addition of livePlotFlag
%   06/15/20    14:21   Thad Haines     addition of rlmod_con
%   06/17/20    10:21   Thad Haines     addition of exc_con
%   06/18/20    14:21   Thad Haines     addition of mac_con and ibus_con
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
%% inifite bus
if exist('ibus_con','var')
    g.mac.ibus_con = ibus_con;
    clear ibus_con
else
    g.mac.ibus_con = [];
end
%% Global for plot flag
if exist('livePlotFlag','var')
    g.sys.livePlotFlag = livePlotFlag;
    clear livePlotFlag
else
    % default behavior
    g.sys.livePlotFlag = 1;
end