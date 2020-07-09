% file to attempt to batch compare all file from known directories.
% NOTE: uses MATLABs visual comparison tool -  MANY tabs may open 
%       depending on number of files that are compared.

clear;close all; clc

% specify folders
dirNameA = 'pstV2P3';
dirNameB = 'pstV3';

dirA = dir(dirNameA);
dirB = dir(dirNameB);

% find longer directory
if length(dirA)> length(dirB)
    longDir = dirA;
    shortDir = dirB;
else
    longDir = dirB;
    shortDir = dirA;
end
clear dirA dirB

% Find common File names
commonName ={};
% for each directory entry in the long directory
for LdirN=3:length(longDir) % start at 3 to skip '.' and '..'
    % collect name to search for
    searchName = longDir(LdirN).name;
    % search short directory for match
    for searchN = 1:length(shortDir)
        if strcmp(searchName, shortDir(searchN).name );
            % match found
            commonName{end+1}= shortDir(searchN).name;
            continue
        end
    end
end
clear searchN LdirN serachName        

% Compare
for comNDX=1:length(commonName)
    fPath1 = ['./',dirNameA,'/',commonName{comNDX}];
    fPath2 = ['./',dirNameB,'/',commonName{comNDX}];
    visdiff(fPath1, fPath2)
    disp('hold for debug...')
end

%% Initial v2 v3 results of possile note
%{

dcr_sud - different handling of 's' and setting of number of states
exc_dc12 - saturation block alteration by JHC in 2015
i_simu - debug warning around line 80 notifying where bus_v is logged
inv_lf - output suppression differences only?
loadflow - possibly just format cleanup
! mac_ind - updated version in pst3 - added deep bar, double cage, saturation...
mac_indx - more logic checks - ivm in v2
mac_tra - previously noted x'd==x'q, add pm_sig in v3
mdc_sig - example of dc modulation?
nc_load - dif mainly due to pwrmod... I think
p_cont - scaling of tcsc, pwrmod, c_state indexing?!
red_ybus - pwrmod, spacing and code line dif otherwise
rlmod - handling of rlmod_st during anti-windup
tcsc_indx - more initialization of globals - incorporated into SETO version

Essentially: mac_ind will probably require 2 models - similar to pss, but
 maybe not - just a thought
%}