function handleNetworkSln(k, flag)
%HANDLENETWORKSLN  save and restore network solution data
% HANDLENETWORKSLN saves or restores the network solution at data index k
%
% Syntax: handleNetworkSln(k, flag)
%
%   NOTES: Used to reset the newtork values to the initial solution in VTS.
%
%   Input:
%   k - data index to log from and restore to
%   flag - choose funtion operation
%       0 - initialize globals used to store data
%       1 - collect newtork solution values from index k into a global vector
%       2 - write stored network solution vector to network globals at data index k
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   08/03/20    09:03   Thad Haines     Version 1
%   08/12/20    20:55   Thad Haines     Version 1.0.1 - begun adding values changed in model call with flag == 1

%%
global g

if flag == 0
    %% initialize cell of global names for solution values
    netSlnCell = { ...
        % field, state, number
        'bus', 'bus_v', 0;
        'bus', 'theta', 0;
        'mac', 'cur_re', 0;
        'mac', 'cur_im', 0;
        % added 08/12/20-thad
        'mac', 'mac_ang', 0;
        'mac', 'psi_re', 0;
        'mac', 'psi_im', 0;
        %'mac', 'psidpp', 0; %single non time series data...
        %'mac', 'psiqpp', 0;
        'mac', 'vex', 0;
        'exc', 'V_A', 0;
        'exc', 'V_B', 0;
        };
    % induction motor
    if g.ind.n_mot ~= 0
        netSlnCell = [netSlnCell; {'ind', 'idmot',0}];
        netSlnCell = [netSlnCell; {'ind', 'iqmot',0}];
        netSlnCell = [netSlnCell; {'ind', 's_mot',0}];
        netSlnCell = [netSlnCell; {'ind', 'p_mot',0}];
        netSlnCell = [netSlnCell; {'ind', 'q_mot',0}];
        % added 08/12/20-thad
        netSlnCell = [netSlnCell; {'ind', 'vdmot',0}];
        netSlnCell = [netSlnCell; {'ind', 'vqmot',0}];
    end
    
    % induction generator
    if g.igen.n_ig ~=0
        netSlnCell = [netSlnCell; {'igen', 'idig',0}];
        netSlnCell = [netSlnCell; {'igen', 'iqig',0}];
        netSlnCell = [netSlnCell; {'igen', 's_igen',0}];
        netSlnCell = [netSlnCell; {'igen', 'pig',0}];
        netSlnCell = [netSlnCell; {'igen', 'qig',0}];
    end
    
    if g.dc.n_conv ~=0
        netSlnCell = [netSlnCell; {'dc', 'Vdc',0}];
        netSlnCell = [netSlnCell; {'dc', 'i_dc',0}];
        % added 08/12/20-thad
        netSlnCell = [netSlnCell; {'dc', 'alpha',0}];
        netSlnCell = [netSlnCell; {'dc', 'gamma',0}];
    end
    
    % added 08/12/20-thad
    if g.dc.ndcr_ud~=0
        netSlnCell = [netSlnCell; {'dc', 'angdcr',0}];
        netSlnCell = [netSlnCell; {'dc', 'dcrd_sig',0}];
        netSlnCell = [netSlnCell; {'dc', 'dcr_dsig',0}];
    end
    % added 08/12/20-thad
    if g.dc.ndci_ud~=0
        netSlnCell = [netSlnCell; {'dc', 'angdci',0}];
        netSlnCell = [netSlnCell; {'dc', 'dcid_sig',0}];
        netSlnCell = [netSlnCell; {'dc', 'dci_dsig',0}];
    end
    
    % count number of values to track
    totSlns = 0;
    for n=1:size(netSlnCell,1)
        f1 = netSlnCell{n,1};
        f2 = netSlnCell{n,2};
        netSlnCell{n,3} = size( g.(f1).(f2), 1);
        totSlns = totSlns+netSlnCell{n,3};
    end
    
    % save network solution cell to global
    g.vts.netSlnCell = netSlnCell;
    % initialize network solution place holder
    g.vts.netSlnVec = zeros( totSlns, 1);
    
    
elseif flag == 1
    %% collect solution values
    startN = 1; % used track index of network solution vector
    
    for ndx=1:size(g.vts.netSlnCell,1)
        
        % If the number of allocated solutions is larger than zero
        if g.vts.netSlnCell{ndx, 3} ~= 0
            % collect solution loction
            f1 = g.vts.netSlnCell{ndx, 1};
            f2 = g.vts.netSlnCell{ndx, 2};
            %fprintf('collecting... %s.%s \n',f1, f2) % DEBUG output
            % place solution into netSlnVec using dynamic field names
            g.vts.netSlnVec(startN:startN+g.vts.netSlnCell{ndx,3}-1) = g.(f1).(f2)(1:g.vts.netSlnCell{ndx,3}, k);
            
            startN = startN + g.vts.netSlnCell{ndx,3}; % increment starting index
        else
            % probably not executed...
            disp(['skpping... ', g.vts.netSlnCell{ndx,1}, g.vts.netSlnCell{ndx,2} ] );
        end
    end
    
    
elseif flag == 2
    %% write original collected values back to globals
    startN = 1; % used track index of network solution vector
    
    for ndx=1:size(g.vts.netSlnCell,1)
        % If the number of allocated solutions is larger than zero
        if g.vts.netSlnCell{ndx, 3} ~= 0
            % collect solution loction
            f1 = g.vts.netSlnCell{ndx, 1};
            f2 = g.vts.netSlnCell{ndx, 2};
            
            % place solution back into globals using dynamic field names
            for i=1:g.vts.netSlnCell{ndx, 3}
                g.(f1).(f2)(i,k) = g.vts.netSlnVec(startN+i-1);
            end
            
            startN = startN + g.vts.netSlnCell{ndx,3}; % increment starting index
        else
            % probably not executed...
            disp(['skpping... ', g.vts.netSlnCell{ndx,1}, g.vts.netSlnCell{ndx,2} ] );
        end
    end
else
    %invalid flag
    fprintf('*** Invalid Flag\n')
end
end% end function