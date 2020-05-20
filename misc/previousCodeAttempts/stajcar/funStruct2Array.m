% Matt Stajcar

% This function converts bus, branch, and gen parameters from type
% struct to type array. This enables the user to take advantage of PSTs
% load flow solver.

function [PSTbus,PSTbranch] = funStruct2Array(Bus,Branch,Gen,n,m)

    PSTbus = zeros(n,15);
                                        % PST bus matrix variables.
    for ii = 1:n
        PSTbus(ii,1) = Bus(ii).Num;     % bus.
        PSTbus(ii,2) = Bus(ii).Vmag;    % voltage.
        PSTbus(ii,3) = Bus(ii).AngD;    % angle.
        PSTbus(ii,4) = Bus(ii).Pgen;	% p_gen.
        PSTbus(ii,5) = Bus(ii).Qgen;	% q_gen.
        PSTbus(ii,6) = Bus(ii).Pload;	% p_load.
        PSTbus(ii,7) = Bus(ii).Qload;	% q_load.
        PSTbus(ii,8) = 0;               % G_shunt.
        PSTbus(ii,9) = 0;               % B_shunt.
        PSTbus(ii,10) = Bus(ii).Type;   % type.
        % PSTbus(ii,11) --> q_max.
        % PSTbus(ii,12) --> q_min.
%         if Bus(ii).Type == 1 || Bus(ii).Type == 2
%             for jj = 1:length(Gen)
%                 if Gen(jj).BusNum == Bus(ii).Num
%                     PSTbus(ii,11) = Gen(jj).Qmax;
%                     PSTbus(ii,12) = Gen(jj).Qmin;
%                 end
%             end
%         else
%             PSTbus(ii,11) = 0;
%             PSTbus(ii,12) = 0;
%         end
        PSTbus(ii,11) = 0;
        PSTbus(ii,12) = 0;
        PSTbus(ii,13) = Bus(ii).Vrated; % v_rated.
        PSTbus(ii,14) = Bus(ii).Vmax;   % v_max.
        PSTbus(ii,15) = Bus(ii).Vmin;   % v_min. 
    end
    clear ii jj
    
    PSTbranch = zeros(m,10);
                                                % PST line matrix
                                                % variables.
    for ii = 1:m
        PSTbranch(ii,1) = Branch(ii).FmBus;     % bus.
        PSTbranch(ii,2) = Branch(ii).ToBus;     % bus.
        PSTbranch(ii,3) = real(Branch(ii).Zpu); % r.
        PSTbranch(ii,4) = imag(Branch(ii).Zpu); % x.
        PSTbranch(ii,5) = real(Branch(ii).Bpu); % y.
        PSTbranch(ii,6) = 1.0;                  % tapratio.
        PSTbranch(ii,7) = 0;                    % tapphase.
        PSTbranch(ii,8) = 0;                    % tapmax.
        PSTbranch(ii,9) = 0;                    % tapmin.
        PSTbranch(ii,10) = 0;                   % tapsize.
    end  
    clear ii
    
end