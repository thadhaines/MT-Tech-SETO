%%  dlab_01_subt.m
%   Thad Haines         EELE 5550
%   Program Purpose:    Create bus, line, mac_con and sw_con for sub transient PST Model. 

%   History:
%   03/11/18    09:31   init

%% Create setting matricies
bus = [ ...
%   1   2   3       4       5       6       7       8   9   A   B       C       D        E   F
%busNum V delta     P_gen   Q_gen P_load   Q_load   G   B  type Q_max  Q_minVrated_(kv) Vmax Vmin
    1   1   0       -0.25   0       0       0       0   0   1   100     -100  	20       2   0.5; %slack / infinite bus
    2   1   0       0       0       0       0       0   0   3   0       0       20       2   0.5;
    3   1   0       0       0       0       0       0   0   3   0       0       20       2   0.5;
    4   1   0       0.25    0       0       0       0   0   2   100     -100  	20       2   0.5];    %gen
    
line = [ ...
%   1   2   3   4       5   6   7       8   9   A   
%From   To  R   X       LC tp/ Phase   max min size  
    1   2   0   5       0   1   0       0   0   0;
    1   3   0   3       0   1   0       0   0   0;
    2   3   0   1       0   1   0       0   0   0;
    3   4   0   0.1     0   1   0       0   0   0];
    
mac_con = [ ...
%   1   2   3 [MVA] 4       5       6       7       8       9       A       B       C       D       E       F       16      17      18      19
% Mac_# b_# Sbase   x_l     r_a     x_d     x'_d    x"_d    T'_do   T"_do   x_q     x'_q    x"_q	T'_qo   T"_qo   H       d_o     d_1     b_#
    1   1   100     0       0       0       0.0005  0       0       0       0       0       0       0       0       2370    0       0       1;  %infinite bus
    2   4   100     0.15    0.0011  1.7     0.245   0.185   5.9     0.03    1.64    1.64    0.185   0.082   0.082   2.37    0       0       4];   % G2
    
sw_con = [ ...
%   1                   2   3   4   5   6   7   
    0                   0   0   0   0   0   1/(60*4);
    1                   3   2   0   0   0   1/(60*10);
    1+(2/60)            0   0   0   0   0   1/(60*10);
    1+(2+.1)/60         0   0   0   0   0   1/(60*10);
    3                   0   0   0   0   0   1/(60*2);
    10                  0   0   0   0   0   1/(60*2)];
  