PST Version 4 - current version: 4.0.0a3

(c) Montana Technological University / Thad Haines 2020

An updated/upgraded version of PST with fixes, new models (AGC, PWRMOD, IVM), 
functionality (Average frequency calculation, in-simulation line monitoring),
examples, noticeable speed up, and an experimental variable time step solution.

Uses a global structure 'g' to collect all values. Most to all files related to 
linear and non-linear simulation via s_simu  and sv_mgen_Batch have been udpated 
to use this new structure, however there may be some corner cases or older 
undocumented models left unchanged (they will likely no longer operate as desired 
without modification).

'Funtionalization' of the non-linear simulation allows for easier future 
alterations. For example, if there was a desire to implement a custom ODE solver 
(variable or fixed time step), there is clear place to add that functionality 
and code already in place to accommodate for user specified solution methods.

A number of bug fixes are imagined to be required while new code base is tested.

%   History:  
%   Date        Time    Engineer        Description  
%   08/13/20    12:09   Thad Haines     4.0.0a1 - initial ALPHA collection of version 4 code base  
% 	08/14/20 	11:49 	Thad Haines 	4.0.0a2 - alteration to initTblock, addition to resetORIG of mac_trip_logic, default fix of mexc_sig  
% 	08/15/20 	09:33	Thad Haines 	4.0.0a3 - fix of pwrmod to include multiple cell states, update of initTblock time vector creation, cleanup of s_simu and PSTv4 directory.  
