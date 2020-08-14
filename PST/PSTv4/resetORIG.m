%% File resets all user defined or specified files to original state.
% format of command:  copyfile(source, destination);

% NOTE: must be run from root PST directory

%   History:
%   Date        Time    Engineer        Description
%   07/07/20    15:09   Thad Haines     Init
%   07/17/20    13:31   Thad Haines     Added live plotting
%   08/14/20    11:19   Thad Haines     Added machine trip logic

%% Model Dynamic Configuration
% IVM MOD
copyfile('ivmmod_dyn_ORIG.m', 'ivmmod_dyn.m');

% PWR MOD
copyfile('pwrmod_dyn_ORIG.m', 'pwrmod_dyn.m');

%% Alternate Models
% MAC_SUB
copyfile('mac_sub_ORIG.m','mac_sub.m');

% PSS
copyfile('pss2.m','pss.m');

%% generator trip logic
copyfile('mac_trip_logic_ORIG.m', 'mac_trip_logic.m');

%% Moduation Files
% MEXC_SIG - Exciter Signal
copyfile('mexc_sig_ORIG.m', 'mexc_sig.m');

% ML_SIG - Real load
copyfile('ml_sig_ORIG.m', 'ml_sig.m');

% MPM_SIG - Mechanical Power
copyfile('mpm_sig_ORIG.m', 'mpm_sig.m');

% MSVC_SIG - SVC Signal
copyfile('msvc_sig_ORIG.m', 'msvc_sig.m');

% MTCSC_SIG - TCSC Signal
copyfile('mtcsc_sig_ORIG.m', 'mtcsc_sig.m');

% MTG_SIG - Turbine Governor Signal (Pref)
copyfile('mtg_sig_ORIG.m', 'mtg_sig.m');

% RML_SIG - Reactive load
copyfile('rml_sig_ORIG.m', 'rml_sig.m');

%% Live plotting
copyfile('livePlot_ORIG.m', 'livePlot.m'); % use altered as default

%%
display('Models, Modulation files, and life plotting restored')