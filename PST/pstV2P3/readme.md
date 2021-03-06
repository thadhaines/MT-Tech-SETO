Contents (may be incomplete- from contents.m)

Temporary/Rewritten Files  
DataFile.m - Contains system information for transient stability runs   
mac_trip_logic - trip generators...   


Load Flow Functions  
calc - calculates load flow mismatch, checks convergence  
chq_lim - checks for generator VAR limits  
dc_lf - performs HVDC load flow calculations  
form_jac - forms the load flow Jacobian  
inv_lf - load flow computations for inverter  
lfdc - ac and HVDC load flow driver  
lfdcs - ac and HVDC load flow function  
lfdemo - ac load flow driver  
lftap - modifies transformer tap settings  
loadflow - performs ac load flow calculations  
rec_lf - load flow computations for the rectifiers  
vsdemo - voltage stability driver  
y_sparse - forms the sparse admittance matrix of the network  

Simulation Models  
dc_cont - dc converter, models hvdc pole controls  
dc_cur - calculates dc line currents - used in nc_load  
dc_indx - checks dc data for consistency and sets numbers and options  
dc_line - Models HVDC line dynamics  
dc_load - used in nc_load, associated with HVDC links  
exc_dc12 - dc exciter and AVR  
exc_indx - checks exciter data and presets exciter numbers and options  
exc_st3 - IEEE ST3 static exciter and AVR  
line_cur - calculates currents in lines from transient voltage records  
line_pq - calculates real and reactive powers in lines  
lm_indx - index for real load modulation  
lmod - modulates specified real loads  
pwrmod_p - real-power injection at spectified buses  
pwrmod_q - reac-power injection at spectified buses  
mac_em - 'Classical' generator   
mac_ib - infinite bus generator model  
mac_igen - induction generator model  
mac_ind - induction motor  
mac_indx - checks machine\data and sets numbers and options  
mac_sub - subtransient generator  
mac_tra - transient generator  
mdc_sig - modulation function for HVDC*  
mexc_sig - modulation function for exciter Vref*  
ml_sig - active load modulation*  
mpm_sig - generator mechanical power modulation*  
msvc_sig - modulation function for SVC reference input*  
mtg_sig - modulation function for turbine governor reference*  
nc_load - non-conforming loads, HVDC, SVC, load modulation network interface  
ns_file - determines total number of states in small signal stability  
p_cont - controls perturbations for linear model development  
p_exc - perturbs exciter models  
p_file - forms columns of state matrix, b, c, and d matrices  
p_m_file - forms permutation matrix for state matrix calculation  
p_pss - perturbs pss models  
p_tg - perturbs turbine/governor models  
pss - power system stabilizer  
pss_des - power system stabilizer design  
pss_indx - checks pss data and sets numbers and options  
pss_phse - calculates phase shift through pss  
pst_var - contains all global variables required for simulation  
rbus_ang - computes bus angle changes  
red_ybus - calculates reduced y matrices  
rlm_indx - index for reactive load modulation  
rlmod - modulates selected reactive loads  
rltf - calculates root locus for transfer function feedback around state space system  
rml_sig - forms modulation signal for reactive loads  
rootloc - calculates rootlocus for scalar feedback around state space system  
s_simu - transient simulation driver (replaced by s_simu_Batch)  
s_simu_Batch - modified transient simulation driver for batch runs  
sd_torque - calculates generator synchronizing and damping torques  
smpexc - simple exciter model  
stab_d - interactive pss design  
stab_f - pss frequency response (renamed stabf?)  
statef - frequency response from state space  
step_res - step response from state space  
svc - static VAR compensator  
svc_indx - index for svc  
svm_gen - small signal stability driver  
tg - turbine/governor  
tg_hydro - hydraulic turbine/governor  
tg_indx - turbine/governor index  
y_switch - organizes reduced y matrices  

'New', Undocumented, or Mystery functions...  
dc_sim -    
dc_sud - rectifier user defined damping control  
dc_vidc - updates Vdc and i_dc assuming ac bus voltage remains constant
dci_sud - inverter user defined damping control  
dcr_sud - rectifier user defined damping control   
dpwf - filter model for deltaP/w stabilizer  
dpwf_indx - Forms indexes for the deltaP/w filter  
hvdc_sud - HVDC user defined damping control   
i_simu - forms the network interface variables  
ind_ldto - Template for Induction Motor Load Torque Calculation as a function of slip   
insimit - simultaneous iteration on inverse A  
ivmmod_dyn - Implement state or output variables to model power injection   
mac_ivm - Internal Voltage Model type generator   
mtcsc_sig - modulation signal for tcsc control*   
nm_if - network-machine interface  
p_dpw - perturb the deltaP/omega filter variables  
pss_des_gain16 - used to select the PSS gain for the 16 machine system  
pwrm_indx - determines the relationship between pwrmod and nc loads  
pwrmod_dyn - Implement state or output variables to model power injection  
pwrmod_indx - determines the relationship between pwrmod and nc loads (duplicate?)  
smppi - simple excitation system with pi avr  
stabf - calculates the phase lag through the exciter and generator (same as stab_f?)  
svc_Open_Loop - static var system, built from svc?  
svc_sud - svc user defined damping control  
svm_mgen_Batch - modified svn_mgen   
swcap - no comments - useful?  
switch - switching point generation   
tcsc - thyristor controlled series capacitor  
tcsc_indx - determines the relationship between tcsc and nc loads  
tcsc_sud - tcsc user defined damping control  
time_stamp - returns parsed time string  
ybus - build admittance matrix Y from the line data (loadflow function)   