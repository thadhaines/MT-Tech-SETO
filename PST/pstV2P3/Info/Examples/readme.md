PST Examples  

Running comments about PST:  

s_simu_Batch modified to show warnings at various simulation points - the warning also supplies line numbers which are useful  

When k == 50, i.e. the 50th simulation loop interation, any call to the  network solution i_simu is identified with a warning.  
This shows that the network is solved twice...   

[sInjF,sInjT]=line_pq(.... added post simulation to produce line flow values  


Current issues:

exclusive use of global variables  
batch run assumes 60 Hz, 100 MVA base  
switching of fault status via replacing Y matricies requires further study (seems kinda janky)  
no identifier for multiple lines connecting to same bus  
no area definitions (will be required for AGC)  
no string identifier for any object  
no running log of bus load powers, current injections*, branch flows* (*done post simulation)  
seperation of data from object identifier -> labeling/identification requires cross referencing  
matrix definitions get large/clunky/unweildy  
Transformers just reactance?  
Tap changers exist...   

easy fix:  
angles require a reference to the slack bus angle and be adjusted appropriately to produce 'standard reference' plot

Things possibly to do with integration routine:  
    i_simu called multiple times per time step (shown if DEBUG = 1 in d_ file)  
    bus_v written multiple times per time step (part of i_simu)  
    time step larger than 1 cycle causes system to 'blow up' i.e. too large of an integration step  