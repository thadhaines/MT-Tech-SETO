Partial lit review  

SAND2015-11084R - Dynamic Simulation over Long time periods with 100% solar generation  
https://prod-ng.sandia.gov/techlib-noauth/access-control.cgi/2015/1511084r.pdf  
Summary: Uses eigen value analysis to investigate stability of integration methods applied to dynamic simulation of power systems with photovoltaic sources.  
Relevance: Uses PST to investigate similar topics as SETO project.  
Provides suggestions on integration methods and introduces a procedure to test integration method applicability based on analysis of system eigen values.

PST Manual and Software availbe from:  
https://www.ecse.rpi.edu/~chowj/  
Relevance: Open source transient simulation software  

Sub-Hour Solar Data for Power System Modeling from Static Spatial Variability Analysis  
56204 - NREL  
https://www.nrel.gov/docs/fy13osti/56204.pdf  
Summary:  Focus on statistically modelling sub-hour (minute) solar data and comparing to measured data.  
Provides some characteristics of variable solar irradiance.  
Relevance: Possible source for relevant event data to simulate  

MAFRIT == Multi-Area Frequency Response Integration Tool  
https://github.com/NREL/MAFRIT  
Summary:  Based on MATPOWER load flow solver.  
Models machines, governors, wind turbines and AGC in the long-term.  
Meant to function in the ms to minute range of dynamic simulation.  
Relevance:  Similar aim as this project with a focus on economic scheduling. Possible replacement/update of MIDAS?

64637 - NREL Investigating Power System Primary and Secondary Reserve Interaction under High Wind Power Penetration  
https://www.nrel.gov/docs/fy17osti/64637.pdf    
Summary:  Uses Flexible Energy Scheduling Tool for Integrating Variable Generation (FESTIV) and MAFRIT to investigate primary and secondary frequency response in a multi-area system.   
Relevance: Use case example of the MARFRIT NREL software package to do long-term simulation.  
While focusing on wind generation, similar concepts may apply to PV generation.  

PSLTDSim - Power System Long-Term Dynamic Simulator  
Github code source: https://github.com/thadhaines/PSLTDSim  
Master thesis location: https://github.com/thadhaines/Thesis-Release/blob/master/200501-haines-thesis.pdf  
Summary: Uses time sequence power flow, combined system frequency, governors, and AGC to model long-term power system dynamics in Python.  
Relies on PSLF for system dynamic and topographic information, as well load flow solver algorithm.  
Does not focus on transients / sub second system responses.   
Relevance:  Long-term simulation of power system dynamics.   
Shows that time-sequenced power flow can be used to model primary and secondary frequency response.  

POWER SYSTEM SIMULATION USING AN ADAPTIVE MODELING FRAMEWORK  
https://digitalcommons.mtech.edu/grad_rsch/76/  
Summary: Masters thesis describing software that switches between classical transient simulation and long-term time sequenced power flow simulation.  
Relevance: Possible approach/idea to consider for long-term simulation if variable time step/ multi-step integration proves unsatisfactory.  

Fast Frequency Response Concepts and Bulk Power System Reliability Needs  
NERC Inverter-Based Resource Performance Task Force  
https://www.nerc.com/comm/PC/InverterBased%20Resource%20Performance%20Task%20Force%20IRPT/Fast_Frequency_Response_Concepts_and_BPS_Reliability_Needs_White_Paper.pdf   
Summary: Provides background into basic frequency response and control, factors in rate of change of frequency (ROCOF), inertia effects, technology-specific FFR capabilities (wind turbine, solar, and battery...).   
Relevance: Provides information and illustrations of system impacts due to fast frequency response from various sources.  
Models of these technologies may be useful to consider for this project.  

The following slide decks may not provide much 'substantive' information, but do provide references with more detail of topics that may prove to be of interest.

Integrating High Levels of Variable Renewable Energy into Electric Power Systems  
68349 - NREL  
https://www.nrel.gov/docs/fy17osti/68349.pdf  
Summary: Overview of where things are, and are going (relative to publication date), in relation to integration of renewables/inverter based energy.   
Relevance: Provides challenges, solutions, and references related to variable renewable energy topics.  

Grid Integration of Variable Renewable Generation: Reliability Challenges and Solutions  
72615- NREL  
https://www.nrel.gov/docs/fy19osti/72615.pdf  
Summary:  Similar to other NREL slide deck  
Relevance: Info and sources on increasing solar/wind usage.    
Interesting comparison graph of reactive power capabilities of various sources (generators, inverters ... ).  
Places MAFRIT in dynamic simulation time scale spanning ms to multiple minutes  
