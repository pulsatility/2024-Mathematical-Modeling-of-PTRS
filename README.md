#  Origins of Ultrasensitivity and Complex Signaling Dynamics of Cellular Hydrogen Peroxide and Peroxiredoxin
Shengnan Liu<sup>1,2,3</sup>, Jingbo Pi<sup>1,2,3</sup> and Qiang Zhang<sup>4</sup>

Published in Antioxidants (2025) https://doi.org/10.3390/antiox14020235

1	Key Laboratory of Environmental Stress and Chronic Disease Control & Prevention Ministry of Education, China Medical University, Shenyang, 110122, China

2	Key Laboratory of Liaoning Province on Toxic and Biological Effects of Arsenic, China Medical University, Shenyang, 110122, China

3	Program of Environmental Toxicology, School of Public Health, China Medical University, Shenyang, 110122, China

4	Gangarosa Department of Environmental Health, Rollins School of Public Health, Emory University, Atlanta, GA, 30322, USA

**Abstract:**
Hydrogen peroxide (H2O2) plays a crucial role in cell signaling in response to physiological and environmental perturbations. H2O2 can oxidize typical 2-Cys peroxire-doxin (PRX) first into a sulfenic acid, which resolves into a disulfide that can be reduced by thioredoxin (TRX)/TRX reductase (TR). At high levels, H2O2 can also hyperoxidize sul-fenylated PRX into a sulfinic acid that can be reduced by sulfiredoxin (SRX). Therefore, PRX, TRX, TR, and SRX (abbreviated as PTRS system here) constitute the coupled sul-fenylation and sulfinylation cycle (CSSC), where certain oxidized PRX and TRX forms also function as redox signaling intermediates. Earlier studies have revealed that the PTRS system is capable of rich signaling dynamics, including linearity, ultrasensitivi-ty/switch-like response, nonmonotonicity, circadian oscillation, and possibly, bistability. However, the origins of ultrasensitivity, which is fundamentally required for redox signal amplification, have not been adequately characterized, and their roles in enabling complex nonlinear dynamics of the PTRS system remain to be determined. Through in-depth mathematical modeling analyses, here we revealed multiple sources of ultrasensitivity that are intrinsic to the CSSC, including zero-order kinetic cycles, multistep H2O2 signaling, and a mechanism arising from diminished H2O2 removal at high PRX hyperoxidation state. The CSSC, structurally a positive feedback loop, is capable of bistability under cer-tain parameter conditions, which requires embedding multiple sources of ultrasensitivity identified. Forming a negative feedback loop with cytosolic SRX as previously observed in energetically active cells, the mitochondrial PTRS system (where PRX3 is expressed) can produce sustained circadian oscillations through supercritical Hopf bifurcations. In con-clusion, our study provided novel quantitative insights into the dynamical complexity of the PTRS system and improved appreciation of intracellular redox signaling.

**Keywords:** H2O2; peroxiredoxin; ultrasensitivity; bistability; feedback; circadian rhythm

#  MATLAB Code
- PTRS_ultrasensitivity_cmd.m: Main MATLAB code to run the Ultrasensitivity Model.
- PTRS_bistability_cmd.m: Main MATLAB code to run the Bistability Model.
- PTRS_ode.m: MATLAB ODE code to be called by PTRS_ultrasensitivity_cmd.m, PTRS_bistability_cmd.m.
- PTRS_ode_Fig3E.m: MATLAB ODE code to be called by PTRS_ultrasensitivity_cmd.m to generate Fig. 3E.
- Figure_8.m: MATLAB code to generate Fig. 8.
- Figure_9_and_S6.m: MATLAB code to generate Fig. 9 and Fig. S6.
- P3TRS_ode: MATLAB ODE code to be called by Figure_9_and_S6.m.
- PTRS_oscillation_cmd.m: Main MATLAB code to run the Oscillation Model.
- PTRS_oscillation_ode.m: MATLAB ODE code to be called by PTRS_oscillation_cmd.m.
- All *.xlsx files: XPP-AUT-generated bifurcation results to be called by PTRS_ultrasensitivity_cmd.m, PTRS_bistability_cmd.m, Figure_8.m, Figure_9_and_S6.m, and PTRS_oscillation.m to plot the bifurcation diagrams.
