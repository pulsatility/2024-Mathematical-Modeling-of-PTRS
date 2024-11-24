#  MATLAB and XPP-AUT Code for manuscript "Origins of Ultrasensitivity and Complexity of Cellular Redox Signaling Dynamics"

Purpose of files
- Bistability_cmd.m: Main MATLAB code to generate bistability model results.
- Bistability_ode.m: MATLAB ODE code to be called by Bistability_cmd.m
- Relaxation_oscillation_cmd.m: Main MATLAB code to generate relaxation oscillator model results.
- Resonant_oscillation_cmd.m: Main MATLAB code to generate resonant oscillator model results.
- Oscillation_ode.m: MATLAB ODE code to be called by Relaxation_oscillation_cmd.m and Resonant_oscillation_cmd.m
- XPP-AUT/XPPAUT_Bistability.ode: XPP-AUT code to generate saddle-node bifurcation for bistability model
- XPP-AUT/XPPAUT_Relaxation_Oscillation.ode: XPP-AUT code to generate subcritical Hopf bifurcation for relaxation oscillator model
- XPP-AUT/XPPAUT_Resonant_Oscillation.ode: XPP-AUT code to generate supercritical Hopf bifurcation for resonant oscillator model
- XPP-AUT/XPPAUT_Relaxation_Oscillation_H2O2_null_curve.ode: XPP-AUT code to generate H2O2 null curve against Mito SRXtot for relaxation oscillator model
- XPP-AUT/XPPAUT_Resonant_Oscillation_H2O2_null_curve.ode: XPP-AUT code to generate H2O2 null curve against Mito SRXtot for resonant oscillator model
- XPP-AUT/*.xlsx: XPP-AUT outputs used to plot bifurcation or null curves as indicated in MATLAB *_cmd files.
