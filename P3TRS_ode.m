function dydt = P3TRS_ode(t, y, option, param)
%% -------------------------- PARAMETER MAPPING ------------------------------- %%
k0 = param.k0;
k1 = param.k1;
k2a = param.k2a;
k2b = param.k2b;
k2c = param.k2c;
k3 = param.k3;
k4f = param.k4f;
k4b = param.k4b;
k4c = param.k4c;
k5 = param.k5;
PRXtot = param.PRXtot;
TRXtot = param.TRXtot;
TRtot = param.TRtot;
Km2c = param.Km2c;
V_ratio = param.V_ratio;
k6 = param.k6;
k7 = param.k7;
SRXtot = param.SRXtot;

%

%% ------------------------- STATE VARIABLE MAPPING---------------------------- %%
H2O2 = y(1);
H2O2_cyto = y(2);
PRXSH = y(3);
PRXSOH = y(4);
PRXSS = y(5);
PRXSO2H	= y(6);
TRXSS = y(7);

PS = PRXtot - PRXSH - PRXSOH - PRXSO2H - PRXSS;
TRXSH = TRXtot - TRXSS;

%

%% ------------------------------ ODEs----------------------------------------- %%
dydt = zeros(length(y),1); % make dydt as a column vector as required by MatLab ode function

%H2O2
dydt(1) = k0 - k1*H2O2*PRXSH - k3*H2O2*PRXSOH - k5*H2O2 - k6*H2O2 + k6/V_ratio*H2O2_cyto;  

%H2O2_cyto
dydt(2) = k6*V_ratio*H2O2 - k7*H2O2_cyto - k6*H2O2_cyto; 

%PRXSH
dydt(3) = - k1*H2O2*PRXSH + k2b*TRXSH*PRXSS;

%PRXSOH
dydt(4) = k1*H2O2*PRXSH - k2a*PRXSOH - k3*H2O2*PRXSOH + k4c*PS;

%PRXSS
dydt(5) = k2a*PRXSOH - k2b*TRXSH*PRXSS; 

%PRXSO2H
dydt(6) = k3*H2O2*PRXSOH - k4f*(SRXtot - PS)*PRXSO2H + k4b*PS;

%TRXSS
dydt(7) = k2b*TRXSH*PRXSS - k2c*0.5*(TRtot + TRXSS + Km2c - ((TRtot + TRXSS + Km2c)^2 - 4*TRtot*TRXSS)^0.5);

%

end                      