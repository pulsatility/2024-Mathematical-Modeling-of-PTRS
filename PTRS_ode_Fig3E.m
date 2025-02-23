function dydt = PTRS_ode_Fig3E(t, y, option, param)
%% -------------------------- PARAMETER MAPPING ----------------------------------%%
k0 = param.k0;
k1 = param.k1;
k2a = param.k2a;
k2b = param.k2b;
k2c = param.k2c;
Km2c = param.Km2c;
k3 = param.k3;
k4f = param.k4f;
k4b = param.k4b;
k4c = param.k4c;
k5 = param.k5;
PRXtot = param.PRXtot;
SRXtot = param.SRXtot;
TRXtot = param.TRXtot;
TRtot = param.TRtot;
PRXSO2H_switch = param.PRXSO2H_switch;
TRXSS_switch = param.TRXSS_switch;
PRXSS_switch = param.PRXSS_switch;
PS_switch = param.PS_switch;
PRXSH_switch = param.PRXSH_switch;
H2O2_switch = param.H2O2_switch;
k3_switch = param.k3_switch;
%

%% ------------------------- STATE VARIABLE MAPPING----------------------------%%
H2O2 = y(1);
PRXSH = y(2);
PRXSOH = y(3);
PRXSS = y(4);
PRXSO2H	= y(5);
TRXSS = y(6);
% PS = PS_switch * (PRXtot - PRXSH - PRXSOH - PRXSS - PRXSO2H);
PS = (1 - PS_switch) * PRXSO2H * SRXtot / ((k4c + k4b) / k4f + PRXSO2H); % used to generate Fig.3E;
PRXSH = PRXtot - PS - PRXSOH - PRXSS - PRXSO2H;
SRX = SRXtot - PS;
TRXSH = TRXtot - TRXSS;
%

%% ------------------------------ ODEs-------------------------------------%%
dydt = zeros(length(y),1); % make dydt as a column vector as required by MatLab ode function

%H2O2
dydt(1) = H2O2_switch * (k0 - k1*H2O2*PRXSH - k3*(1-k3_switch+k3_switch*H2O2)*PRXSOH - k5*H2O2);   

%PRXSH
dydt(2) = PRXSH_switch * (-k1*H2O2*PRXSH + k2b*TRXSH*PRXSS);

%PRXSOH
dydt(3) = k1*H2O2*PRXSH - k2a*PRXSOH - k3*(1-k3_switch+k3_switch*H2O2)*PRXSOH + k4c*PS;

%PRXSS
dydt(4) = PRXSS_switch * (k2a*PRXSOH - k2b*TRXSH*PRXSS);

%PRXSO2H
dydt(5) = PRXSO2H_switch * (k3*(1-k3_switch+k3_switch*H2O2)*PRXSOH - k4f*SRX*PRXSO2H + k4b*PS);

%TRXSS
dydt(6) = TRXSS_switch * (k2b*TRXSH*PRXSS - k2c*0.5*(TRtot + TRXSS + Km2c - ((TRtot + TRXSS + Km2c)^2 - 4*TRtot*TRXSS)^0.5));
%

end                      