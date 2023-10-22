function dydt = Bistability_ode(t, y, option, param)
%% -------------------------- PARAMETER MAPPING ----------------------------------%%
k0 = param.k0;
k1 = param.k1;
k21 = param.k21;
k22 = param.k22;
k3 = param.k3;
k4f = param.k4f;
k4b = param.k4b;
k4c = param.k4c;
k5 = param.k5;
PRXtot = param.PRXtot;
SRXtot = param.SRXtot;
TRX = param.TRX;
H2O2_null_curve_switch = param.H2O2_null_curve_switch;
PRXSO2H_null_curve_switch = param.PRXSO2H_null_curve_switch;
k3_switch = param.k3_switch;

%% ------------------------- STATE VARIABLE MAPPING----------------------------%%
H2O2 = y(1);
PRX	= y(2);
PRXSOH = y(3);
PRXSS = y(4);
PRXSO2H	= y(5);
PS = PRXtot - PRX- PRXSOH - PRXSS - PRXSO2H;
SRX = SRXtot - PS;

%% ------------------------------ ODEs-------------------------------------%%
dydt = zeros(length(y),1); % make dydt as a column vector as required by MatLab ode function

%H2O2
dydt(1) = PRXSO2H_null_curve_switch * (k0 - k1*H2O2*PRX - k3*(1-k3_switch+k3_switch*H2O2)*PRXSOH - k5*H2O2);   

%PRX
dydt(2) = - k1*H2O2*PRX + k22*TRX*PRXSS;

%PRXSOH
dydt(3) = k1*H2O2*PRX - k21*PRXSOH - k3*(1-k3_switch+k3_switch*H2O2)*PRXSOH + k4c*PS;

%PRXSS
dydt(4) = k21*PRXSOH - k22*TRX*PRXSS;

%PRXSO2H
dydt(5) = H2O2_null_curve_switch * (k3*(1-k3_switch+k3_switch*H2O2)*PRXSOH - k4f*SRX*PRXSO2H + k4b*PS);

end                      