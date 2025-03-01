function dydt = PTRS_oscillation_ode(t, y, option, param)
%% -------------------------- PARAMETERS MAPPING ------------------------- %%
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
k4d = param.k4d;
k5 = param.k5;
k6 = param.k6;
k7 = param.k7;
k8 = param.k8;
k9 = param.k9;
k10 = param.k10;
k11 = param.k11;
k12f = param.k12f;
k12b = param.k12b;
k13 = param.k13;
k14 = param.k14;
k15 = param.k15;
k16 = param.k16;
PRXtot = param.PRXtot;
TRXtot = param.TRXtot;
TRtot = param.TRtot;
HSP90 = param.HSP90;
V_ratio = param.V_ratio;
%

%% ------------------------- STATE NAME MAPPING--------------------------- %%
H2O2_cyto = y(1);
H2O2_mito = y(2);
SRX_cyto = y(3);
SRXOH_cyto = y(4);
SRXOH_HSP90_cyto = y(5);
SRX_mito = y(6);
PRXSO2H_SRX = y(7);
PRXSH = y(8);
PRXSOH = y(9);
PRXSO2H = y(10);
TRXSS = y(11);
PRXSS = PRXtot - PRXSH - PRXSOH - PRXSO2H - PRXSO2H_SRX;
TRXSH = TRXtot - TRXSS;
%

%% ------------------------------ ODEs------------------------------------ %%
dydt = zeros(length(y),1); %make dydt as a column vector as required by MatLab ode function

%H2O2_cyto
dydt(1) = k6*V_ratio*H2O2_mito - k7*H2O2_cyto - k6*H2O2_cyto - k10*SRX_cyto*H2O2_cyto;       

%H2O2_mito
dydt(2) = k0 - k1*H2O2_mito*PRXSH - k3*H2O2_mito*PRXSOH - k5*H2O2_mito - k6*H2O2_mito + k6/V_ratio*H2O2_cyto;  

%SRX_cyto
dydt(3) = k8 - k9*SRX_cyto - k10*SRX_cyto*H2O2_cyto+ k11*SRXOH_cyto;

%SRXOH_cyto
dydt(4) = k10*SRX_cyto*H2O2_cyto - k11*SRXOH_cyto - k12f*SRXOH_cyto*HSP90 + k12b* SRXOH_HSP90_cyto - k13*SRXOH_cyto;   

%SRXOH_HSP90_cyto
dydt(5) = k12f*SRXOH_cyto*HSP90 - k12b* SRXOH_HSP90_cyto - k14*SRXOH_HSP90_cyto - k15*SRXOH_HSP90_cyto;
   
%SRX_mito
dydt(6) = k14/V_ratio*SRXOH_HSP90_cyto - k16*SRX_mito - k4f*SRX_mito*PRXSO2H + k4b*PRXSO2H_SRX + k4c*PRXSO2H_SRX;   

%PRXSO2H_SRX
dydt(7) = k4f*SRX_mito*PRXSO2H - k4b*PRXSO2H_SRX - k4c*PRXSO2H_SRX - k4d*PRXSO2H_SRX;
 
%PRXSH
dydt(8) = -k1*H2O2_mito*PRXSH + k2b*TRXSH*PRXSS;  

%PRXSOH 
dydt(9) = k1*H2O2_mito*PRXSH - k2a*PRXSOH - k3*H2O2_mito*PRXSOH + k4c*PRXSO2H_SRX; 

%PRXSO2H 
dydt(10) = k3*H2O2_mito*PRXSOH - k4f*SRX_mito*PRXSO2H + k4b*PRXSO2H_SRX + k4d*PRXSO2H_SRX; 

%TRXSS
dydt(11) = k2b*TRXSH*PRXSS  -  k2c*0.5*(TRtot + TRXSS + Km2c - ((TRtot + TRXSS + Km2c)^2 - 4*TRtot*TRXSS)^0.5);

%
end

                                 