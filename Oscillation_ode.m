function dydt = Oscillation_ode(t, y, option, param)
%% -------------------------- PARAMETER MAPPING ----------------------------------%%
k0 = param.k0;
k1 = param.k1;
k21 = param.k21;
k22 = param.k22;
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
TRX = param.TRX;
k14 = param.k14;
k15 = param.k15;
k16 = param.k16;
HSP90 = param.HSP90;
Vratio = param.Vratio;
PRXtot = param.PRXtot;
%

%% ------------------------- STATE VARIABLE MAPPING----------------------------%%
H2O2cyto = y(1);
H2O2mito = y(2);
SRXcyto = y(3);
SRXOH = y(4);
HSPSRXOH = y(5);
SRXmito = y(6);
PS = y(7);
PRX	= y(8);
PRXSOH = y(9);
PRXSO2H = y(10);
PRXSS = PRXtot-PRX-PRXSOH-PRXSO2H-PS;
%

%% ------------------------------ ODEs-------------------------------------%%
dydt = zeros(length(y),1); %make dydt as a column vector as required by MatLab ode function

%H2O2cyto
dydt(1) = k6*Vratio*H2O2mito - k7*H2O2cyto - k6*H2O2cyto;       

%H2O2mito
dydt(2) = k0 - k1*H2O2mito*PRX - k3*H2O2mito*PRXSOH - k5*H2O2mito - k6*H2O2mito + k6/Vratio*H2O2cyto;  

%SRXcyto
dydt(3) = k8 - k9*SRXcyto - k10*SRXcyto*H2O2cyto+ k11*SRXOH;

%SRXOH
dydt(4) = k10*SRXcyto*H2O2cyto - k11*SRXOH - k12f*SRXOH*HSP90 + k12b*HSPSRXOH - k13*SRXOH;   

%HSPSRXOH
dydt(5) = k12f*SRXOH*HSP90 - k12b*HSPSRXOH - k14*HSPSRXOH - k15*HSPSRXOH;
   
%SRXmito
dydt(6) = k14/Vratio*HSPSRXOH - k16*SRXmito - k4f*SRXmito*PRXSO2H + k4b*PS + k4c*PS;   

%PS
dydt(7) = k4f*SRXmito*PRXSO2H - k4b*PS - k4c*PS - k4d*PS;
 
%PRX
dydt(8) = -k1*H2O2mito*PRX + k22*TRX*PRXSS;  

%PRXSOH 
dydt(9) = k1*H2O2mito*PRX - k21*PRXSOH - k3*H2O2mito*PRXSOH + k4c*PS; 

%PRXSO2H 
dydt(10) = k3*H2O2mito*PRXSOH - k4f*SRXmito*PRXSO2H + k4b*PS + k4d*PS; 
%

end

                                 