
% Copyright  Wang Y and Zhu X-G, 2014. 
% CAS-MPG Partner Institute for Computational Biology, Shanghai Institutes for Biological Sciences, CAS, Shanghai,200031 
% A mixed pathway model of C4 photosynthesis
% The rate of each reaction and transport process. 
 
function Enz_v=C4Vel(t,s)
 
global KValue;
% structure to store the Michaelis-Menten kinetic parameters
KmCO2_1=KValue(1,1);  Ke_1=KValue(1,2);
KmHCO3_2=KValue(2,1);  KmPEP_2=KValue(2,2);   Kimal_2=KValue(2,3);
KmNADPH_3=KValue(3,1);  KmOAA_3=KValue(3,2);  KmNADP_3=KValue(3,3);  Kmmal_3=KValue(3,4);  Ke_3=KValue(3,5);
KmCO2_4=KValue(4,1);  KmNADP_4=KValue(4,2);  KmNADPH_4=KValue(4,3);  KmPyr_4=KValue(4,4);  Kmmal_4=KValue(4,5);  Ke_4=KValue(4,6);
KiPEP_5=KValue(5,1);  KmATP_5=KValue(5,2);  KmPyr_5=KValue(5,3);
KmCO2_6=KValue(6,1);  KmO2_6=KValue(6,2);  KmRuBP_6=KValue(6,3);  KiPGA_6=KValue(6,4);  KiFBP_6=KValue(6,5);  KiSBP_6=KValue(6,6);  KiPi_6=KValue(6,7);  KiNADPH_6=KValue(6,8);
KmADP_7=KValue(7,1);  KmATP_7=KValue(7,2);  KmPGA_7=KValue(7,3);
KmDPGA_8=KValue(8,1);  KmNADPH_8=KValue(8,2);
Ke_9=KValue(9,1); 
KmDHAP_10=KValue(10,1);  KmFBP_10=KValue(10,2);  KmGAP_10=KValue(10,3);  Ke_10=KValue(10,4);
KiF6P_11=KValue(11,1);  KiPi_11=KValue(11,2);  KmFBP_11=KValue(11,3);  Ke_11=KValue(11,4);
KmDHAP_12=KValue(12,1);  KmE4P_12=KValue(12,2);  Ke_12=KValue(12,3);
KiPi_13=KValue(13,1);  KmSBP_13=KValue(13,2);  Ke_13=KValue(13,3);
KmE4P_14=KValue(14,1);  KmF6P_14=KValue(14,2);  KmGAP_14=KValue(14,3);  KmXu5P=KValue(14,4);  Ke_14=KValue(14,5);
KmGAP_15=KValue(15,1);  KmRi5P_15=KValue(15,2);  KmS7P_15=KValue(15,3);  KmXu5P_15=KValue(15,4);  Ke_15=KValue(15,5);
Ke_16=KValue(16,1);
Ke_17=KValue(17,1);
KiADP_18=KValue(18,1);  Ki_ADP_18=KValue(18,2);  KiPGA_18=KValue(18,3);  KiPi_18=KValue(18,4);  KiRuBP_18=KValue(18,5);  KmATP_18=KValue(18,6);  KmRu5P_18=KValue(18,7);  Ke_18=KValue(18,8);

KmADP_7Mchl=KValue(19,1);  KmATP_7Mchl=KValue(19,2);  KmPGA_7Mchl=KValue(19,3);
KmDPGA_8Mchl=KValue(20,1);  KmNADPH_8Mchl=KValue(20,1);

KiADP_Starch=KValue(21,1);  KmATP_Starch=KValue(21,2);  KmG1P_Starch=KValue(21,3);  KaF6P_Starch=KValue(21,4);  KaFBP_Starch=KValue(21,5);  KaPGA_Starch=KValue(21,6);  Ke_Starch1=KValue(21,7);   Ke_Starch2=KValue(21,8);
KmPGA_PGASink=KValue(22,1);
KmDHAP_Suc1=KValue(23,1);  KmGAP_Suc1=KValue(23,2);  KmFBP_Suc1=KValue(23,3);  Ke_Suc1=KValue(23,4);
KiF26BP_Suc2=KValue(24,1);  KiF6P_Suc2=KValue(24,2);  KiPi_Suc2=KValue(24,3);  KmFBP_Suc2=KValue(24,4);  Ke_Suc2=KValue(24,5);
Ke_Suc5=KValue(25,1);  Ke_Suc6=KValue(25,2);
KmG1P_Suc7=KValue(26,1);  KmPPi_Suc7=KValue(26,2);  KmUDPG_Suc7=KValue(26,3);  KmUTP_Suc7=KValue(26,4);  Ke_Suc7=KValue(26,5);
KiFBP_Suc8=KValue(27,1);  KiPi_Suc8=KValue(27,2);  KiSuc_Suc8=KValue(27,3);  KiSucP_Suc8=KValue(27,4);  KiUDP_Suc8=KValue(27,5);  KmF6P_Suc8=KValue(27,6);  KmUDPG_Suc8=KValue(27,7);  Ke_Suc8=KValue(27,8); 
KmSuc_Suc9=KValue(28,1);  KmSucP_Suc9=KValue(28,2);  Ke_Suc9=KValue(28,3);
KmSuc_Suc10=KValue(29,1);
KiADP_Suc3=KValue(30,1);  KIDHAP_Suc3=KValue(30,2);  KmATP_Suc3=KValue(30,3);  KmF26BP_Suc3=KValue(30,4);  KmF6P_Suc3=KValue(30,5);  Ke_Suc3=KValue(30,6);
KiF6P_Suc4=KValue(31,1);  KiPi_Suc4=KValue(31,2);  KmF26BP_Suc4=KValue(31,3);
KePi=KValue(36,1);
KmADP_ATPM=KValue(32,1);  KmATP_ATPM=KValue(32,2);  KmPi_ATPM=KValue(32,3);  X=KValue(32,4);  Y=KValue(32,5);  F=KValue(32,6);  Q=KValue(32,7);  D=KValue(32,8); Ke_ATPM=KValue(32,9);
KmNADP_NADPHM=KValue(33,1);  KmNADPH_NADPHM=KValue(33,2); Ke_NADPHM=KValue(33,3);  E=KValue(33,4);
KmADP_ATPB=KValue(34,1);  KmPi_ATPB=KValue(34,2);   KmATP_ATPB=KValue(34,3);  Ke_ATPB=KValue(34,4);   G=KValue(34,5);
KmNADP_NADPHB=KValue(37,1);  KmNADPH_NADPHB=KValue(37,2); Ke_NADPHB=KValue(37,3);
Voaa=KValue(35,1);  Vmal=KValue(35,2);  Vpyr=KValue(35,3);  Vpep=KValue(35,4);  Vt=KValue(35,5);  Vleak=KValue(35,6); Vpga=KValue(35,7);
KmCO2_PR1=KValue(38,1); KmO2_PR1=KValue(38,2);  KmRuBP_PR1=KValue(38,3);  KiPGA_PR1=KValue(38,4);  KiFBP_PR1=KValue(38,5);  KiSBP_PR1=KValue(38,6);  KiPi_PR1=KValue(38,7);  KiNADPH_PR1=KValue(38,8);
KmPGCA_PR2=KValue(39,1);  KiPI_PR2=KValue(39,2);  KiGCA_PR2=KValue(39,3);
KmGCA_PR3=KValue(40,1);
Ke_PS4=KValue(41,1);  KmGOA_PS4=KValue(41,2);  KmGLU_PS4=KValue(41,3);  KiGLY_PS4=KValue(41,4);
KmGLY_PS5=KValue(42,1);  KiSER_PS5=KValue(42,2);
Ke_PR6=KValue(43,1);  KmGOA_PR6=KValue(43,2);  KmSER_PR6=KValue(43,3);  KmGLY_PR6=KValue(43,4);
Ke_PR7=KValue(44,1);  KiHPR_PR7=KValue(44,2);  KmHPR_PR7=KValue(44,3);
Ke_PR8=KValue(45,1);  KmATP_PR8=KValue(45,2);  KmGCEA_PR8=KValue(45,3);  KiPGA_PR8=KValue(45,4);
KmGCA_PR9=KValue(46,1);  KiGCEA_PR9=KValue(46,2);
KmGCEA_PR10=KValue(47,1);  KiGCA_PR10=KValue(47,2);
KmPGA_62=KValue(48,1); KmPEP_62=KValue(48,2); Ke_62=KValue(48,3);




MC_HCO3= s(1);
MC_OAA=s(2);
MC_PEP=s(3);
MC_malate=s(4);
MC_pyruvate=s(5);
MC_PGA=s(6);
MC_FBP=s(7);
MC_UDPG=s(8);
MC_SUCP=s(9);
MC_SUC=s(10);
MC_F26BP=s(11);
MC_ATP=s(12);
MC_T3P=s(13);
MC_HexP=s(14);
MC_Sucrose=s(15);
Mchl_OAA= s(16);
Mchl_malate =s(17);
Mchl_PEP =s(18);
Mchl_pyruvate= s(19);
Mchl_NADPH= s(20);
Mchl_ATP= s(21);
Mchl_PGA= s(22);
Mchl_DPGA= s(23);
Mchl_T3P= s(24);
BSC_T3P= s(25);
BSC_PGA= s(26);
BSC_malate= s(27);
BSC_pyruvate= s(28);
BSC_CO2=s(29);
Bchl_CO2= s(30);
Bchl_RuBP= s(31);
Bchl_PGA= s(32);
Bchl_DPGA= s(33);
Bchl_ATP=s(34);
Bchl_NADPH= s(35);
Bchl_SBP= s(36);
Bchl_S7P= s(37);
Bchl_FBP= s(38);
Bchl_E4P= s(39);
Bchl_Starch= s(40);
Bchl_Rubisco= s(41);
Bchl_T3P= s(42);
Bchl_HexP= s(43);
Bchl_Pent =s(44);
Bchl_malate= s(45);
Bchl_pyruvate= s(46);

Bchl_PGCA=s(47);
Bchl_GCA=s(48);
Bchl_GCEA=s(49);

Bper_GCA=s(50);
Bper_GOA=s(51);
Bper_GLY=s(52);
Bper_SER=s(53);
Bper_HPR=s(54);
Bper_GCEA=s(55);
MC_CO2=s(56);

Bchl_PPi=s(57);
Bchl_ADPG=s(58);

MC_Glu=s(59);
MC_OxoG=s(60);
MC_Asp=s(61);
MC_Ala=s(62);
BSC_OxoG=s(63);
BSC_Glu=s(64);
BSC_Asp=s(65);
BSC_Ala=s(66);
BSC_OAA=s(67);
BSC_PEP=s(68);
BSC_ATP=s(69);
Bchl_OAA=s(70);

MC_O2=s(71);
Mchl_O2=s(72);
BSC_O2=s(73);
Bchl_O2=s(74);



global Velocity_s;
 
Vm_1=Velocity_s(1);
Vm_2=Velocity_s(2);          
Vm_3=Velocity_s(3);
Vm_4=Velocity_s(4);        
Vm_5=Velocity_s(5);
Vm_6=Velocity_s(6);
Vm_78=Velocity_s(7);
Vm_8=Velocity_s(8);
Vm_10=Velocity_s(9);
Vm_11=Velocity_s(10);
Vm_12=Velocity_s(11);
Vm_13=Velocity_s(12);
Vm_14=Velocity_s(13);
Vm_15=Velocity_s(14);
Vm_18=Velocity_s(15);
Vm_78Mchl=Velocity_s(16);
Vm_8Mchl=Velocity_s(17);
Vm_Starch=Velocity_s(18);
Vm_PGASink=Velocity_s(19);
Vm_Suc1=Velocity_s(20);
Vm_Suc2=Velocity_s(21);
Vm_Suc7=Velocity_s(22);
Vm_Suc8=Velocity_s(23);
Vm_Suc9=Velocity_s(24);
Vm_Suc10=Velocity_s(25);
Vm_Suc3=Velocity_s(26);
Vm_Suc4=Velocity_s(27);
I=Velocity_s(28);
Jmax=Velocity_s(29);
Vm_ATPM=Velocity_s(30);
Vm_NADPHM=Velocity_s(31); 
Vm_ATPB=Velocity_s(32);
Vm_NADPHB=Velocity_s(33);
Vm_PR1=Velocity_s(34);
Vm_PR2=Velocity_s(35);
Vm_PR3=Velocity_s(36);
Vm_PR4=Velocity_s(37);
Vm_PR5=Velocity_s(38);
Vm_PR6=Velocity_s(39);
Vm_PR7=Velocity_s(40);
Vm_PR8=Velocity_s(41);
VTgca_PR9=Velocity_s(42);
VTgcea_PR10=Velocity_s(43);
Vm_62=Velocity_s(44);
Vtp_Bchl=Velocity_s(45);
Vtp_Mchl=Velocity_s(46);
Vm_Sta1=Velocity_s(47);
Vm_Sta2=Velocity_s(48);
Vm_Sta3=Velocity_s(49);
Vm_OAA_M=Velocity_s(50);
Vm_PYR_B=Velocity_s(51);
Vm_PYR_M=Velocity_s(52);
Vm_PEP_M=Velocity_s(53);
Pmal=Velocity_s(54);
Ppyr=Velocity_s(55);
Pco2=Velocity_s(56);
PC3P=Velocity_s(57);
Pco2_B=Velocity_s(58);
Vm_MAL_B=Velocity_s(59);
Vm_MAL_M=Velocity_s(60);


global Bchl_CA;%assume constant concentration
global Bchl_CN;
global Bchl_CP;
global MC_CU;
global MC_CA;
global MC_CP;
global MC_UTP;
global Mchl_CP;
global Mchl_CA;
global Mchl_CN;
Mchl_NADP = Mchl_CN-Mchl_NADPH;
Mchl_Pi = Mchl_CP-Mchl_PGA-2*Mchl_DPGA-Mchl_T3P-Mchl_ATP-Mchl_PEP;
Mchl_GAP = Ke_9*Mchl_T3P/(1+Ke_9);
Mchl_DHAP = Mchl_T3P/(1+Ke_9);
Mchl_ADP = Mchl_CA-Mchl_ATP;

MC_UDP = MC_CU-MC_UTP-MC_UDPG;
MC_PiT = MC_CP-2*MC_FBP-2*MC_F26BP-MC_PGA-MC_T3P-MC_HexP-MC_SUCP-MC_UTP-MC_ATP-MC_PEP;%1
MC_Pi = (sqrt(KePi^2+4*KePi*MC_PiT)-KePi)/2;
MC_PPi = MC_PiT-MC_Pi;
MC_GAP = Ke_9*MC_T3P/(1+Ke_9);
MC_DHAP = MC_T3P/(1+Ke_9);
MC_G6P = MC_HexP/(1/Ke_Suc5+Ke_Suc6+1);
MC_G1P = Ke_Suc6*MC_HexP/(1/Ke_Suc5+Ke_Suc6+1);
MC_F6P = (MC_HexP/Ke_Suc5)/(1/Ke_Suc5+Ke_Suc6+1);
MC_ADP = MC_CA-MC_ATP;

BSC_GAP = Ke_9*BSC_T3P/(1+Ke_9);
BSC_DHAP = BSC_T3P/(1+Ke_9);

Bchl_NADP = Bchl_CN-Bchl_NADPH;
Bchl_Xu5P = (Bchl_Pent/Ke_17)/(1/Ke_16+1/Ke_17+1);
Bchl_Ru5P = Bchl_Pent/(1/Ke_16+1/Ke_17+1);
Bchl_Ri5P = (Bchl_Pent/Ke_16)/(1/Ke_16+1/Ke_17+1);

Bchl_Pi = Bchl_CP-Bchl_PGA-2*Bchl_DPGA-Bchl_T3P-2*Bchl_FBP-Bchl_HexP-Bchl_E4P-2*Bchl_SBP-Bchl_S7P-Bchl_Pent-2*Bchl_RuBP-Bchl_ATP-Bchl_PGCA-Bchl_PPi;
global Pim;
Pim=Bchl_Pi;
Bchl_GAP = Ke_9*Bchl_T3P/(1+Ke_9);
Bchl_DHAP = Bchl_T3P/(1+Ke_9);
Bchl_G6P = Bchl_HexP/(1/Ke_Starch1+Ke_Starch2+1);
Bchl_G1P = Ke_Starch2*Bchl_HexP/(1/Ke_Starch1+Ke_Starch2+1);
Bchl_F6P = (Bchl_HexP/Ke_Starch1)/(1/Ke_Starch1+Ke_Starch2+1);
Bchl_ADP = Bchl_CA-Bchl_ATP-Bchl_ADPG;

global CI;
global Bper_GLU;
global Bper_KG;
global Bper_NADH;
global Bper_NAD;
global U;
global V;

global v1;
global v2;
global v3;
global v4;
global v5;
global v6;
global v78;
global v10;
global v11;
global v12;
global v13;
global v14;
global v15;
global v18;
global v78Mchl;
global vSta1;
global vStarch;
global vHexP;
global vpgasink;
global vSuc1;
global vOAA_M;
global vMAL_M;
global vMAL;
global vMAL_B;
global vPYR;

global vPGA_B;
global vDHAP_B;
global vGAP_B;
global vPGA;
global vDHAP;
global vGAP;
global vPGA_M;
global vDHAP_M;
global vGAP_M;
global vleakage;
global vleakage2;
global vpr1;
global vAsp;
global vtO2;
global vtO2_B;
global vtO2_M;
global enzyme_flux;


global Ratio;

%%parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%vinf
gm=0.7;%molm-2 bar-1 
Sc=3*10^4;%;%ubarL/mmol
%v2%
KAG6P_2=0.29;%
KAT3P_2=3;%
%v78
KmPGA_78=5;
KmATP_78=0.3;
KmNADPH_78=0.1;
KmADP_78=0.5;
KmNADP_78=0.5;
%v78Mchl
KmPGA_78Mchl=5;
KmATP_78Mchl=0.3;
KmNADPH_78Mchl=0.1;
KmADP_78Mchl=0.5;
KmNADP_78Mchl=0.5;
%vSta1  2.7.7.27
KaPGA_Sta1=0.2;%
KmG1P_Sta1=0.06;%
KmATP_Sta1=0.12;%
KIAPi_ATP_Sta1=2.96;
KmPPi_Sta1=0.033;
KICPP1_ATP_Sta1=13.8E-4;
KmADPG_Sta1=0.24;
KIAADP_ATP_Sta1=2.0;
Ke_Sta1=1.1;
%vSta2  3.6.1.1
KmPPi_Sta2=0.154;
Ke_Sta2=15700.0;
%vSta3 %  2.4.1.21
KmADPG_Sta3=0.077;
%vhexp
Kmpi_hexp=1.5;
Kmhexp_hexp=1;
Vm_hexp=0.0005;
%Transport
Km_OAA_M=0.053;
Kimal_OAA_M=7.5;
Km_MAL_M=0.5;
KiOAA_MAL_M=0.3;%0.3;
Km_MAL_B=1;
Km_PYR_B=0.1;
Km_PYR_M=0.1;
Km_PEP_M=0.5;
KmPGA_B = 2; 
KmGAP_B =2; 
KmDHAP_B =2; 
KmPGA =2; 
KmGAP =2; 
KmDHAP = 2; 

%%%%%%%%%%
%PCK
%%%%%%%%%%
%2.6.1.1M   PCK1
KmAsp_PCK1=2.5;KmOxog_PCK1=0.14;KmGlu_PCK1=17;KmOAA_PCK1=0.056; Ke_PCK1=1/0.148;
%2.6.1.1B   PCK2
KmAsp_PCK2=2.5;KmOxog_PCK2=0.14;KmGlu_PCK2=17;KmOAA_PCK2=0.056; Ke_PCK2=0.148;
%4.1.1.49   PCK3
KmOAA_PCK3=0.06; KmATP_PCK3=0.034*3;
%2.6.1.2B   PCK4
KmPyr_PCK4=0.33;KmGlu_PCK4=5;KmAla_PCK4=6.67;KmOxog_PCK4=0.15; Ke_PCK4=1;
%2.6.1.2M   PCK5
KmPyr_PCK5=0.33;KmGlu_PCK5=5;KmAla_PCK5=6.67;KmOxog_PCK5=0.15; Ke_PCK5=1;
%1.1.1.82B  PCK6
KmNADPH_PCK6 =0.05;  KmOAA_PCK6 =0.056;  KmNADP_PCK6 =0.045;  Kmmal_PCK6 =32.0;  Ke_PCK6 =4450.0; % No unit

Vm_PCK1=0.075*Ratio;
Vm_PCK2=0.075*Ratio;
Vm_PCK4=0.075*Ratio;
Vm_PCK5=0.075*Ratio;
Vm_PCK6=0.02*Ratio; 
Vm_PCK3=0.0075*Ratio;

PAsp=0.0332;
PAla=0.0436;
PPEP=0.0327;
POAA_B=1;

VtATP_B=3;
VtATP=5;

%O2 diffusion
PO2=0.122633;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rate equations of each enzyme and transport process
%%CO2 from intercellular space to mesophyll cell
vinf=gm*Sc*10^(-3)*(CI-MC_CO2);

%C4 Cycle 5
v1=Vm_1*(MC_CO2-MC_HCO3/Ke_1)/(KmCO2_1+MC_CO2);
v2=Vm_2*MC_HCO3*MC_PEP/(MC_PEP+KmPEP_2*(1+MC_malate/Kimal_2)/(1+MC_G6P/KAG6P_2+MC_T3P/KAT3P_2))/(MC_HCO3+KmHCO3_2);
v3=Vm_3*Mchl_OAA*Mchl_NADPH/(KmOAA_3+Mchl_OAA)/(KmNADPH_3+Mchl_NADPH);
v4=Vm_4*(Bchl_malate*Bchl_NADP-Bchl_pyruvate*Bchl_NADPH*Bchl_CO2/Ke_4)/(Kmmal_4*KmNADP_4)/(1+Bchl_malate/Kmmal_4+Bchl_NADP/KmNADP_4+Bchl_pyruvate/KmPyr_4+Bchl_NADPH/KmNADPH_4+Bchl_CO2/KmCO2_4+Bchl_malate*Bchl_NADP/(Kmmal_4*KmNADP_4)+Bchl_pyruvate*Bchl_NADPH/(KmPyr_4*KmNADPH_4)+Bchl_pyruvate*Bchl_CO2/(KmPyr_4*KmCO2_4)+Bchl_NADPH*Bchl_CO2/(KmNADPH_4*KmCO2_4)+Bchl_pyruvate*Bchl_NADPH*Bchl_CO2/(KmPyr_4*KmNADPH_4*KmCO2_4));
v5=Vm_5*Mchl_pyruvate*Mchl_ATP/(Mchl_pyruvate+KmPyr_5*(1+Mchl_PEP/KiPEP_5))/(Mchl_ATP+KmATP_5);

% Calvin cycle 12
v6=Vm_6*Bchl_RuBP*Bchl_CO2/((Bchl_CO2+KmCO2_6*(1+Bchl_O2/KmO2_6))*(Bchl_RuBP+KmRuBP_6*(1+Bchl_PGA/KiPGA_6+Bchl_FBP/KiFBP_6+Bchl_SBP/KiSBP_6+Bchl_Pi/ KiPi_6)));
v7=0;%not used
v8=0;%not used
v78=Vm_78*(Bchl_PGA*Bchl_ATP*Bchl_NADPH)/((Bchl_PGA+KmPGA_78)*(Bchl_ATP+KmATP_78)*(Bchl_NADPH+KmNADPH_78));
v10=Vm_10*(Bchl_GAP*Bchl_DHAP-Bchl_FBP/Ke_10)/(KmGAP_10*KmDHAP_10*(1+Bchl_GAP/KmGAP_10+Bchl_DHAP/KmDHAP_10+Bchl_FBP/KmFBP_10+Bchl_GAP*Bchl_DHAP/(KmGAP_10*KmDHAP_10)));
v11=Vm_11*(Bchl_FBP-Bchl_F6P*Bchl_Pi/Ke_11)/(Bchl_FBP+KmFBP_11*(1+Bchl_F6P/KiF6P_11+Bchl_Pi/KiPi_11));
v12=Vm_12*(Bchl_DHAP*Bchl_E4P-Bchl_SBP/Ke_12)/((Bchl_E4P+KmE4P_12)*(Bchl_DHAP+KmDHAP_12));
v13=Vm_13*(Bchl_SBP-Bchl_Pi*Bchl_S7P/Ke_13)/(Bchl_SBP+KmSBP_13*(1+Bchl_Pi/KiPi_13));
v14=Vm_14*(Bchl_F6P*Bchl_GAP-Bchl_Xu5P*Bchl_E4P/Ke_14)/((Bchl_F6P + KmF6P_14*(1+Bchl_Xu5P/KmXu5P+Bchl_E4P/KmE4P_14))*(Bchl_GAP+KmGAP_14));
v15=Vm_15*(Bchl_GAP*Bchl_S7P-Bchl_Ri5P*Bchl_Xu5P/Ke_15)/((Bchl_GAP+KmGAP_15*(1+Bchl_Xu5P/KmXu5P_15+Bchl_Ri5P/KmRi5P_15))*(Bchl_S7P+KmS7P_15));
v18=Vm_18*(Bchl_ATP*Bchl_Ru5P-Bchl_ADP*Bchl_RuBP/Ke_18)/((Bchl_ATP*(1+Bchl_ADP/KiADP_18)+KmATP_18*(1+Bchl_ADP/Ki_ADP_18))*(Bchl_Ru5P+KmRu5P_18*(1+Bchl_PGA/KiPGA_18+Bchl_RuBP/KiRuBP_18+Bchl_Pi/KiPi_18)));
v7Mchl=0;%not used
v8Mchl=0;%not used
v78Mchl=Vm_78Mchl*(Mchl_PGA*Mchl_ATP*Mchl_NADPH)/((Mchl_PGA+KmPGA_78Mchl)*(Mchl_ATP+KmATP_78Mchl)*(Mchl_NADPH+KmNADPH_78Mchl));
vStarch1=0;%not used
vStarch2=0;%not used
Vm_Sta1=Vm_Sta1*Bchl_PGA/(Bchl_PGA+KaPGA_Sta1);
vSta1=Vm_Sta1*(Bchl_G1P*Bchl_ATP-Bchl_ADPG*Bchl_PPi/Ke_Sta1)/(KmG1P_Sta1*KmATP_Sta1*(1+Bchl_ADP/KIAADP_ATP_Sta1+Bchl_PPi/KICPP1_ATP_Sta1+Bchl_Pi/KIAPi_ATP_Sta1)*(1+Bchl_G1P/KmG1P_Sta1+Bchl_ATP*(1+Bchl_Pi/KIAPi_ATP_Sta1+Bchl_ADP/KIAADP_ATP_Sta1)/(KmATP_Sta1*(1+Bchl_ADP/KIAADP_ATP_Sta1+Bchl_PPi/KICPP1_ATP_Sta1+Bchl_Pi/KIAPi_ATP_Sta1))+Bchl_ADPG/KmADPG_Sta1+Bchl_PPi/KmPPi_Sta1+Bchl_G1P*Bchl_ATP*(1+Bchl_Pi/KIAPi_ATP_Sta1+Bchl_ADP/KIAADP_ATP_Sta1)/(KmG1P_Sta1*KmATP_Sta1*(1+Bchl_ADP/KIAADP_ATP_Sta1+Bchl_PPi/KICPP1_ATP_Sta1+Bchl_Pi/KIAPi_ATP_Sta1))+Bchl_ADPG*Bchl_PPi/(KmADPG_Sta1*KmPPi_Sta1)));
vSta2=Vm_Sta2*(Bchl_PPi-Bchl_Pi*Bchl_Pi/Ke_Sta2)/(Bchl_PPi+KmPPi_Sta2);
vSta3=Vm_Sta3*Bchl_ADPG/(Bchl_ADPG+KmADPG_Sta3);
vStarch=vSta3;
%Starch degradation
vhexp=Vm_hexp*(Bchl_Pi/(Kmpi_hexp*(1+Bchl_HexP/Kmhexp_hexp)+Bchl_Pi));
%vhexp=0;
vHexP=vhexp;

vPGASink=Vm_PGASink*MC_PGA/(MC_PGA+KmPGA_PGASink);%MC
vpgasink=vPGASink;

vSuc1=Vm_Suc1*(MC_GAP*MC_DHAP-MC_FBP/Ke_Suc1)/(KmGAP_Suc1*KmDHAP_Suc1*(1+MC_GAP/KmGAP_Suc1+MC_DHAP/KmDHAP_Suc1+MC_FBP/KmFBP_Suc1+MC_GAP*MC_DHAP/(KmGAP_Suc1*KmDHAP_Suc1)));
vSuc2=Vm_Suc2*(MC_FBP-MC_F6P*MC_Pi/Ke_Suc2)/(KmFBP_Suc2*(1+MC_F26BP/KiF26BP_Suc2)*(1+MC_FBP/(KmFBP_Suc2*(1+MC_F26BP/KiF26BP_Suc2))+MC_Pi/KiPi_Suc2+MC_F6P/KiF6P_Suc2+MC_Pi*MC_F6P/(KiPi_Suc2*KiF6P_Suc2)));
vSuc7=Vm_Suc7*(MC_UTP*MC_G1P-MC_UDPG*MC_PPi/Ke_Suc7)/(KmUTP_Suc7*KmG1P_Suc7*(1+MC_UTP/KmUTP_Suc7+MC_G1P/KmG1P_Suc7+MC_UDPG/KmUDPG_Suc7+MC_PPi/KmPPi_Suc7+MC_UTP*MC_G1P/(KmUTP_Suc7*KmG1P_Suc7)+MC_UDPG*MC_PPi/(KmUDPG_Suc7*KmPPi_Suc7)));
vSuc8=Vm_Suc8*(MC_F6P*MC_UDPG-MC_SUCP*MC_UDP/Ke_Suc8)/((MC_F6P+KmF6P_Suc8*(1+MC_FBP/KiFBP_Suc8))*(MC_UDPG+KmUDPG_Suc8*(1+MC_UDP/KiUDP_Suc8)*(1+MC_SUCP/KiSucP_Suc8)*(1+MC_SUC/KiSuc_Suc8)*(1+MC_Pi/KiPi_Suc8)));
vSuc9=Vm_Suc9*(MC_SUCP-MC_SUC*MC_Pi/Ke_Suc9)/(MC_SUCP+KmSucP_Suc9*(1+MC_SUC/KmSuc_Suc9));
vSuc10=Vm_Suc10*MC_SUC/(MC_SUC+KmSuc_Suc10);
vSuc3=Vm_Suc3*(MC_ATP*MC_F6P-MC_ADP*MC_F26BP/Ke_Suc3)/((MC_F6P+KmF6P_Suc3*(1+MC_F26BP/KmF26BP_Suc3)*(1+MC_DHAP/KIDHAP_Suc3))*(MC_ATP+KmATP_Suc3*(1+MC_ADP/KiADP_Suc3)));
vSuc4=Vm_Suc4*MC_F26BP/(KmF26BP_Suc4*(1+MC_F26BP/KmF26BP_Suc4)*(1+MC_Pi/KiPi_Suc4)*(1+MC_F6P/KiF6P_Suc4));

%ATP&NADPH 3
ETRa=D*((F/2*X*I+Y*Jmax-sqrt((F/2*X*I+Y*Jmax)^2-4*Q*F/2*X*I*Y*Jmax))/(2*Q));
ETRn=E*((F/2*X*I+Y*Jmax-sqrt((F/2*X*I+Y*Jmax)^2-4*Q*F/2*X*I*Y*Jmax))/(2*Q));
vATPM=min(Vm_ATPM,ETRa)*(Mchl_ADP*Mchl_Pi-Mchl_ATP/(Ke_ATPM))/(KmADP_ATPM*KmPi_ATPM*(1+Mchl_ADP/KmADP_ATPM+Mchl_Pi/KmPi_ATPM+Mchl_ATP/KmATP_ATPM+Mchl_ADP*Mchl_Pi/(KmADP_ATPM*KmPi_ATPM)));
vNADPHM=min(Vm_NADPHM,ETRn)*(Mchl_NADP-Mchl_NADPH/Ke_NADPHM)/(KmNADP_NADPHM*(1+Mchl_NADP/KmNADP_NADPHM+Mchl_NADPH/KmNADPH_NADPHM));
vO2_Mchl=vNADPHM/2;

ETRab=G*((F*(1-U)*(1-X)*I+(1-V)*(1-Y)*Jmax-sqrt((F*(1-U)*(1-X)*I+(1-V)*(1-Y)*Jmax)^2-4*Q*F*(1-U)*(1-X)*I*(1-V)*(1-Y)*Jmax))/(2*Q));
ETRabl=D*((F/2*U*(1-X)*I+V*(1-Y)*Jmax-sqrt((F/2*U*(1-X)*I+V*(1-Y)*Jmax)^2-4*Q*F/2*U*(1-X)*I*V*(1-Y)*Jmax))/(2*Q));
ETRnbl=E*((F/2*U*(1-X)*I+V*(1-Y)*Jmax-sqrt((F/2*U*(1-X)*I+V*(1-Y)*Jmax)^2-4*Q*F/2*U*(1-X)*I*V*(1-Y)*Jmax))/(2*Q));

vATPB=min(Vm_ATPB,ETRab+ETRabl)*(Bchl_ADP*Bchl_Pi-Bchl_ATP/Ke_ATPB)/(KmADP_ATPB*KmPi_ATPB*(1+Bchl_ADP/KmADP_ATPB+Bchl_Pi/KmPi_ATPB+Bchl_ATP/KmATP_ATPB+Bchl_ADP*Bchl_Pi/(KmADP_ATPB*KmPi_ATPB)));
vNADPHB=min(Vm_NADPHB,ETRnbl)*(Bchl_NADP-Bchl_NADPH/Ke_NADPHB)/(KmNADP_NADPHB*(1+Bchl_NADP/KmNADP_NADPHB+Bchl_NADPH/KmNADPH_NADPHB));
vO2_Bchl=vNADPHB/2;

%Transport 17
vOAA_M=Vm_OAA_M*(MC_OAA-Mchl_OAA)/(MC_OAA+Km_OAA_M*(1+MC_malate/Kimal_OAA_M));
vMAL_M=Vm_MAL_M*(Mchl_malate-MC_malate)/(Mchl_malate+Km_MAL_M*(1+Mchl_OAA/KiOAA_MAL_M));% 2  Mchl_malate -> MC.malate
vMAL=Pmal*(MC_malate-BSC_malate);% 3 MC_malate -> BSC_malate
vMAL_B=Vm_MAL_B*(BSC_malate-Bchl_malate)/(BSC_malate+Km_MAL_B);% 4 BSC_malate -> Bchl_malate
vPYR_B=Vm_PYR_B*(Bchl_pyruvate-BSC_pyruvate/10)/(Km_PYR_B+Bchl_pyruvate);
vPYR=Ppyr*(BSC_pyruvate-MC_pyruvate);% 6 BSC_pyruvate -> MC_pyruvate

vPYR_M=Vm_PYR_M*(MC_pyruvate-Mchl_pyruvate/10)/(Km_PYR_M+MC_pyruvate);
vPEP_M=Vm_PEP_M*(Mchl_PEP-MC_PEP)/(Km_PEP_M+Mchl_PEP); % 8 Mchl_PEP + MC_Pi -> MC_PEP + Mchl_Pi

vPGA_B = Vtp_Bchl * Bchl_PGA/(Bchl_PGA + KmPGA_B * ( 1 + Bchl_DHAP/KmDHAP_B) * ( 1 + Bchl_GAP/KmGAP_B))-Vtp_Bchl * BSC_PGA/(BSC_PGA + KmPGA_B * ( 1 + BSC_DHAP/KmDHAP_B) * ( 1 + BSC_GAP/KmGAP_B));  
vGAP_B = Vtp_Bchl * BSC_GAP/(BSC_GAP + KmGAP_B * ( 1 + BSC_PGA/KmPGA_B) * ( 1 + BSC_DHAP/KmDHAP_B))-Vtp_Bchl * Bchl_GAP/(Bchl_GAP + KmGAP_B * ( 1 + Bchl_PGA/KmPGA_B) * ( 1 + Bchl_DHAP/KmDHAP_B));
vDHAP_B= Vtp_Bchl * BSC_DHAP/(BSC_DHAP + KmDHAP_B * ( 1 + BSC_PGA/KmPGA_B) * ( 1 + BSC_GAP/KmGAP_B))- Vtp_Bchl * Bchl_DHAP/(Bchl_DHAP + KmDHAP_B * ( 1 + Bchl_PGA/KmPGA_B) * ( 1 + Bchl_GAP/KmGAP_B));

vPGA=PC3P*(BSC_PGA-MC_PGA); % 10 BSC_PGA+MC.Pi -> MC.PGA+BSC_Pi
vGAP=PC3P*(MC_GAP-BSC_GAP); % 13 MC_GAP+BSC_Pi ->BSC_GAP+MC_Pi
vDHAP=PC3P*(MC_DHAP-BSC_DHAP); % 16 MC_DHAP+BSC_Pi ->BSC_DHAP+MC_Pi
vPGA_M=Vtp_Mchl * MC_PGA/(MC_PGA + KmPGA * ( 1 + MC_DHAP/KmDHAP) * ( 1 + MC_GAP/KmGAP))-Vtp_Mchl * Mchl_PGA/(Mchl_PGA + KmPGA * ( 1 + Mchl_DHAP/KmDHAP) * ( 1 + Mchl_GAP/KmGAP));
vGAP_M=Vtp_Mchl * Mchl_GAP/(Mchl_GAP + KmGAP * ( 1 + Mchl_PGA/KmPGA) * ( 1 + Mchl_DHAP/KmDHAP))-Vtp_Mchl * MC_GAP/(MC_GAP + KmGAP * ( 1 + MC_PGA/KmPGA) * ( 1 + MC_DHAP/KmDHAP));
vDHAP_M=Vtp_Mchl * Mchl_DHAP/(Mchl_DHAP + KmDHAP * ( 1 + Mchl_PGA/KmPGA) * ( 1 + Mchl_GAP/KmGAP))-Vtp_Mchl * MC_DHAP/(MC_DHAP + KmDHAP * ( 1 + MC_PGA/KmPGA) * ( 1 + MC_GAP/KmGAP));

vtATP=VtATP*(Mchl_ATP-MC_ATP*2);

vleak_B=Pco2_B*(Bchl_CO2-BSC_CO2);%0.5376*10*(Bchl_CO2-BSC_CO2);%Vleak*(Bchl_CO2-BSC_CO2);%  Bchl_CO2 -> BSC_CO2
vleak=Pco2*(BSC_CO2-MC_CO2);%Vleak*(BSC_CO2-MC_CO2);%BSC_CO2->MC_CO2
vleakage=vleak;
vleakage2=vleak_B;

% Photorespiration

vpr1=Vm_PR1*Bchl_RuBP*Bchl_O2/((Bchl_O2+KmO2_PR1*(1+Bchl_CO2/KmCO2_PR1))*(Bchl_RuBP+KmRuBP_PR1*(1+Bchl_PGA/KiPGA_PR1+Bchl_FBP/KiFBP_PR1+Bchl_SBP/KiSBP_PR1+Bchl_Pi/KiPi_PR1+Bchl_NADPH/KiNADPH_PR1)));
vpr2=Vm_PR2*Bchl_PGCA/(Bchl_PGCA+KmPGCA_PR2*(1+Bchl_GCA/KiGCA_PR2)*(1+Bchl_Pi/KiPI_PR2));
vpr3=Vm_PR3*Bper_GCA/(Bper_GCA+KmGCA_PR3);
vpr4=Vm_PR4*(Bper_GOA*Bper_GLU-Bper_KG*Bper_GLY/Ke_PS4)/((Bper_GOA+KmGOA_PS4)*(Bper_GLU+KmGLU_PS4*(1+Bper_GLY/KiGLY_PS4)));
vpr5=Vm_PR5*Bper_GLY/(Bper_GLY+KmGLY_PS5*(1+Bper_SER/KiSER_PS5));
vpr6=Vm_PR6*(Bper_GOA*Bper_SER-Bper_HPR*Bper_GLY/Ke_PR6)/((Bper_GOA+KmGOA_PR6)*(Bper_SER+KmSER_PR6*(1+Bper_GLY/KmGLY_PR6)));
vpr7=Vm_PR7*(Bper_HPR*Bper_NADH-Bper_NAD*Bper_GCEA/Ke_PR7)/(Bper_HPR+KmHPR_PR7*(1+Bper_HPR/KiHPR_PR7));
vpr8=Vm_PR8*(Bchl_ATP*Bchl_GCEA-Bchl_ADP*Bchl_PGA/Ke_PR8)/((Bchl_ATP+KmATP_PR8*(1+Bchl_PGA/KiPGA_PR8))*(Bchl_GCEA+KmGCEA_PR8));
vpr9=VTgca_PR9*(Bchl_GCA/(Bchl_GCA+KmGCA_PR9*(1+Bchl_GCEA/KiGCEA_PR9))-Bper_GCA/(Bper_GCA+KmGCA_PR9*(1+Bper_GCEA/KiGCEA_PR9)));
vpr10=VTgcea_PR10*(Bper_GCEA/(Bper_GCEA+KmGCEA_PR10*(1+Bper_GCA/KiGCA_PR10))-Bchl_GCEA/(Bchl_GCEA+KmGCEA_PR10*(1+Bchl_GCA/KiGCA_PR10)));

%PGA enolase, phosphoglyceromutase in MC PGA->PEP
v62=Vm_62*(MC_PGA-MC_PEP/Ke_62)/(KmPGA_62*(1+MC_PGA/KmPGA_62+MC_PEP/KmPEP_62));

vgly1B=0;%not used 

%PCK
vPCK1=Vm_PCK1*(MC_OAA*MC_Glu-MC_OxoG*MC_Asp/Ke_PCK1)/(KmOAA_PCK1*KmGlu_PCK1*(1+MC_OAA/KmOAA_PCK1+MC_Glu/KmGlu_PCK1+MC_OxoG/KmOxog_PCK1+MC_Asp/KmAsp_PCK1+MC_OAA/KmOAA_PCK1*MC_Glu/KmGlu_PCK1+MC_OxoG/KmOxog_PCK1*MC_Asp/KmAsp_PCK1));
vPCK2=Vm_PCK2*(BSC_OxoG*BSC_Asp-BSC_OAA*BSC_Glu/Ke_PCK2)/(KmOxog_PCK2*KmAsp_PCK2*(1+BSC_OAA/KmOAA_PCK2+BSC_Glu/KmGlu_PCK2+BSC_OxoG/KmOxog_PCK2+BSC_Asp/KmAsp_PCK2+BSC_OAA/KmOAA_PCK2*BSC_Glu/KmGlu_PCK2+BSC_OxoG/KmOxog_PCK2*BSC_Asp/KmAsp_PCK2));
vPCK3=Vm_PCK3*BSC_OAA*BSC_ATP/((BSC_OAA+KmOAA_PCK3)*(BSC_ATP+KmATP_PCK3));
vPCK4=Vm_PCK4*(BSC_Glu*BSC_pyruvate-BSC_Ala*BSC_OxoG/Ke_PCK4)/(KmGlu_PCK4*KmPyr_PCK4*(1+BSC_Glu/KmGlu_PCK4+BSC_pyruvate/KmPyr_PCK4+BSC_Ala/KmAla_PCK4+BSC_OxoG/KmOxog_PCK4+BSC_Glu/KmGlu_PCK4*BSC_pyruvate/KmPyr_PCK4+BSC_Ala/KmAla_PCK4*BSC_OxoG/KmOxog_PCK4));
vPCK5=Vm_PCK5*(MC_Ala*MC_OxoG-MC_Glu*MC_pyruvate/Ke_PCK5)/(KmAla_PCK5*KmOxog_PCK5*(1+MC_Glu/KmGlu_PCK5+MC_pyruvate/KmPyr_PCK5+MC_Ala/KmAla_PCK5+MC_OxoG/KmOxog_PCK5+MC_Glu/KmGlu_PCK5*MC_pyruvate/KmPyr_PCK5+MC_Ala/KmAla_PCK5*MC_OxoG/KmOxog_PCK5));
vPCK6=Vm_PCK6*(Bchl_OAA*Bchl_NADPH-Bchl_NADP*Bchl_malate/Ke_PCK6)/( KmOAA_PCK6* KmNADPH_PCK6*(1+Bchl_OAA/KmOAA_PCK6+ Bchl_NADPH/KmNADPH_PCK6+ Bchl_NADP/KmNADP_PCK6+ Bchl_malate/Kmmal_PCK6+ Bchl_OAA*Bchl_NADPH/(KmOAA_PCK6* KmNADPH_PCK6)+ Bchl_NADP*Bchl_malate/(KmNADP_PCK6* Kmmal_PCK6)));

%%%%%%%%%%%%%
%Transport for PCK
%%%%%%%%%%%%%
%Tasp
vAsp=PAsp*(MC_Asp-BSC_Asp);
%TAla
vAla=PAla*(BSC_Ala-MC_Ala);
%TPEP
vPEP=PPEP*(BSC_PEP-MC_PEP);
%TOAAB
vOAA_B=POAA_B*(BSC_OAA-Bchl_OAA);
%ATPB
vATP_B=VtATP_B*(Bchl_ATP-2*BSC_ATP);

%%%%%O2diffusion%%%%
vtO2=PO2*(BSC_O2-MC_O2);
vtO2_B=Pco2_B*(Bchl_O2-BSC_O2);
vtO2_M=Pco2_B*(Mchl_O2-MC_O2);


Enz_v=zeros(87,1);   
Enz_v(1)=v1;
Enz_v(2)=v2;
Enz_v(3)=v3;
Enz_v(4)=v4;
Enz_v(5)=v5;
Enz_v(6)=v6;
Enz_v(7)=v7;
Enz_v(8)=v8;
Enz_v(9)=v10;
Enz_v(10)=v11;
Enz_v(11)=v12;
Enz_v(12)=v13;
Enz_v(13)=v14;
Enz_v(14)=v15;
Enz_v(15)=v18;
Enz_v(16)=v7Mchl;
Enz_v(17)=v8Mchl;
Enz_v(18)=vStarch1;
Enz_v(19)=vPGASink;
Enz_v(20)=vSuc1;
Enz_v(21)=vSuc2;
Enz_v(22)=vSuc3;
Enz_v(23)=vSuc4;
Enz_v(24)=vSuc7;
Enz_v(25)=vSuc8;
Enz_v(26)=vSuc9;
Enz_v(27)=vSuc10;
Enz_v(28)=vATPM;
Enz_v(29)=vNADPHM;
Enz_v(30)=vATPB;
Enz_v(31)=vOAA_M;
Enz_v(32)=vMAL_M;
Enz_v(33)=vMAL;
Enz_v(34)=vMAL_B;
Enz_v(35)=vPYR_B;
Enz_v(36)=vPYR;
Enz_v(37)=vPYR_M;
Enz_v(38)=vPEP_M;
Enz_v(39)=vPGA_B;
Enz_v(40)=vPGA;
Enz_v(41)=vPGA_M;
Enz_v(42)=vGAP_M;
Enz_v(43)=vGAP;
Enz_v(44)=vGAP_B;
Enz_v(45)=vDHAP_M;
Enz_v(46)=vDHAP;
Enz_v(47)=vDHAP_B;
Enz_v(48)=vleak_B;
Enz_v(49)=vleak;
Enz_v(50)=vNADPHB;
Enz_v(51)=vpr1;
Enz_v(52)=vpr2;
Enz_v(53)=vpr3;
Enz_v(54)=vpr4;
Enz_v(55)=vpr5;
Enz_v(56)=vpr6;
Enz_v(57)=vpr7;
Enz_v(58)=vpr8;
Enz_v(59)=vpr9;
Enz_v(60)=vpr10;
Enz_v(61)=vinf;
Enz_v(62)=v62;
Enz_v(63)=v78;
Enz_v(64)=v78Mchl;
Enz_v(65)=vStarch2;
Enz_v(66)=vgly1B;
Enz_v(67)=vhexp;
Enz_v(68)=vSta1;
Enz_v(69)=vSta2;
Enz_v(70)=vSta3;
Enz_v(71)=vtATP;
Enz_v(72)=vPCK1;
Enz_v(73)=vPCK2;
Enz_v(74)=vPCK3;
Enz_v(75)=vPCK4;
Enz_v(76)=vPCK5;
Enz_v(77)=vPCK6;
Enz_v(78)=vAsp;
Enz_v(79)=vAla;
Enz_v(80)=vPEP;
Enz_v(81)=vOAA_B;
Enz_v(82)=vATP_B;
Enz_v(83)=vO2_Mchl;
Enz_v(84)=vO2_Bchl;
Enz_v(85)=vtO2;
Enz_v(86)=vtO2_B;
Enz_v(87)=vtO2_M;
enzyme_flux=Enz_v;
