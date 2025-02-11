% Copyright  Wang Y and Zhu X-G, 2014. 
% CAS-MPG Partner Institute for Computational Biology, Shanghai Institutes for Biological Sciences, CAS, Shanghai,200031 
% A mixed pathway model of C4 photosynthesis
% The code for differential equations

function EnzMB = C4MB(t,s)  
global vrpd; 
Enz_v=C4Vel(t,s);
v1=Enz_v(1);
v2=Enz_v(2);
v3=Enz_v(3);
v4=Enz_v(4);
v5=Enz_v(5);
v6=Enz_v(6);
v7=Enz_v(7);
v8=Enz_v(8);
v10=Enz_v(9);
v11=Enz_v(10);
v12=Enz_v(11);
v13=Enz_v(12);
v14=Enz_v(13);
v15=Enz_v(14);
v18=Enz_v(15);
v7Mchl=Enz_v(16);
v8Mchl=Enz_v(17);
vStarch1=Enz_v(18);
vPGASink=Enz_v(19);
vSuc1=Enz_v(20);
vSuc2=Enz_v(21);
vSuc3=Enz_v(22);
vSuc4=Enz_v(23);
vSuc7=Enz_v(24);
vSuc8=Enz_v(25);
vSuc9=Enz_v(26);
vSuc10=Enz_v(27);
vATPM=Enz_v(28);
vNADPHM=Enz_v(29);
vATPB=Enz_v(30);
vOAA_M=Enz_v(31);
vMAL_M=Enz_v(32);
vMAL=Enz_v(33);
vMAL_B=Enz_v(34);
vPYR_B=Enz_v(35);
vPYR=Enz_v(36);
vPYR_M=Enz_v(37);
vPEP_M=Enz_v(38);
vPGA_B=Enz_v(39);
vPGA=Enz_v(40);
vPGA_M=Enz_v(41);
vGAP_M=Enz_v(42);
vGAP=Enz_v(43);
vGAP_B=Enz_v(44);
vDHAP_M=Enz_v(45);
vDHAP=Enz_v(46);
vDHAP_B=Enz_v(47);
vleak_B=Enz_v(48);
vleak=Enz_v(49);
vNADPHB=Enz_v(50);
vpr1=Enz_v(51);
vpr2=Enz_v(52);
vpr3=Enz_v(53);
vpr4=Enz_v(54);
vpr5=Enz_v(55);
vpr6=Enz_v(56);
vpr7=Enz_v(57);
vpr8=Enz_v(58);
vpr9=Enz_v(59);
vpr10=Enz_v(60);
vinf=Enz_v(61);
v62=Enz_v(62);
v78=Enz_v(63);
v78Mchl=Enz_v(64);
vStarch2=Enz_v(65);
vgly1B=Enz_v(66);
vhexp=Enz_v(67);
vSta1=Enz_v(68);
vSta2=Enz_v(69);
vSta3=Enz_v(70);
vtATP=Enz_v(71);

vPCK1=Enz_v(72);
vPCK2=Enz_v(73);
vPCK3=Enz_v(74);
vPCK4=Enz_v(75);
vPCK5=Enz_v(76);
vPCK6=Enz_v(77);
vAsp=Enz_v(78);
vAla=Enz_v(79);
vPEP=Enz_v(80);
vOAA_B=Enz_v(81);
vATP_B=Enz_v(82);

vO2_Mchl=Enz_v(83);
vO2_Bchl=Enz_v(84);
vtO2=Enz_v(85);
vtO2_B=Enz_v(86);
vtO2_M=Enz_v(87);

VolMC=0.01;
VolMchl=0.02;
VolBSC=0.0045;
VolBchl=0.009;
Volper=0.00045;

global pathway_option;

if pathway_option==0;
    vPCK1=0;
    vPCK2=0;
    vPCK3=0;
    vPCK4=0;
    vPCK5=0;
    vPCK6=0;
    vAsp=0;
    vAla=0;
    vPEP=0;
    vOAA_B=0;
    vATP_B=0;
end

if pathway_option==1;
   vPCK3=0;
   vATP_B=0;
   vPEP=0;
end

if pathway_option==2;
   vPCK6=0; 
   vOAA_B=0;
end

if pathway_option==4;
   v3=0;
   v4=0;
   v5=0;
   vPCK6=0;
   vPYR_B=0;
   vPYR_M=0;
   vOAA_B=0;
   vOAA_M=0;
   vMAL_M=0;
   vMAL_B=0;
   vPEP_M=0;
end

if pathway_option==6;
    vMAL_B=0;
end


%global pathway_option;
global option2;
option2=1; % 0 or 1; if 0, assume O2 concentration is constant during simulation.
   
Delta_MC_CO2=(vinf-v1+vleak+vrpd)/VolMC;
Delta_MC_HCO3= (v1 - v2)/VolMC;
Delta_MC_Malate= (vMAL_M - vMAL)/VolMC;
Delta_MC_PGA= (vPGA - vPGA_M -v62- vPGASink)/VolMC;
Delta_MC_FBP= (vSuc1- vSuc2)/VolMC;
Delta_MC_UDPG= (vSuc7 - vSuc8)/VolMC;
Delta_MC_SUCP= (vSuc8 - vSuc9)/VolMC;
Delta_MC_SUC= (vSuc9 - vSuc10)/VolMC;
Delta_MC_F26BP= (vSuc3 - vSuc4)/VolMC;
Delta_MC_ATP= (vtATP-vSuc7-vSuc3)/VolMC;
Delta_MC_T3P= (vGAP_M + vDHAP_M - vGAP - vDHAP - 2*vSuc1)/VolMC;
Delta_MC_HexP= (vSuc2 + vSuc4 - vSuc3 - vSuc7 - vSuc8)/VolMC;
Delta_MC_Sucrose= (vSuc10)/VolMC;

Delta_Mchl_OAA= (vOAA_M - v3)/VolMchl;
Delta_Mchl_Malate= (v3 - vMAL_M)/VolMchl;
Delta_Mchl_PEP= (v5 - vPEP_M)/VolMchl;
Delta_Mchl_Pyruvate= (vPYR_M - v5)/VolMchl;
Delta_Mchl_NADPH= (vNADPHM - v3 - v78Mchl)/VolMchl;
Delta_Mchl_ATP= (vATPM - 2*v5 - v78Mchl-vtATP)/VolMchl;
Delta_Mchl_PGA= (vPGA_M - v78Mchl)/VolMchl;
Delta_Mchl_DPGA= 0;
Delta_Mchl_T3P= (v78Mchl - vGAP_M - vDHAP_M)/VolMchl;

Delta_BSC_T3P= (vGAP + vDHAP - vGAP_B - vDHAP_B)/VolBSC;
Delta_BSC_PGA= (vPGA_B - vPGA)/VolBSC;
Delta_BSC_Malate= (vMAL - vMAL_B)/VolBSC;

Delta_Bchl_CO2= (v4 - v6 - vleak_B)/VolBchl;
Delta_Bchl_RuBP= (v18 - v6-vpr1)/VolBchl;
Delta_Bchl_PGA= (2*v6 - v78 - vPGA_B + vpr1 + vpr8)/VolBchl;
Delta_Bchl_DPGA=0;
Delta_Bchl_SBP= (v12 - v13)/VolBchl;
Delta_Bchl_S7P= (v13 - v15)/VolBchl;
Delta_Bchl_FBP= (v10 - v11+vgly1B)/VolBchl;
Delta_Bchl_E4P= (v14 - v12)/VolBchl;
Delta_Bchl_Starch= vSta3/VolBchl;
Delta_Bchl_Rubisco= 0;
Delta_Bchl_T3P= (vGAP_B + vDHAP_B + v78 - 2*v10 - v14 - v15 - v12)/VolBchl;
Delta_Bchl_HexP= (v11 - v14 - vSta1 -vgly1B+vhexp)/VolBchl;
Delta_Bchl_Pent= (v14 + 2*v15 - v18)/VolBchl;
Delta_Bchl_Pyruvate= (v4 - vPYR_B)/VolBchl;
Delta_Bchl_PPi= (vSta1-vSta2)/VolBchl;
Delta_Bchl_ADPG= (vSta1-vSta3)/VolBchl;
Delta_Bchl_PGCA= (vpr1 - vpr2)/VolBchl;
Delta_Bchl_GCA= (vpr2 - vpr9)/VolBchl;
Delta_Bchl_GCEA= (vpr10 - vpr8)/VolBchl;

Delta_Bper_GCA= (vpr9 - vpr3)/Volper;
Delta_Bper_GOA= (vpr3 - vpr4 -vpr6)/Volper;
Delta_Bper_GLY= (vpr4 + vpr6 - 2*vpr5)/Volper;
Delta_Bper_SER= (vpr5 - vpr6)/Volper;
Delta_Bper_HPR= (vpr6 - vpr7)/Volper;
Delta_Bper_GCEA= (vpr7 - vpr10)/Volper;

Delta_MC_O2=0;
Delta_Mchl_O2=vO2_Mchl-vtO2_M;
Delta_BSC_O2=vtO2_B-vtO2;
Delta_Bchl_O2=vO2_Bchl-vtO2_B;

Delta_MC_OAA= (v2 - vOAA_M-vPCK1)/VolMC;
Delta_MC_Pyruvate= (vPYR - vPYR_M+vPCK5)/VolMC;
Delta_MC_PEP= (vPEP_M - v2 + v62 + vPEP )/VolMC;
Delta_BSC_Pyruvate= (vPYR_B - vPYR-vPCK4)/VolBSC;
Delta_BSC_OAA=(vPCK2-vOAA_B-vPCK3)/VolBSC;
Delta_BSC_PEP=(vPCK3-vPEP)/VolBSC;
Delta_BSC_ATP=vATP_B-vPCK3;
Delta_Bchl_NADPH=(v4 - v78-vPCK6+vNADPHB)/VolBchl;%+ vNADPHB;
Delta_Bchl_Malate= (vMAL_B - v4+vPCK6)/VolBchl; 
Delta_Bchl_OAA=(vOAA_B-vPCK6)/VolBchl;
Delta_Bchl_ATP=(vATPB - v78 - v18 - vSta1 - vpr8-vgly1B-vATP_B)/VolBchl;
   
Delta_MC_Glu=(vPCK5-vPCK1)/VolMC;
Delta_MC_OxoG=(vPCK1-vPCK5)/VolMC;
Delta_MC_Asp=(vPCK1-vAsp)/VolMC;
Delta_MC_Ala=(vAla-vPCK5)/VolMC;
Delta_BSC_OxoG=(vPCK4-vPCK2)/VolBSC;
Delta_BSC_Glu=(vPCK2-vPCK4)/VolBSC;
Delta_BSC_Asp=(vAsp-vPCK2)/VolBSC;
Delta_BSC_Ala=(vPCK4-vAla)/VolBSC;
   
Delta_BSC_CO2= (vpr5+vleak_B - vleak+vrpd+vPCK3)/VolBSC;
   

if option2 == 0
    Delta_MC_O2=0;
    Delta_Mchl_O2=0;
    Delta_BSC_O2=0;
    Delta_Bchl_O2=0;
end


EnzMB=zeros(58,1);
EnzMB(1)=Delta_MC_HCO3;
EnzMB(2)=Delta_MC_OAA;
EnzMB(3)=Delta_MC_PEP;
EnzMB(4)=Delta_MC_Malate;
EnzMB(5)=Delta_MC_Pyruvate;
EnzMB(6)=Delta_MC_PGA;
EnzMB(7)=Delta_MC_FBP;
EnzMB(8)=Delta_MC_UDPG;
EnzMB(9)=Delta_MC_SUCP;
EnzMB(10)=Delta_MC_SUC;
EnzMB(11)=Delta_MC_F26BP;
EnzMB(12)=Delta_MC_ATP;
EnzMB(13)=Delta_MC_T3P;
EnzMB(14)=Delta_MC_HexP;
EnzMB(15)=Delta_MC_Sucrose;
EnzMB(16)=Delta_Mchl_OAA;
EnzMB(17)=Delta_Mchl_Malate;
EnzMB(18)=Delta_Mchl_PEP;
EnzMB(19)=Delta_Mchl_Pyruvate;
EnzMB(20)=Delta_Mchl_NADPH;
EnzMB(21)=Delta_Mchl_ATP;
EnzMB(22)=Delta_Mchl_PGA;
EnzMB(23)=Delta_Mchl_DPGA;
EnzMB(24)=Delta_Mchl_T3P;
EnzMB(25)=Delta_BSC_T3P;
EnzMB(26)=Delta_BSC_PGA;
EnzMB(27)=Delta_BSC_Malate;
EnzMB(28)=Delta_BSC_Pyruvate;
EnzMB(29)=Delta_BSC_CO2;
EnzMB(30)=Delta_Bchl_CO2;
EnzMB(31)=Delta_Bchl_RuBP;
EnzMB(32)=Delta_Bchl_PGA;
EnzMB(33)=Delta_Bchl_DPGA;
EnzMB(34)=Delta_Bchl_ATP;
EnzMB(35)=Delta_Bchl_NADPH; 
EnzMB(36)=Delta_Bchl_SBP;
EnzMB(37)=Delta_Bchl_S7P;
EnzMB(38)=Delta_Bchl_FBP;
EnzMB(39)=Delta_Bchl_E4P;
EnzMB(40)=Delta_Bchl_Starch;
EnzMB(41)=Delta_Bchl_Rubisco;
EnzMB(42)=Delta_Bchl_T3P;
EnzMB(43)=Delta_Bchl_HexP;
EnzMB(44)=Delta_Bchl_Pent;
EnzMB(45)=Delta_Bchl_Malate;
EnzMB(46)=Delta_Bchl_Pyruvate;
EnzMB(47)=Delta_Bchl_PGCA;
EnzMB(48)=Delta_Bchl_GCA;
EnzMB(49)=Delta_Bchl_GCEA;
EnzMB(50)=Delta_Bper_GCA;
EnzMB(51)=Delta_Bper_GOA;
EnzMB(52)=Delta_Bper_GLY;
EnzMB(53)=Delta_Bper_SER;
EnzMB(54)=Delta_Bper_HPR;
EnzMB(55)=Delta_Bper_GCEA;
EnzMB(56)=Delta_MC_CO2;
EnzMB(57)=Delta_Bchl_PPi;
EnzMB(58)=Delta_Bchl_ADPG;

EnzMB(59)=Delta_MC_Glu;
EnzMB(60)=Delta_MC_OxoG;
EnzMB(61)=Delta_MC_Asp;
EnzMB(62)=Delta_MC_Ala;

EnzMB(63)=Delta_BSC_OxoG;
EnzMB(64)=Delta_BSC_Glu;
EnzMB(65)=Delta_BSC_Asp;
EnzMB(66)=Delta_BSC_Ala;
EnzMB(67)=Delta_BSC_OAA;
EnzMB(68)=Delta_BSC_PEP;
EnzMB(69)=Delta_BSC_ATP;
EnzMB(70)=Delta_Bchl_OAA;
EnzMB(71)=Delta_MC_O2;
EnzMB(72)=Delta_Mchl_O2;
EnzMB(73)=Delta_BSC_O2;
EnzMB(74)=Delta_Bchl_O2;
