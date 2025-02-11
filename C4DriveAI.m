% Copyright  Wang Y and Zhu X-G, 2014. 
% CAS-MPG Partner Institute for Computational Biology, Shanghai Institutes for Biological Sciences, CAS, Shanghai,200031 
% A mixed pathway model of C4 photosynthesis
% light curve simulation


clear all;

Input=importdata('C4model_input.txt');
InputData=Input.data;
Ca_Input=InputData(:,ismember(Input.colheaders,'CO2')==1);
Light_Input=InputData(:,ismember(Input.colheaders,'Light')==1);
Pathway_Input=InputData(:,ismember(Input.colheaders,'Pathway_type')==1);


for i=1:50
i
global phi;
global Lpd;
global I;
global U;
global V;
global Ratio;
global CI
global O2;
global Bchl_CP;
global MC_CP;
global Mchl_CP;
I = 0.05*i;%mmol m-2 s-1 light intensity input
Ii(i,1)=I;
phi=0.03; %¦Õ=0.03 plasmodesmata proportion
Lpd=400;  % plasmodesmata length   lPD=0.4¦Ìm
U=0;% light partition coefficient
V=0;% Jmax partition coefficient
Ratio=4; % Enezyme activity variation factor for PCK pathway

CI=Ca_Input*0.45/(3 * 10^4);%0.006%intercellular CO2 concentration
O2= 0.2646;% O2 concentration
Bchl_CP= 25.0;%Total phosphate concentration in bundle sheath chloroplast
MC_CP=15.0;%Total phosphate concentration in mesophyll cell cytosol
Mchl_CP=15.0;%Total phosphate concentration in mesophyll cell chloroplast

global pathway_option;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pathway_option=Pathway_Input;
%%% 0 normal NADP-ME type 
%%% 1 Asp+Mal transport and ME type 
%%% 2 Asp+Mal and PCK type 
%%% 3 Asp+Mal and PCK+ME type 
%%% 4 Asp and PCK only type
%%% 6 DiT2 mutant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Ini=C4Ini;
time=80000; %simulation time
[Tt, d]=ode15s(@C4MB, [0, time], Ini);
global Result;
Result =[Tt,d]; 

% Calculate the CO2 uptake rate
global Velocity_s;
n=size(Tt);
global v6;
global v2;
global vpr1;
global vc;
global vpr;
global vleakage;
global vMAL;
global vPGA;
global vpr1;
global vAsp;

Vm_6= Velocity_s(6);
vpr(i,1)=vpr1;
vc(i,1)= v6-0.5*vpr1;
leakiness(i,1)=vleakage/v2;
vPGAi(i,1)=vPGA;
vPGAratioi(i,1)=vPGA/v6;
MalRatio(i,1)=vMAL/v2;
AspRatio(i,1)=vAsp/v2;

Meta_con(i,1:70)=d(n(1),1:70);
MalGradient(i,1)=Meta_con(i,4)-Meta_con(i,27);
AspGradient(i,1)=Meta_con(i,61)-Meta_con(i,65);
PyrGradient(i,1)=Meta_con(i,28)-Meta_con(i,5);
AlaGradient(i,1)=Meta_con(i,66)-Meta_con(i,62);
PEPGradient(i,1)=Meta_con(i,68)-Meta_con(i,3);
BSCCO2(i,1)=Meta_con(i,29);
BSCChlCO2(i,1)=Meta_con(i,30);

end
Rd=1; % Respiration rate
A=1*1000*vc-Rd;% net co2 uptake
figure;
plot (Ii*1000,A,'k');
xlabel('PPFD (\mumolm^-^2 s^-^1)');
ylabel('A (\mumolm^-^2 s^-^1)');

DataI=zeros(i,15);
DataI(:,1)=Ii/20*1000;
DataI(:,2)=vPGAi*1000;
DataI(:,3)=vPGAratioi;
DataI(:,4)=A;
DataI(:,5)=leakiness;%leakiness
DataI(:,6)=MalGradient;
DataI(:,7)=AspGradient;
DataI(:,8)=PyrGradient;
DataI(:,9)=PEPGradient;
DataI(:,10)=AlaGradient;
DataI(:,11)=MalRatio;
DataI(:,12)=AspRatio;
DataI(:,13)=vpr*1000;%photorespiration
DataI(:,14)=BSCCO2;
DataI(:,15)=BSCChlCO2;

% Calculate the acid concentration in leaf
VolMC=0.01;
VolMchl=0.02;
VolBSC=0.0045;
VolBchl=0.009;
Volper=0.00045;
Volt=0.25;
Malate=(Meta_con(:,4)*VolMC+Meta_con(:,17)*VolMchl+Meta_con(:,27)*VolBSC+Meta_con(:,45)*VolBchl)*1000;
Pyruvate=(Meta_con(:,5)*VolMC+Meta_con(:,19)*VolMchl+Meta_con(:,28)*VolBSC+Meta_con(:,46)*VolBchl)*1000;
PEP=(Meta_con(:,3)*VolMC+Meta_con(:,18)*VolMchl+Meta_con(:,68)*VolBSC)*1000; 
Alanine=(Meta_con(:,62)*VolMC+Meta_con(:,66)*VolBSC)*1000;
Aspartate=(Meta_con(:,61)*VolMC+Meta_con(:,65)*VolBSC)*1000;
Datamata=zeros(i,5);
Datamata(:,1)=Malate;
Datamata(:,2)=Aspartate;
Datamata(:,3)=Pyruvate;
Datamata(:,4)=PEP;
Datamata(:,5)=Alanine;

% figure;%each metabolite concentrations
%    subplot(5,10,1); plot(Ii,Meta_con(:,1));
%    title('MC.HCO3');
%    subplot(5,10,2); plot(Ii,Meta_con(:,2));
%    title('MC.OAA');
%    subplot(5,10,3); plot(Ii,Meta_con(:,3));
%    title('MC.PEP');
%    subplot(5,10,4); plot(Ii,Meta_con(:,4));
%    title('MC.Malate');
%    subplot(5,10,5); plot(Ii,Meta_con(:,5));
%    title('MC.Pyruvate');
%    subplot(5,10,6); plot(Ii,Meta_con(:,6));
%    title('MC.PGA');
%    subplot(5,10,7); plot(Ii,Meta_con(:,7));
%    title('MC.FBP');
%    subplot(5,10,8); plot(Ii,Meta_con(:,8));
%    title('MC.UDPG');
%    subplot(5,10,9); plot(Ii,Meta_con(:,9));
%    title('MC.SUCP');
%    subplot(5,10,10); plot(Ii,Meta_con(:,10));
%    title('MC.SUC');
%    subplot(5,10,11); plot(Ii,Meta_con(:,11));
%    title('MC.F26BP');
%    subplot(5,10,12); plot(Ii,Meta_con(:,12));
%    title('MC.ATP');
%    subplot(5,10,13); plot(Ii,Meta_con(:,13));
%    title('MC.T3P');
%    subplot(5,10,14); plot(Ii,Meta_con(:,14));
%    title('MC.HexP');
%    subplot(5,10,15); plot(Ii,Meta_con(:,15));
%    title('MC.Sucrose');
%    subplot(5,10,16); plot(Ii,Meta_con(:,16));
%    title('Mchl.OAA');
%    subplot(5,10,17); plot(Ii,Meta_con(:,17));
%    title('Mchl.Malate');
%    subplot(5,10,18); plot(Ii,Meta_con(:,18));
%    title('Mchl.PEP');
%    subplot(5,10,19); plot(Ii,Meta_con(:,19));
%    title('Mchl.Pyruvate');
%    subplot(5,10,20); plot(Ii,Meta_con(:,20));
%    title('Mchl.NADPH');
%    subplot(5,10,21); plot(Ii,Meta_con(:,21));
%    title('Mchl.ATP');
%    subplot(5,10,22); plot(Ii,Meta_con(:,22));
%    title('Mchl.PGA');
%    subplot(5,10,23); plot(Ii,Meta_con(:,23));
%    title('Mchl.DPGA');
%    subplot(5,10,24); plot(Ii,Meta_con(:,24));
%    title('Mchl.T3P');
%    subplot(5,10,25); plot(Ii,Meta_con(:,25));
%    title('BSC.T3P');
%    subplot(5,10,26); plot(Ii,Meta_con(:,26));
%    title('BSC.PGA');
%    subplot(5,10,27); plot(Ii,Meta_con(:,27));
%    title('BSC.Malate');
%    subplot(5,10,28); plot(Ii,Meta_con(:,28));
%    title('BSC.Pyruvate');
%    subplot(5,10,29); plot(Ii,Meta_con(:,29));
%    title('BSC.CO2');
%    subplot(5,10,30); plot(Ii,Meta_con(:,30));
%    title('Bchl.CO2');
%    subplot(5,10,31); plot(Ii,Meta_con(:,31));
%    title('Bchl.RuBP');
%    subplot(5,10,32); plot(Ii,Meta_con(:,32));
%    title('Bchl.PGA');
%    subplot(5,10,33); plot(Ii,Meta_con(:,33));
%    title('Bchl.DPGA');
%    subplot(5,10,34); plot(Ii,Meta_con(:,34));
%    title('Bchl.ATP');
%    subplot(5,10,35); plot(Ii,Meta_con(:,35));
%    title('Bchl.NADPH');
%    subplot(5,10,36); plot(Ii,Meta_con(:,36));
%    title('Bchl.SBP');
%    subplot(5,10,37); plot(Ii,Meta_con(:,37));
%    title('Bchl.S7P');
%    subplot(5,10,38); plot(Ii,Meta_con(:,38));
%    title('Bchl.FBP');
%    subplot(5,10,39); plot(Ii,Meta_con(:,39));
%    title('Bchl.E4P');
%    subplot(5,10,40); plot(Ii,Meta_con(:,40));
%    title('Bchl.Starch');
%    subplot(5,10,41); plot(Ii,Meta_con(:,41));
%    title('Bchl.Rubisco');
%    subplot(5,10,42); plot(Ii,Meta_con(:,42));
%    title('Bchl.T3P');
%    subplot(5,10,43); plot(Ii,Meta_con(:,43));
%    title('Bchl.HexP');
%    subplot(5,10,44); plot(Ii,Meta_con(:,44));
%    title('Bchl.Pent');
%    subplot(5,10,45); plot(Ii,Meta_con(:,45));
%    title('Bchl.Malate');
%    subplot(5,10,46); plot(Ii,Meta_con(:,46));
%    title('Bchl.Pyruvate'); 
   

   
   
   