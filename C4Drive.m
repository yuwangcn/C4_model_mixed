
% Copyright  Wang Y and Zhu X-G, 2014. 
% CAS-MPG Partner Institute for Computational Biology, Shanghai Institutes for Biological Sciences, CAS, Shanghai,200031 
% A mixed pathway model of C4 photosynthesis
% The CO2 uptake and metabolite concentraion changes with time.


clear all;

Input=importdata('C4model_input.txt');
InputData=Input.data;
Ca_Input=InputData(:,ismember(Input.colheaders,'CO2')==1);
Light_Input=InputData(:,ismember(Input.colheaders,'Light')==1);
Pathway_Input=InputData(:,ismember(Input.colheaders,'Pathway_type')==1);

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

I = Light_Input/1000;%mmol m-2 s-1 light intensity input
CI=Ca_Input*0.45/(3 * 10^4);%0.006;%intercellular CO2 concentration

phi=0.03; %¦Õ=0.03 plasmodesmata proportion
Lpd=400;  % plasmodesmata length   lPD=0.4¦Ìm
U=0;% light partition coefficient
V=0;% Jmax partition coefficient
Ratio=4; % % Enezyme activity variation factor for PCK pathway


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
time=3000;%simulation time
[Tt, d]=ode15s(@C4MB, [0, time], Ini);
global Result;
Result =[Tt,d]; 

gm=0.7;%molm-2 bar-10.2
Sc=3*10^4;%3.36*10^4;%ubarL/mmol
vinf=gm*Sc*10^(-3)*(CI-d(:,56));%vinf=gm*Sc*10^(-3)*(CI-MC_CO2);
figure;
plot(Tt,vinf*1000) %CO2 influx of the model
xlabel('TIme (s)');
ylabel('CO_2 uptake rate (\mumolm^-^2 s^-^1)');
ylim([0,60]);
figure
subplot(5,10,1); plot(Tt,d(:,1));
title('MC.HCO3');
subplot(5,10,2); plot(Tt,d(:,2));
title('MC.OAA');
subplot(5,10,3); plot(Tt,d(:,3));
 title('MC.PEP');
 subplot(5,10,4); plot(Tt,d(:,4));
 title('MC.Malate');
 subplot(5,10,5); plot(Tt,d(:,5));
 title('MC.Pyruvate');
 subplot(5,10,6); plot(Tt,d(:,6));
 title('MC.PGA');
 subplot(5,10,7); plot(Tt,d(:,7));
 title('MC.FBP');
 subplot(5,10,8); plot(Tt,d(:,8));
 title('MC.UDPG');
 subplot(5,10,9); plot(Tt,d(:,9));
 title('MC.SUCP');
 subplot(5,10,10); plot(Tt,d(:,10));
 title('MC.SUC');
 subplot(5,10,11); plot(Tt,d(:,11));
 title('MC.F26BP');
 subplot(5,10,12); plot(Tt,d(:,12));
 title('MC.ATP');
 subplot(5,10,13); plot(Tt,d(:,13));
 title('MC.T3P');
 subplot(5,10,14); plot(Tt,d(:,14));
 title('MC.HexP');
 subplot(5,10,15); plot(Tt,d(:,15));
 title('MC.Sucrose');
 subplot(5,10,16); plot(Tt,d(:,16));
 title('Mchl.OAA');
 subplot(5,10,17); plot(Tt,d(:,17));
 title('Mchl.Malate');
 subplot(5,10,18); plot(Tt,d(:,18));
 title('Mchl.PEP');
 subplot(5,10,19); plot(Tt,d(:,19));
 title('Mchl.Pyruvate');
 subplot(5,10,20); plot(Tt,d(:,20));
 title('Mchl.NADPH');
 subplot(5,10,21); plot(Tt,d(:,21));
 title('Mchl.ATP');
 subplot(5,10,22); plot(Tt,d(:,22));
 title('Mchl.PGA');
 subplot(5,10,23); plot(Tt,d(:,24));
 title('Mchl.T3P');
 subplot(5,10,24); plot(Tt,d(:,25));
 title('BSC.T3P');
 subplot(5,10,25); plot(Tt,d(:,26));
 title('BSC.PGA');
 subplot(5,10,26); plot(Tt,d(:,27));
 title('BSC.Malate');
 subplot(5,10,27); plot(Tt,d(:,28));
 title('BSC.Pyruvate');
 subplot(5,10,28); plot(Tt,d(:,29));
 title('BSC.CO2');
 subplot(5,10,29); plot(Tt,d(:,30));
 title('Bchl.CO2');
 subplot(5,10,30); plot(Tt,d(:,31));
 title('Bchl.RuBP');
 subplot(5,10,31); plot(Tt,d(:,32));
 title('Bchl.PGA');
 subplot(5,10,32); plot(Tt,d(:,34));
 title('Bchl.ATP');
 subplot(5,10,33); plot(Tt,d(:,35));
 title('Bchl.NADPH');
 subplot(5,10,34); plot(Tt,d(:,36));
 title('Bchl.SBP');
 subplot(5,10,35); plot(Tt,d(:,37));
 title('Bchl.S7P');
 subplot(5,10,36); plot(Tt,d(:,38));
 title('Bchl.FBP');
 subplot(5,10,37); plot(Tt,d(:,39));
 title('Bchl.E4P');
 subplot(5,10,38); plot(Tt,d(:,40));
 title('Bchl.Starch');
 subplot(5,10,39); plot(Tt,d(:,41));
 title('Bchl.Rubisco');
 subplot(5,10,40); plot(Tt,d(:,42));
 title('Bchl.T3P');
 subplot(5,10,41); plot(Tt,d(:,43));
 title('Bchl.HexP');
 subplot(5,10,42); plot(Tt,d(:,44));
 title('Bchl.Pent');
 subplot(5,10,43); plot(Tt,d(:,45));
 title('Bchl.Malate');
 subplot(5,10,44); plot(Tt,d(:,46));
 title('Bchl.Pyruvate'); 

