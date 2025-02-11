Model description
---
The C4 metabolic model is a MATLAB program developed for simulating metabolic fluxes and metabolite concentrations for C4 photosynthesis. 

Contact Author
---
Yu Wang - https://github.com/yuwangcn

Runing steps:
---
**1.	Copy all the files to a directory, then open MATLAB and set the working directory to the directory where the code files are saved.**

**2.	Parameter setting：**

Before simulation, the following three parameters need to be defined.

Table 1 Input parameters of the C4 model

Input parameter| Description| unit
---- | ----- | ----- |
CO2(Ca)| Air CO2 concentration| µmol mol-1
Light| photosynthetic photon flux density| µmol m-2 s-1
C4 subtype| NADP-ME or different mixtures| NA

C4 subtype：

0: normal NADP-ME subtype 

1: Asp+Mal transport and ME subtype 

2: Asp+Mal and PCK subtype 

3: Asp+Mal and PCK+ME subtype 

4: Asp and PCK only subtype


**3.Output：**

Output of this model can be dynamic metabolite concentrations and reaction rate changes or steady state photosynthesis rate at different environmental conditions. 

Run "C4Drive.m" to simulate the photosynthesis and metabolite concentration changes with time; 

Run "C4DriveAI.m" to simulate the steady-state light response curve.

