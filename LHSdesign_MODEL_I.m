% This code generates the sampling data for the
% Manuscript: Ex vivo experiments shed light on the innate immune response from influenza virus
% Authors: Olmos, Nunes & Saenz
% Journal: Bulletin of Mathematical Biology (BMAB)

% This code corresponds to Model I. Basic Viral Dynamics

clear Mparam
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

Nsamples = 20000; 
Npar = 2;%number of parameters

%%%%%
% Varied parameter values (baseline values to make sampling regions)
%%%%%
% Infection rate, ie. ~probability (depending on free virus in the patch)
% that a healthy cell gets infected (per unit of time)
% Secretion rate; amount of released virus from infected-secreting cells (per unit of time)
%%%%%

% latin hipercube sampling
% in a unit (0,1)^Npar hipercube
Mparam = lhsdesign(Nsamples,Npar);
%disp(Mparam)

%infection rate
Mparam(:,1)=(10.^(10*Mparam(:,1)-5))*0.266; % 0.266/tick=8/hour (8/h / 30ticks) --Beauchemin et al 2006
%virus secretion rate
Mparam(:,2)=(10.^(10*Mparam(:,2)-5))*0.00166;% 0.00166 /tick=0.05 /hour (0.05/h / 30ticks) --Beauchemin et al 2006


%create and open text files to save results
fidparam=fopen('LHS_MODEL_I.txt','w');
fprintf(fidparam,'%10.12f\t %10.12f\t \n',Mparam');
fclose('all');


