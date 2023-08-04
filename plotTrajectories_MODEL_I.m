% This code is used to plot the simulations obtained from
% Simulations_ComputeErrorandTrajectories_MODEL_I.m
%
% Manuscript: Ex vivo experiments shed light on the innate immune response from influenza virus
% Authors: Olmos, Nunes & Saenz
% Journal: Bulletin of Mathematical Biology (BMAB)
%
% This code corresponds to Model I. Basic Viral Dynamics

load trajectories_MODEL_I.mat trajectories


datares=trajectories;

% %data results in a dataresults.mat file
%
% trajectories.inf_rate = INF_RATE;
% trajectories.sec_rate = SEC_RATE;
% trajectories.clearance_rate = CLEARANCE_RATE;
% trajectories.errorCells = ERROR_CELLS;
% trajectories.errorVirus = ERROR_VIRUS;
% trajectories.timeatmax = TIMEATMAX;
% trajectories.virusfree=VirusTotalMat;
% trajectories.healthycells=HealthyCellsMat;
% trajectories.eclipsecells=EclipseCellsMat;
% trajectories.secretingcells=SecretingCellsMat;

errorsum=trajectories.errorCells/0.4198+trajectories.errorVirus/6.19;
%(1/(30*24))*trajectories.timeatmax
[trajectories.inf_rate,trajectories.sec_rate,errorsum]

% Experimental data
time_data=1:4;
infectedcells_data=[8.35, 18.58, 41.98, 11.20]/100; % proportion of infected cells in experimental data
freevirus_data=[4.65, 5.87, 6.19, 5.97];%(log10(PT40)+log10(PT41))/2 in exp data

Ntot=length(trajectories.inf_rate)

%
% Plot trajectories for individual combinations
%
time=(1/(24*30))*(0:2879);%time in days (2880 points)
day4=2880;% in #ticks

for indiv=1:Ntot

    %cells alive
    liveCells=datares.healthycells(indiv,:)+datares.eclipsecells(indiv,:)+...
              datares.secretingcells(indiv,:);
    %proportion of infected cells relative to cells alive
    infectedcells=(datares.eclipsecells(indiv,:)+datares.secretingcells(indiv,:))./liveCells;

    %figure(2)
    %healthy cells
    subplot(2,2,1)
    plot(time(1:day4),datares.healthycells(indiv,1:day4),'b')
    hold on
    plot(time(1:day4),liveCells(1:day4),'k')
    ylabel("Cells")
    xlabel("Days")
    xlim([0,4.1])
    ylim([0,4.1*10^5])
    legend('Susceptible','Alive','Location','SW')
    
    %secreating cells
    subplot(2,2,2)
    plot(time(1:day4),datares.eclipsecells(indiv,1:day4),'m')
    hold on
    plot(time(1:day4),datares.secretingcells(indiv,1:day4),'r')
    ylabel("Cells")
    xlabel("Days")
    xlim([0,4.1])
    legend('Eclipse','Secreting','Location','NW')
    
    %total free virus
    subplot(2,2,3)
    plot(time,log10(datares.virusfree(indiv,:)),'b')
    hold on
    plot(time_data,freevirus_data,'o','MarkerEdgeColor','r','MarkerFaceColor','r')
    ylabel("Free virus")
    xlabel("Days")
    xlim([0,4.1])
    %infected cells proportion
    subplot(2,2,4)
    plot(time,infectedcells,'b')
    hold on
    plot(time_data,infectedcells_data,'o','MarkerEdgeColor','r','MarkerFaceColor','r')
    ylabel("Infected cells proportion")
    xlabel("Days")
    xlim([0,4.1])

end


