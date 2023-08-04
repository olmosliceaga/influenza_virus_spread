% Analysis of CA simulations for virus model
% run simulations with parameter values with smallest errors
% save results of trajectories for: 
% 1)free virus, and 
% 2)proportion of infected cells (with respect to alive cells)

% This code runs the automata code to show how a virus propagate on swine
% throat (ex-vivo) tissue.
% Under the basic viral dynamics model by Beauchemin et al (2006)
% Assuming an hexagonal grid

% This code corresponds to Model III. Immune response induces a refractory state

% Manuscript: Ex vivo experiments shed light on the innate immune response from influenza virus
% Authors: Olmos, Nunes & Saenz
% Journal: Bulletin of Mathematical Biology (BMAB)

close all
tinitprocess=cputime;
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

% The matrices A_{} are of size MxN. 
% A_{} contains the appropriate info of the stage of each cell.
% Although assuming an HEXAGONAL GRID, a squared arrangement (a matrix) is
% used to represent each cell (their status, viral load, and so on) by
% "dephasing" the even columns (see Figure 1 in the manuscript)
% M*N=633*633=400,689 is approximately the number of cells in the ex-vivo tissue
M=633;
N=633;

%%%%%
%
% Varied parameter values (baseline values to make sampling regions)
%
% Infection rate, ie. InfectionRate*FreeVirus(on the patch) = probability 
% that a healthy cell gets infected (per unit of time)

% Secretion rate; amount of released virus from infected-secreting cells (per unit of time)
%%%%%

%LHS_MODEL_III.txt is generated on LHSdesign_MODEL_III.m
param=load('LHS_MODEL_III.txt');

[pr,pc]=size(param);
%pr=number of replicates
%pc=number of parameters

%%%%%
%
% (Fixed) parameter values
%
% Number of ticks that the cell remains in eclipse stage
ticks_eclipse=7*30; % 7 hours (Beauchemin et al 2006)
% Number of ticks that the cell remains in viral secretion stage
ticks_secreting=23*30; % 23 hours (Beauchemin et al 2006)
%
%time before infected cell starts producing interferon
ticks_interferon_delay=8*30;
% Dv is the free virus proportion that moves from each cell to its neighbors per unit of time
Dv=0.0126;  % Dv= 0.0126 = 4*D*dt/(dx^2) = (4*3.18*10^(-15)*2*60)/((11*10^(-6))^2)
% beacause from Beauchemin et al (2006):
% D=3.18*10^(-15) m^2/s (diffussion coefficient in PDEs),
% dt=2 minutes, dx=11 micrometers
% 
Di=0.32;%Diffusion interferon value

% initial amount of free virus in each patch
init_amount=0.00125; %this number is chosen so that init_amount*N*M ~= 500 pfu (as in experiments)
% the above implies that units of free virus are PFUs
%
% number of iterations for each run (ie, duration of the simulated infection)
TotalTicks=2880; % each iteration (tick) is 2 minutes, so that simulated time is 4.16 days 
        % (need at least 4 days to compare to experimental data)
%
%%%%%%%%%%%%%%%%%%%%
%
% DATA TO FIT
%
% Experimental data
% Sw_Eng06.xlsx
% 500 pfu (initial amount of virus inoculated)
%
% Time (measurements at days 1,2,3,4 post-inoculation)
time_data=[24*30;48*30;72*30;96*30]; % in number of ticks (each tick=2 minutes)
% virus replication (PFU/ml)
virus_free_data_log_avg=[4.65;5.87;6.19;5.98]; % (log10(PT40)+log10(PT41))/2; only 2-significant decimals
%Proportion of infected cells
per_infect_cells=[8.35;18.58;41.98;11.20]/100; % segment A (proximal trachea) in dataset
%
%%%%%%%%%%%%%%%%%%%%

% number of replicates
TotalRuns=pr;%pr;

% create the space to save results
INF_RATE=zeros(TotalRuns,1);
SEC_RATE=zeros(TotalRuns,1);
IMMUNITY_RATE=zeros(TotalRuns,1);
% CLEARANCE_RATE=zeros(TotalRuns,1);

ERROR_CELLS=zeros(TotalRuns,1);
ERROR_VIRUS=zeros(TotalRuns,1);

for iRun=1:TotalRuns % Run Nruns examples and compute associated errors
    indice=1;
    
    A_healthy=ones(M,N); % Here all the cells are healthy. %1 healthy. 0 Not healthy
    A_eclipse=zeros(M,N);% No cells in eclipse stage
    A_secreting=zeros(M,N); % No cells in secreting stage
    A_death=zeros(M,N); % No cells in death stage
    A_immune=zeros(M,N);% No cells in the immune stage
 
    VirusTot=zeros(1,TotalTicks); % Keep record of total amount of free virus at each time step
    IFNTot=zeros(1,TotalTicks); % Keep record of total amount of interferon at each time step
    HealthyCells=zeros(1,TotalTicks); % Keep record of total amount of healthy cells at each time step
    EclipseCells=zeros(1,TotalTicks); % Keep record of total amount of eclipse cells at each time step
    SecretingCells=zeros(1,TotalTicks); % Keep record of total amount of secreting cells at each time step
  %  DeathCells=zeros(1,TotalTicks);
    ImmuneCells=zeros(1,TotalTicks);
    
%    A_death=zeros(M,N); % No cells in death stage  % UNUSED FOR ERROR COMPUTATION
%            (MAY BE USEFUL FOR PLOTTING TRAJECTORIES)
    
    vf=init_amount*ones(M,N); % This matrix saves the info of the amount of free virus in each cell.
    
    %Infection rate
    inf_rate=param(iRun,1); 
    %Secretion rate 
    sec_rate=param(iRun,2);     
    %Interferon rate
    interferon_sec_rate=1;  %FIXED
    %Immunity rate
    immunity_rate=param(iRun,3);
    %Viral clearance rate
    clearance_rate=0.0028;%per minute (c=4d^-1 => 1/c=6hours)
    clearance_proportion=clearance_rate*2; %clearance per tick (every 2 minutes)
 
    %I use these matrices to record the inner clock of a cell in a particular
    %stage
    
    status_interferon=zeros(M,N);
    interferon_delay_counter=zeros(M,N);
    free_interferon=zeros(M,N);
    Eclipse_counter=zeros(M,N);
    Secreting_counter=zeros(M,N);
      
    %Initialize errors to 0 at each new run
    errorCells=0;
    errorVirus=0;
    
     time_series=[0];
     healthy_cells=M*N;
     eclipse_cells=0;
     secreting_cells=0;
     dead_cells=0;
     virusfree=sum(sum(vf));
    
     for iTick=1:TotalTicks %Iterate the automata
         
         vf_ad=(1-Dv)*vf; %vf_ad virusfree left after diffusion
         fi_ad=(1-Di)*free_interferon;
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %
         % Diffusion of free virus
         %
         % Different patterns for odd and even columns
         % Consider boundary cells explicitely (changes in number of neighbors)
         % (see Figure 1 in the manuscript)
         %
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %
         % Odd columns
         %
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         j=1; % First column
         i=1; % (only 2 neighbors)
         VFtoNeig=(Dv/6)*vf(i,j);
         FItoNeig=(Di/6)*free_interferon(i,j);
         
         %vf_ad(i-1,j-1)=vf_ad(i-1,j-1)+VFtoNeig; % outside of grid
         %vf_ad(i-1,j)=vf_ad(i-1,j)+VFtoNeig;     % outside of grid
         %vf_ad(i-1,j+1)=vf_ad(i-1,j+1)+VFtoNeig; % outside of grid
         %vf_ad(i,j-1)=vf_ad(i,j-1)+VFtoNeig;     % outside of grid
         vf_ad(i+1,j)=vf_ad(i+1,j)+VFtoNeig;
         vf_ad(i,j+1)=vf_ad(i,j+1)+VFtoNeig;
         vf_ad(i,j)=vf_ad(i,j)+4*VFtoNeig; % reflective boundary conditions (dV/di=0)
         
         %fi_ad(i-1,j-1)=fi_ad(i-1,j-1)+FItoNeig; % outside of grid
         %fi_ad(i-1,j)=fi_ad(i-1,j)+FItoNeig;     % outside of grid
         %fi_ad(i-1,j+1)=fi_ad(i-1,j+1)+FItoNeig; % outside of grid
         %fi_ad(i,j-1)=fi_ad(i,j-1)+FItoNeig;     % outside of grid
         fi_ad(i+1,j)=fi_ad(i+1,j)+FItoNeig;
         fi_ad(i,j+1)=fi_ad(i,j+1)+FItoNeig;
         fi_ad(i,j)=fi_ad(i,j)+4*FItoNeig; % reflective boundary conditions (dV/di=0)
         
         
         for i=2:(M-1) % only 4 neighbors
             VFtoNeig=(Dv/6)*vf(i,j);
             FItoNeig=(Di/6)*free_interferon(i,j);
             
             %vf_ad(i-1,j-1)=vf_ad(i-1,j-1)+VFtoNeig; % outside of grid
             vf_ad(i-1,j)=vf_ad(i-1,j)+VFtoNeig;
             vf_ad(i-1,j+1)=vf_ad(i-1,j+1)+VFtoNeig;
             %vf_ad(i,j-1)=vf_ad(i,j-1)+VFtoNeig;     % outside of grid
             vf_ad(i+1,j)=vf_ad(i+1,j)+VFtoNeig;
             vf_ad(i,j+1)=vf_ad(i,j+1)+VFtoNeig;
             vf_ad(i,j)=vf_ad(i,j)+2*VFtoNeig; % reflective boundary conditions (dV/di=0)
             
             %fi_ad(i-1,j-1)=fi_ad(i-1,j-1)+FItoNeig; % outside of grid
             fi_ad(i-1,j)=fi_ad(i-1,j)+FItoNeig;
             fi_ad(i-1,j+1)=fi_ad(i-1,j+1)+FItoNeig;
             %fi_ad(i,j-1)=fi_ad(i,j-1)+FItoNeig;     % outside of grid
             fi_ad(i+1,j)=fi_ad(i+1,j)+FItoNeig;
             fi_ad(i,j+1)=fi_ad(i,j+1)+FItoNeig;
             fi_ad(i,j)=fi_ad(i,j)+2*FItoNeig; % reflective boundary conditions (dV/di=0)
         end
         i=M; % only 3 neighbors
         VFtoNeig=(Dv/6)*vf(i,j);
         FItoNeig=(Di/6)*free_interferon(i,j);
         
         %vf_ad(i-1,j-1)=vf_ad(i-1,j-1)+VFtoNeig; % outside of grid
         vf_ad(i-1,j)=vf_ad(i-1,j)+VFtoNeig;
         vf_ad(i-1,j+1)=vf_ad(i-1,j+1)+VFtoNeig;
         %vf_ad(i,j-1)=vf_ad(i,j-1)+VFtoNeig;    % outside of grid
         %vf_ad(i+1,j)=vf_ad(i+1,j)+VFtoNeig;    % outside of grid
         vf_ad(i,j+1)=vf_ad(i,j+1)+VFtoNeig;
         vf_ad(i,j)=vf_ad(i,j)+3*VFtoNeig; % reflective boundary conditions (dV/di=0)
         
         %fi_ad(i-1,j-1)=fi_ad(i-1,j-1)+FItoNeig; % outside of grid
         fi_ad(i-1,j)=fi_ad(i-1,j)+FItoNeig;
         fi_ad(i-1,j+1)=fi_ad(i-1,j+1)+FItoNeig;
         %fi_ad(i,j-1)=fi_ad(i,j-1)+FItoNeig;    % outside of grid
         %fi_ad(i+1,j)=fi_ad(i+1,j)+FItoNeig;    % outside of grid
         fi_ad(i,j+1)=fi_ad(i,j+1)+FItoNeig;
         fi_ad(i,j)=fi_ad(i,j)+3*FItoNeig; % reflective boundary conditions (dV/di=0)
         for j=3:2:(N-2) % inner columns
             i=1; % (only 3 neighbors)
             VFtoNeig=(Dv/6)*vf(i,j);
             FItoNeig=(Di/6)*free_interferon(i,j);
             
             %vf_ad(i-1,j-1)=vf_ad(i-1,j-1)+VFtoNeig; % outside of grid
             %vf_ad(i-1,j)=vf_ad(i-1,j)+VFtoNeig;     % outside of grid
             %vf_ad(i-1,j+1)=vf_ad(i-1,j+1)+VFtoNeig; % outside of grid
             vf_ad(i,j-1)=vf_ad(i,j-1)+VFtoNeig;
             vf_ad(i+1,j)=vf_ad(i+1,j)+VFtoNeig;
             vf_ad(i,j+1)=vf_ad(i,j+1)+VFtoNeig;
             vf_ad(i,j)=vf_ad(i,j)+3*VFtoNeig; % reflective boundary conditions (dV/di=0)
             
             %fi_ad(i-1,j-1)=fi_ad(i-1,j-1)+FItoNeig; % outside of grid
             %fi_ad(i-1,j)=fi_ad(i-1,j)+FItoNeig;     % outside of grid
             %fi_ad(i-1,j+1)=fi_ad(i-1,j+1)+FItoNeig; % outside of grid
             fi_ad(i,j-1)=fi_ad(i,j-1)+FItoNeig;
             fi_ad(i+1,j)=fi_ad(i+1,j)+FItoNeig;
             fi_ad(i,j+1)=fi_ad(i,j+1)+FItoNeig;
             fi_ad(i,j)=fi_ad(i,j)+3*FItoNeig; % reflective boundary conditions (dV/di=0)
             
             for i=2:(M-1) % Inner cells (not in the boundary)
                 VFtoNeig=(Dv/6)*vf(i,j);
                 FItoNeig=(Di/6)*free_interferon(i,j);
                 
                 vf_ad(i-1,j-1)=vf_ad(i-1,j-1)+VFtoNeig;
                 vf_ad(i-1,j)=vf_ad(i-1,j)+VFtoNeig;
                 vf_ad(i-1,j+1)=vf_ad(i-1,j+1)+VFtoNeig;
                 vf_ad(i,j-1)=vf_ad(i,j-1)+VFtoNeig;
                 vf_ad(i+1,j)=vf_ad(i+1,j)+VFtoNeig;
                 vf_ad(i,j+1)=vf_ad(i,j+1)+VFtoNeig;
                 
                 fi_ad(i-1,j-1)=fi_ad(i-1,j-1)+FItoNeig;
                 fi_ad(i-1,j)=fi_ad(i-1,j)+FItoNeig;
                 fi_ad(i-1,j+1)=fi_ad(i-1,j+1)+FItoNeig;
                 fi_ad(i,j-1)=fi_ad(i,j-1)+FItoNeig;
                 fi_ad(i+1,j)=fi_ad(i+1,j)+FItoNeig;
                 fi_ad(i,j+1)=fi_ad(i,j+1)+FItoNeig;
             end
             i=M; % only 5 neighbors
             VFtoNeig=(Dv/6)*vf(i,j);
             FItoNeig=(Di/6)*free_interferon(i,j);
                        
             vf_ad(i-1,j-1)=vf_ad(i-1,j-1)+VFtoNeig;
             vf_ad(i-1,j)=vf_ad(i-1,j)+VFtoNeig;
             vf_ad(i-1,j+1)=vf_ad(i-1,j+1)+VFtoNeig;
             vf_ad(i,j-1)=vf_ad(i,j-1)+VFtoNeig;
             %vf_ad(i+1,j)=vf_ad(i+1,j)+VFtoNeig; % outside of grid
             vf_ad(i,j+1)=vf_ad(i,j+1)+VFtoNeig;
             vf_ad(i,j)=vf_ad(i,j)+1*VFtoNeig; % reflective boundary conditions (dV/di=0)
             
             fi_ad(i-1,j-1)=fi_ad(i-1,j-1)+FItoNeig;
             fi_ad(i-1,j)=fi_ad(i-1,j)+FItoNeig;
             fi_ad(i-1,j+1)=fi_ad(i-1,j+1)+FItoNeig;
             fi_ad(i,j-1)=fi_ad(i,j-1)+FItoNeig;
             %fi_ad(i+1,j)=fi_ad(i+1,j)+FItoNeig; % outside of grid
             fi_ad(i,j+1)=fi_ad(i,j+1)+FItoNeig;
             fi_ad(i,j)=fi_ad(i,j)+1*FItoNeig; % reflective boundary conditions (dV/di=0)
           
         end
         j=N; % Last column
         i=1; % (only 2 neighbors)
         VFtoNeig=(Dv/6)*vf(i,j);
         FItoNeig=(Di/6)*free_interferon(i,j);
         
         %vf_ad(i-1,j-1)=vf_ad(i-1,j-1)+VFtoNeig; % outside of grid
         %vf_ad(i-1,j)=vf_ad(i-1,j)+VFtoNeig;     % outside of grid
         %vf_ad(i-1,j+1)=vf_ad(i-1,j+1)+VFtoNeig; % outside of grid
         vf_ad(i,j-1)=vf_ad(i,j-1)+VFtoNeig;
         vf_ad(i+1,j)=vf_ad(i+1,j)+VFtoNeig;
         %vf_ad(i,j+1)=vf_ad(i,j+1)+VFtoNeig;     % outside of grid
         vf_ad(i,j)=vf_ad(i,j)+4*VFtoNeig; % reflective boundary conditions (dV/di=0)
         
         %fi_ad(i-1,j-1)=fi_ad(i-1,j-1)+FItoNeig; % outside of grid
         %fi_ad(i-1,j)=fi_ad(i-1,j)+FItoNeig;     % outside of grid
         %fi_ad(i-1,j+1)=fi_ad(i-1,j+1)+FItoNeig; % outside of grid
         fi_ad(i,j-1)=fi_ad(i,j-1)+FItoNeig;
         fi_ad(i+1,j)=fi_ad(i+1,j)+FItoNeig;
         %fi_ad(i,j+1)=fi_ad(i,j+1)+FItoNeig;     % outside of grid
         fi_ad(i,j)=fi_ad(i,j)+4*FItoNeig; % reflective boundary conditions (dV/di=0)
         
         for i=2:(M-1) % only 4 neighbors
             VFtoNeig=(Dv/6)*vf(i,j);
             FItoNeig=(Di/6)*free_interferon(i,j);
             
             vf_ad(i-1,j-1)=vf_ad(i-1,j-1)+VFtoNeig;
             vf_ad(i-1,j)=vf_ad(i-1,j)+VFtoNeig;
             %vf_ad(i-1,j+1)=vf_ad(i-1,j+1)+VFtoNeig; % outside of grid
             vf_ad(i,j-1)=vf_ad(i,j-1)+VFtoNeig;
             vf_ad(i+1,j)=vf_ad(i+1,j)+VFtoNeig;
             %vf_ad(i,j+1)=vf_ad(i,j+1)+VFtoNeig;     % outside of grid
             vf_ad(i,j)=vf_ad(i,j)+2*VFtoNeig; % reflective boundary conditions (dV/di=0)
             
             fi_ad(i-1,j-1)=fi_ad(i-1,j-1)+FItoNeig;
             fi_ad(i-1,j)=fi_ad(i-1,j)+FItoNeig;
             %fi_ad(i-1,j+1)=fi_ad(i-1,j+1)+FItoNeig; % outside of grid
             fi_ad(i,j-1)=fi_ad(i,j-1)+FItoNeig;
             fi_ad(i+1,j)=fi_ad(i+1,j)+FItoNeig;
             %fi_ad(i,j+1)=fi_ad(i,j+1)+FItoNeig;     % outside of grid
             fi_ad(i,j)=fi_ad(i,j)+2*FItoNeig; % reflective boundary conditions (dV/di=0)
             
         end
         i=M; % only 3 neighbors
         VFtoNeig=(Dv/6)*vf(i,j);
         FItoNeig=(Di/6)*free_interferon(i,j);
         
         vf_ad(i-1,j-1)=vf_ad(i-1,j-1)+VFtoNeig;
         vf_ad(i-1,j)=vf_ad(i-1,j)+VFtoNeig;
         %vf_ad(i-1,j+1)=vf_ad(i-1,j+1)+VFtoNeig; % outside of grid
         vf_ad(i,j-1)=vf_ad(i,j-1)+VFtoNeig;
         %vf_ad(i+1,j)=vf_ad(i+1,j)+VFtoNeig;     % outside of grid
         %vf_ad(i,j+1)=vf_ad(i,j+1)+VFtoNeig;     % outside of grid
         vf_ad(i,j)=vf_ad(i,j)+3*VFtoNeig; % reflective boundary conditions (dV/di=0)
         
         fi_ad(i-1,j-1)=fi_ad(i-1,j-1)+FItoNeig;
         fi_ad(i-1,j)=fi_ad(i-1,j)+FItoNeig;
         %fi_ad(i-1,j+1)=fi_ad(i-1,j+1)+FItoNeig; % outside of grid
         fi_ad(i,j-1)=fi_ad(i,j-1)+FItoNeig;
         %fi_ad(i+1,j)=fi_ad(i+1,j)+FItoNeig;     % outside of grid
         %fi_ad(i,j+1)=fi_ad(i,j+1)+FItoNeig;     % outside of grid
         fi_ad(i,j)=fi_ad(i,j)+3*FItoNeig; % reflective boundary conditions (dV/di=0)
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %
         % Even columns
         %
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         for j=2:2:(N-1)
             i=1; %(only 5 neighbors)
             VFtoNeig=(Dv/6)*vf(i,j);
             FItoNeig=(Di/6)*free_interferon(i,j);
             
             vf_ad(i,j-1)=vf_ad(i,j-1)+VFtoNeig;
             %vf_ad(i-1,j)=vf_ad(i-1,j)+VFtoNeig; % outside of grid
             vf_ad(i,j+1)=vf_ad(i,j+1)+VFtoNeig;
             vf_ad(i+1,j-1)=vf_ad(i+1,j-1)+VFtoNeig;
             vf_ad(i+1,j)=vf_ad(i+1,j)+VFtoNeig;
             vf_ad(i+1,j+1)=vf_ad(i+1,j+1)+VFtoNeig;
             vf_ad(i,j)=vf_ad(i,j)+1*VFtoNeig; % reflective boundary conditions (dV/di=0)
             
             fi_ad(i,j-1)=fi_ad(i,j-1)+FItoNeig;
             %fi_ad(i-1,j)=fi_ad(i-1,j)+FItoNeig; % outside of grid
             fi_ad(i,j+1)=fi_ad(i,j+1)+FItoNeig;
             fi_ad(i+1,j-1)=fi_ad(i+1,j-1)+FItoNeig;
             fi_ad(i+1,j)=fi_ad(i+1,j)+FItoNeig;
             fi_ad(i+1,j+1)=fi_ad(i+1,j+1)+FItoNeig;
             fi_ad(i,j)=fi_ad(i,j)+1*FItoNeig; % reflective boundary conditions (dV/di=0)
             
             for i=2:(M-1) % Inner cells (not in the boundary)
                 VFtoNeig=(Dv/6)*vf(i,j);
                 FItoNeig=(Di/6)*free_interferon(i,j);
                 
                 vf_ad(i,j-1)=vf_ad(i,j-1)+VFtoNeig;
                 vf_ad(i-1,j)=vf_ad(i-1,j)+VFtoNeig;
                 vf_ad(i,j+1)=vf_ad(i,j+1)+VFtoNeig;
                 vf_ad(i+1,j-1)=vf_ad(i+1,j-1)+VFtoNeig;
                 vf_ad(i+1,j)=vf_ad(i+1,j)+VFtoNeig;
                 vf_ad(i+1,j+1)=vf_ad(i+1,j+1)+VFtoNeig;
                 
                 fi_ad(i,j-1)=fi_ad(i,j-1)+FItoNeig;
                 fi_ad(i-1,j)=fi_ad(i-1,j)+FItoNeig;
                 fi_ad(i,j+1)=fi_ad(i,j+1)+FItoNeig;
                 fi_ad(i+1,j-1)=fi_ad(i+1,j-1)+FItoNeig;
                 fi_ad(i+1,j)=fi_ad(i+1,j)+FItoNeig;
                 fi_ad(i+1,j+1)=fi_ad(i+1,j+1)+FItoNeig;
                 
             end
             i=M; %(only 3 neighbors)
             VFtoNeig=(Dv/6)*vf(i,j);
             FItoNeig=(Di/6)*free_interferon(i,j);
             
             vf_ad(i,j-1)=vf_ad(i,j-1)+VFtoNeig;
             vf_ad(i-1,j)=vf_ad(i-1,j)+VFtoNeig;
             vf_ad(i,j+1)=vf_ad(i,j+1)+VFtoNeig;
             %vf_ad(i+1,j-1)=vf_ad(i+1,j-1)+VFtoNeig;   % outside of grid
             %vf_ad(i+1,j)=vf_ad(i+1,j)+VFtoNeig;       % outside of grid
             %vf_ad(i+1,j+1)=vf_ad(i+1,j+1)+VFtoNeig;   % outside of grid
             vf_ad(i,j)=vf_ad(i,j)+3*VFtoNeig; % reflective boundary conditions (dV/di=0)
             
             fi_ad(i,j-1)=fi_ad(i,j-1)+FItoNeig;
             fi_ad(i-1,j)=fi_ad(i-1,j)+FItoNeig;
             fi_ad(i,j+1)=fi_ad(i,j+1)+FItoNeig;
             %fi_ad(i+1,j-1)=fi_ad(i+1,j-1)+FItoNeig;   % outside of grid
             %fi_ad(i+1,j)=fi_ad(i+1,j)+FItoNeig;       % outside of grid
             %fi_ad(i+1,j+1)=fi_ad(i+1,j+1)+FItoNeig;   % outside of grid
             fi_ad(i,j)=fi_ad(i,j)+3*FItoNeig; % reflective boundary conditions (dV/di=0)
             
         end
         %
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         vf=vf_ad; %Update the virus free matrix info
         vf=(1-clearance_proportion)*vf; %viral clearance
         free_interferon=fi_ad; %update the interferon matrix info
         
         %Infection and immunity processes
         for i=1:M
             for j=1:N
                 if A_healthy(i,j)==1 %Meaning the cell is healthy and susceptible.
                     if inf_rate*vf(i,j)>rand %then...infect         
                    %if inf_rate*vf(i,j)> 0.5     
                         A_healthy(i,j)=0;
                         A_eclipse(i,j)=1;
                         status_interferon(i,j)=1;
                     end
                 end
                 if A_healthy(i,j)==1
%                     if immunity_rate*free_interferon(i,j)> 0.5
                     if immunity_rate*free_interferon(i,j)> rand %Means the cell becomes immune
                         A_healthy(i,j)=0; % the cell is no longer susceptible. It becomes immune
                         A_immune(i,j)=1; % The cell is now immune
                     end
                 end
             end
         end
         
         % In here, I have to update the values for cells in the eclipse and
         %secreting stages.
         
         for i=1:M
             for j=1:N
                 if A_eclipse(i,j)==1 % Check inner clock to jump to the secreting stage
                     if Eclipse_counter(i,j) > ticks_eclipse
                         A_eclipse(i,j)=0; %Remove from the eclipse class
                         A_secreting(i,j)=1; % Move to the secreting class
                     else
                         Eclipse_counter(i,j)=Eclipse_counter(i,j)+1;
                     end
                 end
                 if A_secreting(i,j)==1 %Check inner clock to jump to the death stage
                     %First, produce virus to be released into the medium
                     vf(i,j)=vf(i,j)+sec_rate;
                     if Secreting_counter(i,j) > ticks_secreting
                         A_secreting(i,j)=0; %Remove from the secreting class
                         A_death(i,j)=1; %The cell has died out
                         status_interferon(i,j)=3; %No more secretion of interferon
                     else
                         Secreting_counter(i,j)=Secreting_counter(i,j)+1;
                     end
                 end
                 if status_interferon(i,j)==1 % infected but not yet secretes interferon
                     if interferon_delay_counter(i,j) > ticks_interferon_delay
                         status_interferon(i,j)=2;
                     else
                         interferon_delay_counter(i,j)=interferon_delay_counter(i,j)+1;
                     end
                 elseif status_interferon(i,j)==2 %Interferon secretion stage
                     free_interferon(i,j)=free_interferon(i,j)+interferon_sec_rate;
                 end
             end
         end
         
         
         VirusTot(iTick)=sum(sum(vf));
         IFNTot(iTick)=sum(sum(free_interferon));
         HealthyCells(iTick)=sum(sum(A_healthy));
         EclipseCells(iTick)=sum(sum(A_eclipse));
         SecretingCells(iTick)=sum(sum(A_secreting));
%          DeathCells(iTick)=sum(sum(A_death));
         ImmuneCells(iTick)=sum(sum(A_immune));
         
         if iTick==time_data(indice)
             % count cells in each stage
             Cells_healthy=HealthyCells(iTick);
             Cells_eclipse=EclipseCells(iTick);
             Cells_secreting=SecretingCells(iTick);
             Cells_immune=ImmuneCells(iTick);
             Total_FreeVirus=VirusTot(iTick);

             errorCells=errorCells+abs(per_infect_cells(indice)-(Cells_eclipse+Cells_secreting)./...
                 (Cells_healthy+Cells_eclipse+Cells_secreting+Cells_immune));
             errorVirus=errorVirus+abs(virus_free_data_log_avg(indice)-log10(Total_FreeVirus));
             indice=indice+1;
             %[indice,errorCells,errorVirus]
             %pause
%              if indice>4
%                  break
%              end
         end
         
         INF_RATE(iRun)=inf_rate;
         SEC_RATE(iRun)=sec_rate;
         %INTERF_SEC_RATE(iRun)=interferon_sec_rate;
         IMMUNITY_RATE(iRun)=immunity_rate;
%          CLEARANCE_RATE(iRun)=clearance_rate;         
         ERROR_CELLS(iRun)=errorCells;
         ERROR_VIRUS(iRun)=errorVirus;
     end
     
     
     VirusTotMat(iRun,:)=VirusTot;
     IFNTotMat(iRun,:)=IFNTot;
     HealthyCellsMat(iRun,:)=HealthyCells;
     EclipseCellsMat(iRun,:)=EclipseCells;
     SecretingCellsMat(iRun,:)=SecretingCells;
     ImmuneCellsMat(iRun,:)=ImmuneCells;
     
end

trajectories.inf_rate = INF_RATE;
trajectories.sec_rate = SEC_RATE;
%trajectories.interferon_sec_rate = INTERF_SEC_RATE;
trajectories.immunity_rate=IMMUNITY_RATE;
% trajectories.clearance_rate=CLEARANCE_RATE;
trajectories.errorcells = ERROR_CELLS;
trajectories.errorvirus = ERROR_VIRUS;
trajectories.virusfree=VirusTotMat;
trajectories.interferon=IFNTotMat;
trajectories.healthycells=HealthyCellsMat;
trajectories.eclipsecells=EclipseCellsMat;
trajectories.secretingcells=SecretingCellsMat;
trajectories.immunecells=ImmuneCellsMat;

save trajectories_MODEL_III.mat trajectories

tiempototalproceso=cputime-tinitprocess;
disp(tiempototalproceso)