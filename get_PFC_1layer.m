function spec=get_PFC_1layer(PYtype,Ne,PVtype,Nfs,CBtype,Nrs)
% spec=get_PFC_1layer(PYtype,Ne,PVtype,Nfs,CBtype,Nrs)
% Purpose: retrieve deterministic single layer PFC network with variable
% pyramidal cell (PY) and 1-2 IN populations (PV+ FS and CB+ RSNP cells).
% Inputs:
%   Ne = number of pyramidal cells
%   PYtype = pyramidal cell model (default: 'DS02PYjs') (see get_PFC_cell for options)
%   Nfs = number of fast-spiking (FS) cells
%   PVtype = PV+ IN (FS) cell model (default: 'DS02FSjs') (see get_PFC_cell for options)
%   Nrs = number of regular-spiking non-pyramidal (RSNP) cells
%   CBtype = CB+ IN (RSNP or LTS) cell model (default: 'DS02RSNPjs') (see get_PFC_cell for options)
% Outputs:
%   spec: DynaSim specification structure
% 
% % Example: default network model (note: ~3min to compile first time only)
% solver_options={'tspan',[0 500],'solver','rk2','dt',.01,'compile_flag',1,'verbose_flag',1};
% data=SimulateModel(get_PFC_1layer,'vary',{'Es','Iapp',.1},solver_options{:});
% PlotData(data)
% 
% % Example: PY-FS-RSNP
% spec=get_PFC_1layer('DS02PYjs',8,'DS02FSjs',2,'DS02RSNPjs',2);
% data=SimulateModel(spec,'vary',{'Es','Iapp',.1},solver_options{:});
% PlotData(data);
% 
% % Example: PY-FS
% spec=get_PFC_1layer('DS02PYjs',8,'DS02FSjs',2,[],0);
% 
% % Example: default cell models with variable population sizes
% spec=get_PFC_1layer([],8,[],2,[],2);
% 
% See also: get_PFC_cell, PFC_Nlayers

% Note: connectivity differences b/w IN populations are based on PV+ and
% CB+ cells. Therefore, PV and CB may be better variable identifiers to use
% in this script (instead of FS and RSNP). 

% created by JSS on 11-Apr-2016, contact: sherfey@bu.edu

if nargin<1 || isempty(PYtype), PYtype='DS02PYjs'; end
if nargin<2, Ne=8; end
if nargin<3 || isempty(PVtype), PVtype='DS02FSjs'; end
if nargin<4, Nfs=2; end
if nargin<5 || isempty(CBtype), CBtype='DS02RSNPjs'; end
if nargin<6, Nrs=2; end

% Single layer PFC model size
% Ne=8;  % # pyramidal cells
% Nfs=2; % # fast-spiking (FS) cells
% Nrs=2; % # regular-spiking non-pyramidal (RSNP) cells

% -------------------------------------------------------------------------
% Specify populations: Custom cell models (see also: 'help get_PFC_cell')
% -------------------------------------------------------------------------
% Two-compartment Pyramidal cell model ('Es'=soma, 'Ed'=dendrite)
spec=get_PFC_cell(PYtype,Ne);
% One-compartment FS model ('FS' = PV+ interneuron inhibiting PY soma)
spec.populations(end+1)=getfield(get_PFC_cell(PVtype,Nfs),'populations');
% One-compartment RSNP model ('RSNP' = CB+ interneuron inhibiting PY dendrite)
spec.populations(end+1)=getfield(get_PFC_cell(CBtype,Nrs),'populations');

% -------------------------------------------------------------------------
% Specify network connections (all-to-all)
% -------------------------------------------------------------------------
% [DS02,DG07]: AMPA (taur=.2,taud=1,E=0), NMDA (taur=2.3,taud=95,E=0), GABA (taur=.5,taud=5,E=-75)
tauAMPAr=.2;  % ms, AMPA rise time
tauAMPAd=1;   % ms, AMPA decay time
tauNMDAr=2.3; % ms, NMDA rise time
tauNMDAd=95;  % ms, NMDA decay time
tauGABAr=.5;  % ms, GABAa rise time
tauGABAd=5;   % ms, GABAa decay time
% [DS02]:
Npc=100;  % # principal cells in original [DS02] publication
Nin=37;   % # interneurons in original [DS02] publication
gAMPAee=3e-3;       % uS, PY->PY, maximal AMPA conductance
gNMDAee=gAMPAee/50; % uS, PY->PY, maximal NMDA conductance
gGABAie=.2e-3;      % uS, FS->PY, maximal GABAa conductance
gAMPAei=.74e-3;     % uS, PY->IN
gNMDAei=gAMPAei/50; % uS, PY->IN
gGABAii=.6e-3;      % uS, FS->IN

% [D97] DeFelipe, J. (1997). Types of neurons, synaptic connections and chemical characteristics of cells immunoreactive for calbindin-D28K, parvalbumin and calretinin in the neocortex. Journal of chemical neuroanatomy, 14(1), 1-19.
% - NMDA colocalizes with PV+ but not CB+
% - PV+ INs target PY soma; CB+ INs target PY dendrite
% - PV+ INs are FS; CB+ INs are RSNP or LTS
% [Traub/Whittington] (2010, in "Cortical Oscillations in Health and Disease")
% - CB+ inhibition is longer-lasting than PV+ inhibition (reference needed for primary research)
tauGABAdCB=13; % ms, GABAa decay time for inhibition from RSNP CB+ interneurons

% specify connection kernels that exclude connections to self
Kee=ones(Ne)-eye(Ne);
Kfsfs=ones(Nfs)-eye(Nfs);
Krsrs=ones(Nrs)-eye(Nrs);

% recurrent connections between pyramidal cells (add to existing intercompartmental connections)
index=find(strcmp('Es->Ed',{spec.connections.direction}),1,'first');
spec.connections(index).mechanism_list={'iAMPA','iNMDA',spec.connections(index).mechanism_list{:}};
spec.connections(index).parameters={'gAMPA',gAMPAee*Npc/Ne,'gNMDA',gNMDAee*Npc/Ne,'EAMPA',0,'ENMDA',0,'netcon',Kee,...
  'tauAMPAr',tauAMPAr,'tauAMPA',tauAMPAd,'tauNMDAr',tauNMDAr,'tauNMDA',tauNMDAd,spec.connections(index).parameters{:}};
% pyramidal<->FS connections
spec.connections(end+1).direction='Es->FS';
spec.connections(end).mechanism_list={'iAMPA','iNMDA'};
spec.connections(end).parameters={'gAMPA',gAMPAei*Npc/Ne,'gNMDA',gNMDAei*Npc/Ne,'EAMPA',0,'ENMDA',0,...
  'tauAMPAr',tauAMPAr,'tauAMPA',tauAMPAd,'tauNMDAr',tauNMDAr,'tauNMDA',tauNMDAd};
spec.connections(end+1).direction='FS->Es';
spec.connections(end).mechanism_list={'iGABA'};
spec.connections(end).parameters={'gGABA',gGABAie*Nin/Nfs,'tauGABAr',tauGABAr,'tauGABA',tauGABAd,'EGABA',-75};
% pyramidal<->RSNP connections
spec.connections(end+1).direction='Es->RSNP';
spec.connections(end).mechanism_list={'iAMPA'};
spec.connections(end).parameters={'gAMPA',gAMPAei*Npc/Ne,'EAMPA',0,'tauAMPAr',tauAMPAr,'tauAMPA',tauAMPAd};
spec.connections(end+1).direction='RSNP->Ed';
spec.connections(end).mechanism_list={'iGABA'};
spec.connections(end).parameters={'gGABA',gGABAie*Nin/Nrs,'tauGABAr',tauGABAr,'tauGABA',tauGABAdCB,'EGABA',-75};
% interneuron<->interneuron connections
spec.connections(end+1).direction='FS->FS';
spec.connections(end).mechanism_list={'iGABA'};
spec.connections(end).parameters={'gGABA',gGABAii*Nin/Nfs,'tauGABAr',tauGABAr,'tauGABA',tauGABAd,'EGABA',-75,'netcon',Kfsfs};
spec.connections(end+1).direction='RSNP->RSNP';
spec.connections(end).mechanism_list={'iGABA'};
spec.connections(end).parameters={'gGABA',gGABAii*Nin/Nrs,'tauGABAr',tauGABAr,'tauGABA',tauGABAdCB,'EGABA',-75,'netcon',Krsrs};
spec.connections(end+1).direction='FS->RSNP';
spec.connections(end).mechanism_list={'iGABA'};
spec.connections(end).parameters={'gGABA',gGABAii*Nin/Nfs,'tauGABAr',tauGABAr,'tauGABA',tauGABAd,'EGABA',-75};
spec.connections(end+1).direction='RSNP->FS';
spec.connections(end).mechanism_list={'iGABA'};
spec.connections(end).parameters={'gGABA',gGABAii*Nin/Nrs,'tauGABAr',tauGABAr,'tauGABA',tauGABAdCB,'EGABA',-75};

% remove populations with size=0
spec=dsCheckSpecification(spec);

if 0
  % -------------------------------------------------------------------------
  % Run simulations
  % -------------------------------------------------------------------------

  % simulation controls (using rk2)
  tspan=[0 500];  % [beg end], ms
  dt=.01;         % fixed time step, ms
  compile_flag=1; % whether to compile simulation (note: takes ~minutes to compile the first time, but then provides 20x speed-up on subsequent simulations compared to uncompiled sims; eg., compiled: 13sec, uncompiled: 254sec)
  solver_options={'tspan',tspan,'solver','rk2','dt',dt,'compile_flag',compile_flag,'verbose_flag',1};

  % optionally specify things to vary in model and across simulations
  % (set vary=[] to simulate default model; see 'help Vary2Modifications' for more details)
  vary={'Es','Iapp',.1};

  % run simulation
  data=dsSimulate(spec,'vary',vary,solver_options{:});

  % plot results
  dsPlot(data);
end
