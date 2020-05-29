% 2-layer PFC model with variable populations and Poisson-based inputs.
% Script created by JSS on 11-Apr-2016, contact: sherfey@bu.edu
% For network model details, see 'help get_PFC_1layer'
% For cell model details, see 'help get_PFC_cell'

Ne=8;       % number of E-cells per layer
Ni=.25*Ne;  % number of I-cells per layer (FS + RSNP)

% PY taum is slower in superficial layers (L2/3 25+/-20ms vs L5/6 14+/-7ms)
% L2/3: http://www.neuroelectro.org/neuron/110/
% L5/6: http://www.neuroelectro.org/neuron/111/
% Tip: increase Cm or Rm=1/gleak in superficial layer wrt deep layer

% --------------------------------------------------------
% TODO: use ProbeCellProperties and CalcCellProperties to check taum in 
% L1 PY (target: 26+/-20ms) and L2 PY (target: 14+/-7ms). May need to 
% decrease Cm in L2 instead of increasing it in L1. See Example 3 in
% get_PFC_cell.m for instructions on calculating PY taum.
% --------------------------------------------------------

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Superficial layer (PY,FS,RSNP): two-compartment pyramidal cells (Es,Ed) and one-compartment interneurons (FS and RSNP)
L1=get_PFC_1layer('DS02PYjs',Ne,'DS02FSjs',Ni/2,'DS02RSNPjs',Ni/2);
% note: PY model has two compartments: soma (Es) and dendrite (Ed); single-compartment interneuron populations are named 'FS' and 'RSNP'
% get default value for membrane capacitance in PY soma
Cm_index=2*find(strcmp('Cm',L1.populations(1).parameters(1:2:end))); % index to 'Cm' parameter in PY soma (first population in specification L1)
Cm_soma=L1.populations(1).parameters{Cm_index}; % Cm value in PY soma
% specify connectivity matrix for PY->PY (can be used to specify network structure)
Kee=ones(Ne)-eye(Ne); % PY->PY connectivity matrix; e.g., all-to-all excluding self connections, [N_pre x N_post]
% update PY Cm and connection parameters for superificial layer model
modifications={...
  'Es','Cm',1.5*Cm_soma;      % increase by 50% the value of Cm soma in superficial layer
  'Ed','Cm',1.5*Cm_soma*1.92; % preserve relationship b/w soma and dendritic capacitance in [DS02]
  '(RSNP->FS','gGABA',1;      % RSNP inhibit FS cells (true for CR+ RSNP cells, less evidence for CB+ RSNP cells)
  'FS->RSNP','gGABA',0;       % FS do not inhibit RSNP cells
  '(FS->FS,RSNP->RSNP)','gGABA',0; % remove within-population inhibition
  'Es->Ed','netcon',Kee       % set connectivity matrix from PY soma to dendrite
  };
L1=dsApplyModifications(L1,modifications);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deep layer: PY and FS cells
L2=get_PFC_1layer('DS02PYjs',Ne,'DS02FSjs',Ni,[],0);
modifications={...
  'FS->FS','gGABA',0 % remove within-population inhibition
  };
L2=dsApplyModifications(L2,modifications);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepend layer-distinguising layer_id to names in populations and connections
layer_id='L1';
for i=1:length(L1.populations) 
  oldname=L1.populations(i).name;
	L1=dsApplyModifications(L1,{oldname,'name',[layer_id oldname]});
end
layer_id='L2';
for i=1:length(L2.populations) 
  oldname=L2.populations(i).name;
	L2=dsApplyModifications(L2,{oldname,'name',[layer_id oldname]});
end

% Combine network specifications to form a two-layer model
spec=[];
spec.populations=cat(2,L1.populations,L2.populations);
spec.connections=cat(2,L1.connections,L2.connections);

% Specify feedforward connectivity matrix from superficial to deep pyramidal cells
Kee=ones(L1.populations(1).size,L2.populations(1).size);  % [N_pre x N_post]
% Add feedforward connections to model specification:
spec.connections(end+1).direction='L1Es->L2Ed'; % source -> target
spec.connections(end).mechanism_list={'iAMPA'}; % connection mechanisms
spec.connections(end).parameters={'gAMPA',3e-3,'netcon',Kee}; % set max synaptic conductance and connectivity matrix

% % Add Poisson-based inputs to superficial PY dendrite (20kHz w/gAMPA=1e-3uS)
% generic input and state equations
input_def={'input(V)=iAMPA(V); monitor input; onset=50; offset=inf;';
           'iAMPA(V)=-gAMPA.*sAMPA(k,:).*(V-EAMPA); EAMPA=0; gAMPA=0;';
           'sAMPA=getPoissonGating(0,dcAMPA,0,0,0,onset,offset,tauAMPA,T,Npop); dcAMPA=0; tauAMPA=2;';
          }; % {'gAMPA',1e-3 to 1e-5,'dcAMPA',20e3,'Iapp',.1}; % [gAMPA]=uS, [dcAMPA]=Hz
state_equations=['dV/dt=(@current+input(V)+Iapp*(t>onset&t<offset))./Cm; Cm=1; Iapp=0; V(0)=-65;' input_def{:}];
spec=dsApplyModifications(spec,{'L1Ed','equations',state_equations});

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run simulations

% 1. tonic applied current to PY soma
% vary={'L1Es','Iapp',[.1 .2];'L1Es','gAMPA',0;'L1Es','dcAMPA',0;'L1Es','offset',1500};
% 2. poisson-based AMPAergic input to PY dendrite
vary={'L1Ed','Iapp',0;'L1Ed','gAMPA',1e-4;'L1Ed','dcAMPA',20e3;'L1Ed','offset',1000};
% vary={'(L1Es->L1Ed,L2Es->L2Ed)','gNMDA',[5e-4 5e-3 5e-2 5e-1];'(L1Es->L1Ed,L2Es->L2Ed)','gAMPA',0;'L1Ed','Iapp',0;'L1Ed','gAMPA',1e-4;'L1Ed','dcAMPA',20e3;'L1Ed','offset',1000};
solver_options={'tspan',[0 1500],'solver','rk1','dt',.01,'compile_flag',1,'verbose_flag',1};
data=dsSimulate(spec,'vary',vary,solver_options{:});

dsPlot(data,'plot_type','rastergram','threshold',-10);
dsPlot(data,'plot_type','waveform');


