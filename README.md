# PFC Models
This repository contains a collection of neural models for prefrontal cortical (PFC) cells and networks.

Download the repo on [github](https://github.com/jsherfey/PFC_models) or using git: `git clone https://github.com/jsherfey/PFC_models.git`

All models are implemented in the DynaSim Matlab toolbox: https://github.com/DynaSim/DynaSim

To get started with individual PFC cell models, see [`PFC_cells`](https://github.com/jsherfey/PFC_models/blob/master/PFC_cells.m). To get started with PFC network models, see [`PFC_1layer`](https://github.com/jsherfey/PFC_models/blob/master/PFC_1layer.m) and [`PFC_2layers`](https://github.com/jsherfey/PFC_models/blob/master/PFC_2layers.m).

The PFC network model is based on a simpler one with pyramidal cells and fast-spiking interneurons by Durstewitz et al. Matlab implementations of 2002 and 2007 versions of the original model and references can be found in the "[Durstewitz](https://github.com/jsherfey/PFC_models/tree/master/Durstewitz)" folder in [`DS02_PFC_deep`](https://github.com/jsherfey/PFC_models/blob/master/Durstewitz/DS02_PFC_deep.m) and [`DG07_PFC_deep`](https://github.com/jsherfey/PFC_models/blob/master/Durstewitz/DG07_PFC_deep.m), respectively.

The model by Jason Sherfey consists of cells from the Durstewitz model as well as others. DynaSim implementations of each cell model, references, and supplemental info can be found in [`PFC_cells`](https://github.com/jsherfey/PFC_models/blob/master/PFC_cells.m). [`get_PFC_cell`](https://github.com/jsherfey/PFC_models/blob/master/get_PFC_cell.m) is a function that can be used to retrieve DynaSim specifications of individual cell models. Similarly, [`PFC_1layer`](https://github.com/jsherfey/PFC_models/blob/master/PFC_1layer.m) and [`get_PFC_1layer`](https://github.com/jsherfey/PFC_models/blob/master/get_PFC_1layer.m) are script and function forms of a single-layer PFC network model with pyramidal cells, PV+ FS cells, and CB+ RSNP cells. 

[`PFC_2layers`](https://github.com/jsherfey/PFC_models/blob/master/PFC_2layers.m) is a script that constructs a two-layer PFC model representing minimal superficial and deep layers.
`PFC_competition` ([work-in-progress](https://github.com/jsherfey/PFC_simulations/blob/master/PFC_competition.m)) defines two assemblies in a single-layer PFC network and probes how they compete when they are driven by inputs with varying rhythmicity and synchrony.

------------------------------------------------------------

To install DynaSim using Git: `git clone https://github.com/dynasim/dynasim.git`.
Add the dynasim toolbox to Matlab path and run model scripts from the directory containing them.

