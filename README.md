# Parameter-wise-prediction
This repository contains several Julia codes associated with "Parameter-wise predictions and sensitivity analysis for random walk models in the life sciences" by Yihan Liu, David J. Warne, and Matthew J. Simpson. The preprint is available at 

# Case1
 This file contains four sub-files: Data1, Data2, Analysis1 and Analysis2.
## Data1 from Case1
 DiscreteModel.jl is responsible for generating all synthetic data from discrete model at early time point k=100.
 
 solutions.jl is responsible for generating all data from PDE model at early time point k=100.
 
 Plot.jl is responsible for creating subfigures in Figure 2 and Figure S1 uisng data generated from DiscreteModel.jl and solutions.jl under the same file.
## Data2 from Case1
 DiscreteModel.jl is responsible for generating all synthetic data from discrete model at late time point k=200.
 
 solutions.jl is responsible for generating all data from PDE model at late time point k=200.
 
 Plot.jl is responsible for creating subfigures in Figure S4 uisng data generated from DiscreteModel.jl and solutions.jl under the same file.
## Analysis1 from Case1  
 Possion.jl is responsible for creating subfigures in Figure 4 using the count data generated in DiscreteModel.jl in Data1 (data.csv).
## Analysis2 from Case1  
 Possion.jl is responsible for creating subfigures in Figure S5 using the count data generated in DiscreteModel.jl in Data2 (data.csv).

# Case2
 This file contains four sub-files: Data1, Data2, Analysis1 and Analysis2.
## Data1 from Case2
 DiscreteModel.jl is responsible for generating all synthetic data from discrete model at early time point k=100.
 
 ContinuumModel.jl is responsible for generating all data from PDE model at early time point k=100.
 
 Plot.jl is responsible for creating subfigures in Figure 3 uisng data generated from DiscreteModel.jl and solutions.jl from both Data1 and Data2 in Case2.
## Data2 from Case2
 DiscreteModel.jl is responsible for generating all synthetic data from discrete model at late time point k=600.
 
 ContinuumModel.jl is responsible for generating all data from PDE model at late time point k=600.
 
 Plot.jl is responsible for creating subfigures in Figure S2 and Figure S3 uisng data generated from DiscreteModel.jl and solutions.jl from Data2.
## Analysis1 from Case2  
 Binomial.jl is responsible for creating subfigures in Figure 5 using the count data generated in DiscreteModel.jl in Data1 under Case2 (data.csv).
## Analysis2 from Case2  
 Binomial.jl is responsible for creating subfigures in Figure 6 using the count data generated in DiscreteModel.jl in Data2 under Case2 (data.csv).

