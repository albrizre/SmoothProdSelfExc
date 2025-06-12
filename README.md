R codes and files for fitting the self-exciting model proposed in the paper *A self-exciting spatio-temporal model with a smooth space-time varying productivity parameter*, by Álvaro Briz-Redón and Jorge Mateu.

There is a main file, *Simulation study.R*, that allows:

1 - Simulating patterns from a self-exciting spatio-temporal model with a space-time varying offspring parameter function (see Algorithm 1 in the paper) according to the three scenarios considered in the paper (the user can easily modify the spatio-temporal window or the offspring parameter function).

2 - Fitting a standard self-exciting model and the proposed self-exciting model with a smooth space-time varying productivity parameter to the simulated patterns. For this, the NIMBLE codes included in folder *Model codes* are needed.

We note that, for a given pattern, the standard model needs to be fit first since the estimates of this model are used for building the priors in the proposed model.

Furthermore, file *Figures simulation study.R* allows reproducing the figures included in the paper about the simulation study. In particular:

a) Graphical representation of the offspring parameter functions, kappa(x,t), in space and time.

b) Graphical representation of the average estimate of b(x) and b(t) for the 100 replicates created under each of the three scenarios considered. For this, the models need to be fitted first and storaged in *Models/Sim_study/*. These are not provided in the Github repository due to their size.
