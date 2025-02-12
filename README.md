R codes and files for fitting the self-exciting model proposed in the paper *A self-exciting spatio-temporal model with a smooth space-time varying productivity parameter*, by Álvaro Briz-Redón and Jorge Mateu.

There is a single main file, *Simulation study.R*, that allows:

1 - Simulating patterns from a self-exciting spatio-temporal model with a space-time varying offspring function (see Algorithm 1 in the paper) according to the three scenarios considered in the paper (the user can easily modify the spatio-temporal window or the offspring parameter function).

2 - Fitting a standard self-exciting model and the proposed self-exciting model with a smooth space-time varying productivity parameter to the simulated patterns. For this, the NIMBLE codes included in folder *Model codes* are needed.
