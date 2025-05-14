# HyperTraPS(-CT)

Hypercubic transition path sampling: Flexible inference of accumulation pathways, in discrete or continuous time, under different model structures, using combinations of longitudinal, cross-sectional, and phylogenetically-linked observations (Aga et al., 2024).

R package branch: install with
`remotes::install_github("StochasticBiology/hypertrapsct")`
then load with
`library(hypertrapsct)`

If you're looking for the codebase originally associated with the Aga et al. paper, it's here https://github.com/StochasticBiology/hypertraps-ct

![image](https://github.com/StochasticBiology/hypertraps-ct/assets/50171196/2c0fac84-76bf-41a6-9688-a4e429efed20)
An example inferred hypercubic transition graph (right) showing probable transitions during the evolution of multidrug resistance in tuberculosis, using phylogenetically-embedded original data (left) (Casali et al. 2014).

General content
=========

Requirements
------

HyperTraPS-CT is an R package with dependencies: `Rcpp`, `ggplot2`, `ggpubr`, `ggraph`, `ggwordcloud`, `igraph`, `stringr`, `stringdist`, `phangorn`, `phytools`, `ggtree`, `parallel`.

Demonstration
-----------
A good place to start is `demos/hypertraps-demos.R`, where the basic form of R commands for HyperTraPS, and most of the more interesting arguments that can be provided, are demonstrated. 

Input
------

HyperTraPS deals fundamentally with *transitions between states labelled by binary strings*. This allows us to work with cross-sectional, longitudinally, and phylogenetically-coupled samples.

The fundamental data element that goes into HyperTraPS is an observed transition from a "before" state to an "after" state. In the case of cross-sectional data, the "before" state is assumed to be the state of all zeroes (0000...), corresponding to a state which has not acquired any features. For longitudinal and/or phylogenetic observations, "before" and "after" states must be specified.

HyperTraPS requires at least a matrix describing "after" states -- this is a required argument. 

If a matrix supplying "before" states is absent, the data are assumed to be cross-sectional. For example, the matrix

`0 0 1`  
`0 1 1`

would reflect cross-sectional observations of states 001 and 011, implicitly corresponding to transitions 000->001 and 000->011.

A matrix of "before" states may be specified as a matrix using `initialstates`. For example, including the initial states

`0 0 1`  
`0 0 1`

with the above observations would now reflect the transitions 001->001 (i.e. remaining in state 001) and 001->011.

If you have phylogenetic data, we can curate it into transition format. `curate.tree` takes two arguments: a tree describing a phylogeny and a dataframe describing the features for each tip. These can be provided either as a rooted tree and a dataframe, or a filename for a tree in Newick format, and a filename for a CSV datafile to be read. The dataframe/CSV file should have labels in its first column that correspond to tip labels in the Newick tree, and the subsequent columns should give the L binary features for that tip. The list returned by `curate.tree` contains `srcs`, `dests`, and `times` elements which can be passed to HyperTraPS; the tree and data can be plotted with `plotHypercube.curated.tree`. There's an example in the tuberculosis case study below.

*If you're not interested in continuous time, uncertain data, priors, or old data formats, skip to the next section now.* 

For continuous-time inference, HyperTraPS works with a time window for each observed transition, specified via a start time and an end time. If the start time and end time are equal, the transition is specified as taking exactly that time. If start time = 0 and end time = Inf, the transition can take any amount of time, which is mathematically equivalent to the case without continuous time. In general, the start and end times specify an allowed set of durations for the given transition, allowing uncertain timings to be accounted for.

These start and end times are vectors specified by `starttimes` and `endtimes`. In both cases, absent start times means that all start times are assumed to be zero; absent end times means that all end times are assumed to be Inf.

The digit 2 can be used to reflect uncertainty in a state. For example,

`0 2 1`

corresponds to an observation where the first feature is absent, the second feature may be present or absent, and the third feature is present.

HyperTraPS also accepts descriptions of prior distributions on parameters. For the moment these are assumed to be uniform in log parameter space, and are specified by the min and max for each distribution. These should be provided as a matrix `priors` with two columns and N rows. Remember that N, the number of parameters, will depend on the model structure chosen and the number of features in the system. For example, the 3-feature system above and the default model structure (2, referring to L^2 parameters) would have N = 9.


Arguments to HyperTraPS
------

HyperTraPS needs at least a set of observations. *This is the only essential input.* This should take the form of a matrix. The table below gives other inputs that can be provided, with more technical points appearing in lower sections.

| Argument | R | Default | Notes |
|----------|---|---------|-------|
| Input data | obs=*matrix* |  None (required) | |
| Timings and initial states:||||
| Precursor states | initialstates=*matrix* | None ||
| Time window start | starttimes=*vector* |  0 ||
| Time window end | endtimes=*vector* | starttimes if present (i.e. precisely specified times); otherwise Inf ||
| More technical content: ||||
| Model structure | model=*N* |  2 | -1: every transition in dependently parameterised. 1: every feature independently acquired. 2: each feature acquisition independently influence other feature rates |(pairwise). 3: pairs of feature acquisitions can non-additively influence other acquisitions. 4: triplets of acquisitions can non-additively influence other acquisitions. |
| Prior mins and maxs | priors=*matrix* |  -10 to 10 in log space for each parameter (i.e. very broad range over orders of magnitude)|Specified as an n x 2 matrix, where n is number of parameters for the chosen model and columns give upper and lower bounds of a uniform distribution in log space.|
| Number of walkers | walkers=*N* |  200 | See "Assessing performance" |
| Inference chain length | length=*N* |  3 | Length will be 10^value (so 3 -> 10^3) |
| Perturbation kernel | kernel=*N* | 5 | See "Assessing performance". 1: with probability 0.1 per parameter, apply kernel width 0.005. 2-7: apply kernel to each parameter, with width 0.01 (2), 0.05 (3), 0.1 (4), 0.25 (5), 0.5 (6), 0.75 (7) |
| Random seed | seed=*N* |  1 ||
| Gains (0) or losses (1) | losses=*N* |  0 ||
| Use APM (0/1) | apm=*N* |  0 | Auxiliary pseudo-marginal MCMC (see Greenbury et al. 2020). Attempts to account for the random nature of the likelihood sampler by conditioning on random seed, which also varies throughout the MCMC run. Computationally more involved, but may prevent chains getting stuck.|
| Use SA (0/1) | sa=*N* |  0 | Simulated annealing. Rather than running MCMC, optimises parameter set against a decreasing temperature profile (inverse quarter power of the step). Will give a final point parameter estimate rather than a posterior distribution.|
| Use SGD (0/1) | sgd=*N* |  0 | Stochastic gradient descent. Rather than running MCMC, estimate likelihood gradient at each step and use this (scaled) to propose the next move. Will give a final point parameter estimate rather than a posterior distribution.|
| Scaling factor for gradient in SGD | sgd_scale=*N* |  0.01 | If using SGD, how much to scale the local likelihood gradient to determine the next parameter step.|
| Use PLI (0/1) | pli=*N* |  0 |  Phenotype landscape inference (see Williams et al. 2014). PLI applies no bias to the sampling walkers, reducing efficiency. Included for back-compatibility. |
| Gap between posterior samples | samplegap=*N* | 1e3 (>1e4 steps), 1e2 (>1e2 steps), or 1 | Number of MCMC steps between parameter samples that are treated as posterior samples. Higher values guard against correlated sampling.|
| Regularisation by penalised likelihood | penalty=*X* |  0 | Penalty per non-zero parameter |
| Regularisation by LASSO | lasso=*X* |  0 | LASSO penalty scaling summed absolute parameter values |
| Number of simulations per parameter sample | samples_per_row=*N* | 10 | Number of samples to use for each parameter set when simulating routes and times for output |
| Stepwise regularise model after best parameterisation (0/1) | regularise=*N* |  0 | |
| Output exact transitions (0/1) | output_transitions=*N* |  Switched off for L > 15 to avoid large output; consider switching off for CT ||
| Limit console output | limited_output=*N* |  0| |

So some example calls are (see the various demo scripts for more):

| Task | R | 
|------|---|
| Run HyperTraPS with default settings | HyperTraPS(*matrix*)  |
| Run HyperTraPS-CT with default settings | HyperTraPS(*matrix*, starttimes=*vector*, endtimes=*vector*) | 
| Run HyperTraPS with all-edges model and regularise by penalised likelihood | HyperTraPS(*matrix*, model=-1, penalty=1) | 
| Run HyperTraPS with all-edges model, then stepwise regularise | HyperTraPS(*matrix*, model=-1, regularise=1) | 

Visualising and using output
--------

The various outputs of HyperTraPS can be used in the R plotting functions below, which summarise (amongst other things) the numerical behaviour of the inference processes, the ordering and timings (where appropriate) of feature acquisitions, the structure of the learned hypercubic transition network, and any outputs from regularisation. *To start, `plotHypercube.summary` gives an overview of what we've learned.*

Probably the first thing to check is the likelihood trace (present in `plotHypercube.summary`, also directly via `plotHypercube.lik.trace`) -- see "Assessing performance" below.

All of these except `plotHypercube.prediction` and `plotHypercube.curated.tree` take a fitted model -- the output of `HyperTraPS` -- as a required argument, and may take other options as the table describes (look at the documentation for each function for a fuller description). 

| Plot function | Description | Options and defaults |
|---------------|-------------|---------|
| `plotHypercube.summary` | Summary plot combining several of the above | *f.thresh*=0.05 (flux threshold for graph plot), *t.thresh*=20 (time threshold for time histograms), *continuous.time*=TRUE (plot continuous time summary information) |
| More specific plots: |||
| `plotHypercube.lik.trace` | Trace of likelihood over inference run, re-calculated twice with different samples (to show consistency or lack thereof), along with current "in use" likelihood | |
| `plotHypercube.bubbles` | "Bubble plot" of probability of acquiring trait *i* at ordinal step *j* | *transpose*=FALSE (horizontal and vertical axis), *reorder*=FALSE (order traits by mean acquisition ordering) |
| `plotHypercube.motifs` | Motif-style plot of probability of acquiring trait *i* at ordinal step *j* |  |
| `plotHypercube.motifseries` | Motif-style plot of probability of specific states at a set of given snapshot times | *t.set*=0 (a set of snapshot times); *thresh*=0.05 (minimum probability for a state to be labelled) |
| `plotHypercube.graph` | Transition graph with edge weights showing probability flux (from full output) | *thresh*=0.05 (minimum threshold of flux for drawing an edge), *node.labels*=TRUE (show state labels on nodes), *node.label.size*=2 (size of those labels), *node.labels.box*=FALSE (draw opaque box around labels) |
| `plotHypercube.sampledgraph2` | Transition graph with edge weights showing probability flux (from sampled paths), with mean and s.d. of absolute timings for each step | *thresh*=0.05 (minimum threshold of flux for drawing an edge), *max.samps*=1000 (maximum number of sampled routes to consider), *no.times*=FALSE (avoid annotating edges with time information), *small.times*=FALSE (include alternative, smaller, offset time labels), *times.offset*=c(0.1,-0.1) (offset for those labels), *use.arc*=TRUE (arc edge format -- looks messier but less prone to overlapping edge labels), *node.labels*=TRUE (binary labels for nodes), *edge.label.size*=2 (font size for edge labels), *edge.label.angle*="across" (angle of edge labels), *edge.label.colour*="#000000# (edge label colour), *edge.check.overlap*=TRUE (avoid label overlaps), *featurenames*=c("") (set of feature names), *truncate*=-1 (truncate graph a given number of steps from root, -1 = don't), *use.timediffs*=TRUE (label with timnes for each transition, not overall time since start) |
| `plotHypercube.timehists` | Histograms of absolute timings for each trait's acquisition | *t.thresh*=20 (threshold time for x-axis), *featurenames*=c("") (set of feature names), *log.time*=TRUE (logarithmic time axis) |
| `plotHypercube.regularisation` | Information criterion vs number of nonzero parameters during regularisation | |
| `plotHypercube.motifs` | Motif plot of feature acquisition probabilities at discrete orderings | *featurenames*=c("") (set of feature names), *label.size*=3 (size of feature labels), *label.scheme*="full" (feature labels every timestep, or more sparsely) |
| `plotHypercube.timeseries` | Time series of acquisitions across sampled routes | *featurenames*=c("") (set of feature names), *log.time*=TRUE (logarithmic time axis) |
| `plotHypercube.motifseries` | Motif plot of state probabilities at a set of given times | *t.set*=0 (set of observation times), *thresh*=0.05 (minimum probability to explicitly label state) |
| `plotHypercube.prediction` | Visualise predictions of unobserved features or future behaviour, given a model fit | *prediction* (required, the output of `predictHiddenVals` or `predictNextStep` (see below)), *max.size*=30 (maximum size for word cloud) |
| `plotHypercube.influences` | For the L^2 model, visualise how each feature acquisition influences the rate of acquisition of other features as a matrix |  *featurenames*=c("") (set of names for features); *use.regularised*=FALSE (use stepwise-regularised param set); *reorder*=FALSE (order features by base rate); *upper.right*=FALSE (control orientation of diagonal); *cv.thresh*=Inf (threshold posterior coefficient of variation, only plot interactions below this)|
| `plotHypercube.influencegraph` | For the L^2 or L^3 model, visualise how each feature acquisition influences the rate of acquisition of other features as a network |  as `plotHypercube.influences`, plus *label.size*=2 (size of node labels) |
| `plotHypercube.curated.tree` | For phylogenetic data curated with `curate.tree`, visualise tree and barcodes |  Takes an object returned by `curate.tree` |

Some useful ones are demonstrated by the `plotHypercube.summary` command. This should work for all model structures (it omits influence plots, which are only supported for L^2 and L^3 models):
![image](https://github.com/StochasticBiology/hypertraps-ct/assets/50171196/c70d69b9-8a79-4aae-ba1b-675c2cd8e0b8)

In addition, we can ask HyperTraPS to make predictions about (a) any unobserved features for a given observation (for example, what value the ?s might take in 01??), and (b) what future evolutionary behaviour is likely given that we are currently in a certain state. These are achieved with functions `predictHiddenVals` and `predictNextStep` respectively. Both of these require a fitted hypercubic model (the output of `HyperTraPS`), which we'll call `fit`.

To predict unobserved values in a given observation, you can invoke

`predictHiddenVals(fit, state)`

where `fit` is the fitted model and `state` is a vector describing the incomplete observations with 0s, 1s, and 2s, the latter of which mark unobserved features. So 0122 has uncertain features at positions 3 and 4. The function will output two dataframes, one giving the probability of observing each specific state that could correspond to the incomplete observation, and one giving the aggregate probability of each uncertain feature being 1.

You can optionally provide an argument `level.weights` which provides probability weightings for different "levels" of the hypercube, that is, the different total number of 1s in the observation. This allows you to specify how likely it is that the true observation has acquired a certain number of features. By default this is uniform between the number of certain 1s and the maximum number of possible 1s.

To predict future evolutionary behaviour, you can invoke

`predictNextStep(fit, state)`

where `fit` is the fitted model and `state` is a given state. This will return a dataframe describing the possible future dynamics from `state` and their probabilities.

You can pass the output of both these prediction functions to `plotHypercube.prediction` for a visualisation of the corresponding prediction.

HyperTraPS will also estimate the probability of a particular state for either the discrete or the continuous time cases. In the discrete case, this is returned in `fit$dynamics`, where the probabilities of different states and different transitions are reported. The state probabilities are reported in normalised form (assuming uniform sampling over all the L+1 states in a complete trajectory) and unnormalised (the probability of a state given that we have acquired that particular number of features).

In the continuous case, the function `prob.by.time` will use the sampled dynamics from the inferred parameterisation(s) to output a set of state probabilities at a given sampling time. For example

`prob.by.time(fit, 0.1)` 

will return a set of states and their associated probabilities at time 0.1.

To ask about the observation probability of given states, you can use

`state.probs(fit)`

which reports the probabilities of different states (i) assuming that exactly the number of features in each state has been acquired (e.g. P(00110 | system has acquired 2 features)) and (ii) conditional on a profile of probabilities of feature acquisitions (e.g. P(00110 | P(0 features acquired), P(1 feature acquired), ..., P(5 features acquired))).

To get the likelihood for an observation given a parameterisation, you can use

`getLikelihood(obs, params, model)`

where `obs` is a vector giving a state, `params` is a vector containing the parameter set, and `model` is the model structure (-1, 1, 2, 3, 4). 

Assessing performance
------

HyperTraPS uses sampling to estimate likelihoods, and optionally (and by default) uses these likelihoods in an MCMC chain. Several aspects of this process can pose challenges. We want the likelihood estimates to be reliable and consistent, and the MCMC chain to be equilibrated and well mixing. Departures from this picture can lead to systematic errors in the parameter estimates and downstream output. 

`plotHypercube.lik.trace` can be used to diagnose some of these issues. We want the output to look more like (D) in the figure below: the estimates (different lines) are overlapping, the chain is quite uncorrelated, there aren't systematic changes in the mean. In (A) the likelihood estimates don't agree: the black dashed lines and the red line are inconsistent. Here this has led to an MCMC issue as well: the red line undergoes long periods of stasis (it's been frozen by rare extreme samples). To fix non-converged estimates, use more random walkers in the estimation (increase `walkers`). You could also try auxiliary pseudo-marginal MCMC (`apm`). In (B) the chain is quite correlated, getting `stuck` at the same values for some time. This is because the kernel size for proposing MCMC moves is too big, so lots of proposed moves get rejected. To fix this, use a smaller kernel (decrease `kernel`). In (C) the chain is not equilibrated; there is still a systematic increase in the mean. To fix this, run a longer chain (increase `length`); you could also try a larger kernel (increase `kernel`) as long as the issue in (B) is avoided. 

![image](https://github.com/user-attachments/assets/ceb842f9-6cc3-4fcb-a70b-cdf008a0c0f6)


Output details
------

HyperTraPS will output information about the progress of a run to the console, and return a named list containing outputs from and descriptions of the inference process.

*If you are just interested in plotting summary outputs from the inference process, skip to the next section now.*

The output structures are

| Information | R | Notes |
|-------------|---|-------|
| Number of features, model type, number of parameters, and likelihood traces over the run | *list*$lik.traces | |
| Best parameterisation found during run | *list*$best | |
| (Posterior) samples of parameterisation | *list*$posterior.samples | |
| Probabilities for individual states | *list*$dynamics$states | Only produced for L < 16 |
| Probabilities for individual transitions | *list*$dynamics$trans |  Only produced for L < 16 |
| Best parameterisation after regularisation | *list*$regularisation$best |  Optional |
| Number of parameters, likelihood, and information criteria during regularisation | *list*$regularisation$reg.process | Optional |
| "Bubble" probabilities of trait *i* acquisition at ordinal time *j* | *list*$bubbles |  |
| Histograms of times of trait *i* acquisition | *list*$timehists |  |
| Individual sampled routes of accumulation | *list*$routes | Matrix with L columns; jth element of a row gives the feature acquired at step j. Each row is a sampled trajectory; there are *samples_per_row* samples per output parameterisation in the (posterior) sample set |
| Transition times for individual sampled routes of accumulation | *list*$times | Matrix with L columns, with elements corresponding to the timings of each of the steps in the routes matrix above |
| Dwelling statistics for individual sampled routes of accumulation | *list*$betas |  Matrix with L columns, with elements corresponding to the "beta" (characteristic rate of departure) for each step in the routes matrix above |


References
=====

Aga, O.N., Brun, M., Dauda, K.A., Diaz-Uriarte, R., Giannakis, K. and Johnston, I.G., 2024. HyperTraPS-CT: Inference and prediction for accumulation pathways with flexible data and model structures. PLOS Computational Biology, 20(9), p.e1012393.

Casali, N., Nikolayevskyy, V., Balabanova, Y., Harris, S.R., Ignatyeva, O., Kontsevaya, I., Corander, J., Bryant, J., Parkhill, J., Nejentsev, S. and Horstmann, R.D., 2014. Evolution and transmission of drug-resistant tuberculosis in a Russian population. Nature genetics, 46(3), pp.279-286.

Greenbury, S.F., Barahona, M. and Johnston, I.G., 2020. HyperTraPS: inferring probabilistic patterns of trait acquisition in evolutionary and disease progression pathways. Cell systems, 10(1), pp.39-51.

Williams, B.P., Johnston, I.G., Covshoff, S. and Hibberd, J.M., 2013. Phenotypic landscape inference reveals multiple evolutionary paths to C4 photosynthesis. elife, 2, p.e00961.

