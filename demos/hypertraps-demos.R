#################

remotes::install_github("StochasticBiology/hypertraps-ct@bioconductor-tmp")
library(hypertrapsct)
library(ggplot2)
library(ggpubr)
library(phytools)

#####
## Example 1: simple cross-sectional data

m.2 = matrix(rep(c(1,0,0,0,0,
                   1,1,0,0,0,
                   1,1,1,0,0,
                   1,1,1,1,0,
                   1,1,1,1,1,
                   0,0,0,0,1,
                   0,0,0,1,1,
                   0,0,1,1,1,
                   0,1,1,1,1,
                   1,1,1,1,1),5), byrow=TRUE, ncol=5)

# simplest example: short chain on a simple dataset of cross-sectional observations
my.post = HyperTraPS(m.2)
plotHypercube.summary(my.post)

#####
## Example 2: phylogenetically related data

# now let's imagine those observations are embedded on a phylogeny
# create a dataframe containing the same data but with unique IDs for each observation
df.1 = data.frame(id = paste("id_", 1:nrow(m.2), sep=""))
df.1 = cbind(df.1, m.2)
# create a random tree with the same set of IDs for the tips. this could equally well be read in from a file
tree.1 = rtree(n = nrow(m.2))
tree.1$tip.label = df.1$id
# "curate.tree" reconstructs ancestral states based on rare, irreversible transitions and provides the input needed for HyperTraPS
ct.1 = curate.tree(tree.1, df.1)
plotHypercube.curated.tree(ct.1, names=TRUE)

# now do the inference accounting for the phylogenetic relationships
my.post.tree = HyperTraPS(ct.1$dests, initialstates = ct.1$srcs)
plotHypercube.summary(my.post.tree)

# now do the inference in continuous time accounting for the phylogenetic relationships, including branch lengths
my.post.tree.ct = HyperTraPS(ct.1$dests, initialstates = ct.1$srcs,
                             endtimes = ct.1$times*1000)
plotHypercube.summary(my.post.tree.ct)

#####
## Example 3: checking for convergence / comparing outputs

# three different inference runs with different random seeds
my.post.1 = HyperTraPS(m.2, length = 3,
                       samplegap = 10,
                       output_transitions = 1,
                       seed = 1,
                       featurenames = c("A", "B", "C", "D", "E")); 
my.post.2 = HyperTraPS(m.2, length = 3,
                       samplegap = 10,
                       output_transitions = 1,
                       seed = 2,
                       featurenames = c("A", "B", "C", "D", "E")); 
my.post.3 = HyperTraPS(m.2, length = 3,
                       samplegap = 10,
                       output_transitions = 1,
                       seed = 3,
                       featurenames = c("A", "B", "C", "D", "E")); 

# examine the likelihood traces. these should show:
# i. overlap between dashed black and red lines.
#.    -> if not, likelihood estimates aren't converging. use more random walkers.
# ii. constant mean likelihood (no systematic up or down)
#.    -> if not, MCMC chain isn't converged. run longer chains.
# iii. uncorrelated (white) noise; no periods of stasis
#.    -> if not, MCMC chain is getting stuck. use a smaller kernel.
ggarrange(plotHypercube.lik.trace(my.post.1),
          plotHypercube.lik.trace(my.post.2),
          plotHypercube.lik.trace(my.post.3))

# compare the "bubble plots" for the three random seeds
plotHypercube.bubbles.compare(list(my.post.1, my.post.2, my.post.3))

#####
## Example 4: longitudinal data

# now an artificial set of before-after states, for example emerging from longitudinal observations
m.1 = matrix(rep(c(0,0,0,0,0,
                   1,0,0,0,0,
                   1,1,0,0,0,
                   1,1,1,0,0,
                   1,1,1,1,0,
                   0,0,0,0,0,
                   0,0,0,0,1,
                   0,0,0,1,1,
                   0,0,1,1,1,
                   0,1,1,1,1),5), byrow=TRUE, ncol=5)
times.cs = rep(c(0.1, 0.2, 0.3, 0.4, 0.5), 10)
times = rep(0.1, 50)

my.post.l = HyperTraPS(m.2, initialstates = m.1, 
                       starttimes = times, endtimes = times, 
                       length = 3.5,
                       samplegap = 10,
                       output_transitions = 1,
                       featurenames = c("A", "B", "C", "D", "E")); 

#####
## Example 4: other visualisations and analysis

my.post = my.post.l
# likelihood trace for checking
plotHypercube.lik.trace(my.post)
# sampled transition graph
plotHypercube.sampledgraph2(my.post, thresh=0.1, use.arc=FALSE, edge.label.size=3) + 
  theme(legend.position="none") + expand_limits(x = c(-0.1, 1.1))
# influences between features
plotHypercube.influences(my.post, cv.thresh = Inf)
plotHypercube.influencegraph(my.post, cv.thresh = 1)
# state probabilities as a function of time
plotHypercube.motifseries(my.post, c(0.001, 0.01, 0.5, 1, 2, 5))

# demonstrate predictions of future behaviour
prediction.step = predictNextStep(my.post, c(1,1,0,0,0))
plotHypercube.prediction(prediction.step)

# demonstrate predictions of hidden values...
# ... with no assumptions about progress on the hypercube
prediction.hidden = predictHiddenVals(my.post, c(1,2,2,2,2))
plotHypercube.prediction(prediction.hidden)
# ... enforcing given belief about progress on the hypercube
prediction.hidden = predictHiddenVals(my.post, c(1,2,2,2,2), level.weight=c(0,0,1,0,0,0))
plotHypercube.prediction(prediction.hidden)

# impose priors -- here disallowing every pairwise effect
# prior format: n_param rows, 2 columns. [i,1] = min i; [i,2] = max i
# here we impose priors on the base rates: feature 1 > feature 5 >> all others
priors = matrix(0, ncol=2, nrow=5*5)
priors[,1] = -10
priors[,2] = 10
for(i in 0:4) {
  priors[i*5+i+1,1] = -10
  priors[i*5+i+1,2] = -10
}
priors[0*5+0+1,1] = 1
priors[0*5+0+1,2] = 1
priors[4*5+4+1,1] = 0
priors[4*5+4+1,2] = 0

my.post.priors = HyperTraPS(m.2, initialstates = m.1, 
                     starttimes = times, endtimes = times, 
                     priors = priors,
                     featurenames = c("A", "B", "C", "D", "E")); 
plotHypercube.summary(my.post.priors)
# compare to results from L model (independent features)  
my.post.model1 = HyperTraPS(m.2, initialstates = m.1, 
                            starttimes = times, endtimes = times, 
                            model = 1, kernel=3,
                            featurenames = c("A", "B", "C", "D", "E")); 
plotHypercube.summary(my.post.model1)


### different levels of uncertainty in timings
# precisely specified timings, as above
my.post.time.precise = HyperTraPS(m.2, initialstates = m.1, 
                               starttimes = times, endtimes = times, 
                               length = 3, outputinput = 1,
                               featurenames = c("A", "B", "C", "D", "E")); 
# infinite width time window for transitions (just inferring ordering)
my.post.time.inf = HyperTraPS(m.2, initialstates = m.1, 
                                starttimes = times*0, endtimes = times*Inf, 
                                length = 3, outputinput = 1,
                                featurenames = c("A", "B", "C", "D", "E"));
# finite time window for each uncertain transition time
my.post.time.uncertain = HyperTraPS(m.2, initialstates = m.1, 
                     starttimes = times*0.25, endtimes = times*4, 
                     length = 3, outputinput = 1,
                     featurenames = c("A", "B", "C", "D", "E")); 
ggarrange(plotHypercube.timehists(my.post.time.precise, t.thresh=3), 
          plotHypercube.timehists(my.post.time.uncertain, t.thresh=3),
          plotHypercube.timehists(my.post.time.inf, t.thresh=3),
          nrow=3)

# write output to files
writeHyperinf(my.post, "simpledemo", my.post$L, postlabel = "simpledemo", fulloutput=TRUE)

# retrieve output from files
my.post.r = readHyperinf("simpledemo", postlabel = "simpledemo", fulloutput=TRUE)
plotHypercube.summary(my.post.r)

# run an example with fewer walkers
my.post.sparse = HyperTraPS(m.2, initialstates = m.1, 
                            starttimes = times, endtimes = times,
                            featurenames = c("A", "B", "C", "D", "E"), walkers = 2); 
plotHypercube.summary(my.post.sparse, t.thresh = 2)

# direct time run (no time window specified)
my.post.dt = HyperTraPS(m.2, initialstates = m.1, featurenames = c("A", "B", "C", "D", "E")); 
plotHypercube.summary(my.post.dt, continuous.time = FALSE)
ggarrange(plotHypercube.timehists(my.post.dt, t.thresh=3),
          plotHypercube.timehists(my.post.time.inf, t.thresh=3),
          nrow=2)

### various other demos
# other plots
plotHypercube.motifs(my.post)
plotHypercube.timeseries(my.post)

# regularisation
my.post.regularise = HyperTraPS(m.2, initialstates = m.1, regularise = 1, walkers = 20)
plotHypercube.regularisation(my.post.regularise)
plotHypercube.summary(my.post.regularise, continuous.time = FALSE)

# simulated annealing output
my.post.sa = HyperTraPS(m.2, initialstates = m.1, sa = 1)
plotHypercube.summary(my.post.sa, continuous.time = FALSE)

# phenotypic landscape inference
my.post.pli = HyperTraPS(m.2, initialstates = m.1, pli = 1)
plotHypercube.summary(my.post.pli, continuous.time = FALSE)

# start with every edge parameterised, then regularise
my.post.bigmodel.regularise = HyperTraPS(m.2, initialstates = m.1, model = -1, regularise = 1, walkers = 20)
plotHypercube.regularisation(my.post.bigmodel.regularise)
