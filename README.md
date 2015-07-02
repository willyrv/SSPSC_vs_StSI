SSPSC vs StSI
=============

These python scripts were used to perform the simulations described in our 
paper *Demographic inference using genetic data from a single individual: 
separating population size variation from population structure* 
(see the original paper
  <a href=http://www.sciencedirect.com/science/article/pii/S0040580915000581>
  Mazet O., Rodríguez W., Chikhi L. (2015)</a> or a
  <a href=http://arxiv.org/abs/1412.1243> Preprint </a> for theoretical details).

The code may be used in order to reproduce the analyses performed in the paper
and in the Supplementary Materials. 

The code works fine under:
* Python 2.7

You will also need:

* numpy (version 1.9.2)
* scipy (version 0.15.1)

The file *experiment_settings.txt* allows to specify the parameters that you may
want to use for simulate the values.

The script *experiment.py* can be modified if you want to change the name of
the output file. It is recommended to run many experiments in parallel in order
to save time. To make an experiment just do 

*./experiments.py*

One single experiment produce two files. One contains the simulated T2 values
(i.e. the values of the coalescence times of two individuals under the actual
  model) and the other file contains the results of the experiment.

  Explanation of the output file
  ------------------------------

  The output file (named by default *experiment_OUT.txt*) contains 10 columns:

  * Original_variable: Depending on the model used to produce the T2 values
  this could be SSPSC (Single Step Population Size Change model) or StSI
  (Structured Symmetrical Island model). The model can be changed in the
  configuration file. If everything goes fine, the method
  should be able to identify the model used to produce data by doing the
  analysis of the simulated values of $T_2$.

  * real_parameters: The parameters of the model used for simulate data.

  * number_of_observations: The number of $T_2$ values simulated by the original
  model.

  * log-likelihood_of_real_params: The likelihood of the real parameters (the
    parameters of the model used for simulate the data), computed from the $T_2$
    values simulated by the model itself.

  * Estim_params_T2_SSPSC: The estimated parameters (by Maximum Lileklihood
    Estimation) assuming the $T_2$ values are coming from a SSPSC model.

  * log-likelihood_T2_SSPSC: The likelihood of the parameters estimated in the
  previous column computed from the $T_2$ values.

  * p-value_T2_SSPSC: Assuming that simulated $T_2$ values come from a SSPSC
  model with the parameters estimated in the previous step, we do a KS-test
  (*Kolmogorov-Smirnov* test) in order to see if the model fits the data.

  Now, the same process is repeated, but assuming data come from a StSI model

  * Estim_params_T2_StSI: The estimated parameters (by Maximum Lileklihood
    Estimation) assuming the $T_2$ values are coming from a StSI model.

  * log-likelihood_T2_StSI: The likelihood of the parameters estimated in the
  previous column computed from the $T_2$ values.

  * p-value_T2_StSI: Assuming that simulated $T_2$ values come from a StSI
  model with the parameters estimated in the previous step, we do a KS-test
  (*Kolmogorov-Smirnov* test) in order to see if the model fits the data.

  At this stage, we expect the KS test will reject the wrong model and will not
  reject the right one. For example, if data where simulated under the SSPSC
  model, the p-value_T2_StSI should be low (let's say, lower than 0.05) while
  the p-value_T2_SSPSC should not be too low (let's say, higher than 0.05).
  Given that in some cases, a KS test will not be enough to distinguish both
  models, an AIC (*Akaike Information Criterion*) approach may be used (you can
  find some explanations in
    <a href=https://en.wikipedia.org/wiki/Akaike_information_criterion>
    Wikipedia</a>). In our case, given that we are comparing two model that have
    the same number of parameters, using the AIC is equivalent to make a choise
    based just in the likelihood of estimated parameters under both models.  

  * AIC_selected_model: If it is 0 means that the selected model was the SSPSC
  and if it is 1, means that the selected model was the StSI.

  * AIC_relative_prob: Give some measure of how the chosen model is more
  likely to explain the data than the other.


You may do as many experiments as you want. Then you can process the results
with you favorite statistics software. I used some python scripts coded by
myself (./lib/results_handler.py). Sorry if they are not well documented, I
will explain how to use it soon.
