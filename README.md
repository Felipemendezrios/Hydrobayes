# Hydrobayes
Repository for developing a methodology for the automatic calibration of 1D hydraulic models, initially focusing on spatially distribute the friction coefficient.

# To Do:

- Implement multi-reach case. Changes in readBIN function in MAGE model (Ok, already implemented but not tested yet)
- Test new error model depending on the output variable (Z,Q,V).
  - It could be estimated using a covariate as water depth, donc WSE minus bed river.
  - Be careful during estimation of loglikelihood with several type of measures. Remember P(Q) + P(H) if we are calculating the probabilities in logarithm scale.
  - For multi-reach case, it will be necessary to take matrix for creating the error models. Save the idea to use a binairy matrix to activate or desactivate the bief and the a vector with gamma parameters depending on the reach
- Finish writting the document about structural error model in MAGE model (almost finished)
- Implement a Python module for calculating the initial condition for each event automatically, instead of using PamHyr to get them (Not a good idea, it is preliminary treatment)

# Plan:
- Hydraulic flume case with spatial distribute friction, either in floodplain or main channel. (already tested)
- Implementation of multi-reach experiment (already implemented, not tested yet)
- Taste the methodology in a real case: Rhone river (A model between Anthon and Port Chazey must be created, in progress)
- Structural error spatialisation 

# To remember:

- In multi-event case, some files are mandatory and must be equal to all the events in order to run model: `.RUG`, `.NET`. Some optional files are advised to be specified by event: `.HYD`,`.LIM`,`.PAR`,`.INIT`,
- Time coordinates to specify in calibration data must be in seconds (Model requirement)
- Parameters of structural error for all output variable must be specified, event if any calibration data is available. That's let avoid to accept all values of the parameters, which will tend to infinite value. Same for floodplain friction parameter.
- Assign absolute values to the structural error model (gamma)
- For hydraulic flume, reduce the minimal water surface elevation in `.PAR` file, defining a default value of 10 cm.
- The order of the reach in a Mage model are given in `.NET` file, starting from reach 1. Ensure that observations follow this order to put the observation at the good reach


# Notes to keep in mind after vacations:
- Calibration (+ plots): the most recent file is here: `/Case_studies/HHLab/Compound_channel_spatial_friction/Calibration_compound_channel_transition_friction.r`
- Prediction (+plots): the most recent file is here: `/Case_studies/HHLab/Rectangular_channel/Prediction_rectangular_channel.r`
- Functions folder : must be verified with local function in the scripts: `/Functions` 
- Rhone case study: `/Case_studies/Rhone/latest_advances/`. Not finished yet, model needs to be completed with the Ain tributary