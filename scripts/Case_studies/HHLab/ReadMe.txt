Here is an explanation of each experiment to address some adjacent questions and to help identify the right experiment for understanding the process.

** synthetic_cases: **

Observations are obtained from MAGE simulation, then a sample of simulation are used as observation to introduce in the calibration to estimate the initial scenario.
This case help to answer some questions from the number of observation needed to identify the parameters.

Synthetic cases are performed using a compound rectangular channel with constant longitudinal friction but different in the main channel and the floodplain.

Three cases were performed to answer these questions:

Question:
Can we identify the friction in the main channel and the floodplain using a single event with WSE observations from the floodplain?

Experiment: 1_WSE_floodplain

Question:
Can we identify the friction in the main channel using two events with WSE observations from the floodplain?

Experiment: 2_WSE_floodplain

Question:
Could two events with WSE observations in the main channel and the floodplain lead a better identifiability on the friction in tha main channel?

Experiment: 2_WSE_floodplain_and_main_channel





** Compound_channel_single_friction: **
Papers:
https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.57745/HJKRYH
DOI: 10.1017/jfm.2019.973

Compound channel case with a constant longitudinal friction, but different between the main channel and the floodplain.
The same mage projet is used for modelling the channel as in the synthetic case.

Information:
A single event with WSE observations in uniform flow condition.
Downstream boundary condition is assumed as the last observed WSE.

Question :
- The method could identify the friction in the main channel and the floodplain using a single event in floodplain?
- How to reduce the corelation between kmin/kmoy ?
- Non informative prior information force to get a minimal number of events (2) to identify the friction in the main channel and floodplain? Should one of the events flow only in the main channel?

Notes:
- High corelation between kmin/kmoy. They can be interchangeable and the result will be the same with a weighting or factor to compensate.
- Observation uncertainties were considered closer to zero to get a better identifiability of the friction coefficients



 ** Compound_channel_spatial_friction: **
 Papers:
https://doi.org/10.1007/s10652-017-9525-0
https://theses.fr/2016LYSE1154

 Compound channel case with a constant longitudinal friction in the main channel. In the floodplain, the longitudinal friction changes from wood to grass and vice versa.

Information:
Total discharge is estimated to 162 L/s, but in the real experiments, different discharges have been introduced in the main channel and floodplain.
Downstream boundary condition is assumed as the last observed WSE.

Two types of transition were measured:
-  Transition from wood to grass called 'CWMQ'. Upstream discharge in floodplain vary in 12 and 18 L/s. Measurements were taken either in the main channel or the floodplain.
-  Transition from grass to wood called 'CMWQ'. Upstream discharge in floodplain vary in 18 and 26 L/s. Measurements were taken only in the main channel.

Questions:
- What parameters are identifiables increasing the polynomial degree (n=0,...,n=5?) only in the floodplain?
- What is the behavior of a piecewise function to this case of study? -> break known both cases
- What is the influence of the downstream boundary conditions? -> compare case CMWQ and CWMQ

Notes:
- A high observation uncertainty will led to crash the MCMC sampling.
Because high value of Strickler coefficient could be estimated (low hydralic resistance), lead to estimate the discharge to gravitaional force and advection. Friction are not influencing the calculations
- It seems that Strikler coefficient is trying to represent the water backwash wave, which is controlling the WSE.
- Observation uncertainties were considered closer to zero to get a better identifiability of the friction coefficients (avoid compensation)


