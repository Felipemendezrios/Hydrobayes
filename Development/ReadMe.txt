Folders developement:
- New_feautures_Ben and New_feautures_Ben_spatial are tests with Polynomial degree 0 and 1, respectively.
There are already tested and moved to Verification folder. Now, they are using to understand why parallel computation are not possible, problems with readBin function.
- New_features_Ben_RBaM: Start a script to integrate RBaM for using new model MAGE (MAGE_ZQV)

- run_model_BaM: it is the folder which we are working on.
This file contains events to test MAGE_ZQV model in BaM, but outside the BaM environment. In other words, I get only access to run_model function to run Mage Model.
The likelihood is calculated by myself, using some RBaM functions, but I can modify all that I want.
- New structural error models: spatial regression? different by output variable (sure)?
- I need to integrate new function to get prediction phases
- Need to pass error model parameters for propaging total uncertainty case.

- run_model_multi_reach = file to run a test for multi reach implementation saved in verification folder