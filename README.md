# Dynamic Mode Decomposition (DMD) for Simulator fOr Wind Farm Applications (SWOFA)

This project is inserted in the Master Thesis Dissertation entitled "Economic Model Predictive Control: an Extended Dynamic Mode Decomposition Approach", where Dynamic Mode Decomposition is used with subspace identification methods and the ideias behind the Koopman Operator to derive a tractable model to perform Economic Model Predictive Control.
The Master Thesis is currently being developed, so weekly commitments are expected.

The higher level goal of this code is to perform Dynamic Mode Decomposition in any set of data retrieved from situations in SOWFA, thus allowing one to derive models for wind turbine control. As simulations might differ in terms of control strategy (pitch control, yaw control), the user should thus take the following into consideration: (all changes should only be necessary - ideally - in the main script - MPC_wakesteering):



### Control Strategy

#### BurgersExample 
Runs the Burgers example as explained in the paper, it includes data collection, Extended Dynamic Mode Decomposition (EDMD) for identification of the Koopman linear system, and a run of closed-loop controlled system from some initial condition.
Feel free to play with the paremeters of the code, specially, try different observables, embedding dimension, reference signal, initial condition, etc.
The whole program, with the initial paremeter settings, runs on my personal laptop in under 2 minutes.


#### CavityExample
Runs the lid-driven cavity flow example as explained in the paper,  including  EDMD for identification of the Koopman linear system, and a run of closed-loop controlled system from some initial condition on the limit cycle. There are two options to run this code:
1- ask the code to generate data for EDMD. This is a lengthy process and for the parameter values reported in the paper takes ~10 hours on a powerful desktop (with no parallelization), or 2- go to https://ucsb.box.com/s/367tvkgnzby61x9nrh64q81748ugaw63 and download the data file "Cavity_data_4EDMD_0" (~3GB) which is the data used in the paper. Using the data file, the program  takes about 5 minutes to run on my laptop. 




### before you run the code:

go to "./thehood/" and unzip "qpOASES-3.1.0",
then go to subfolder ".\thehood\qpOASES-3.1.0\interfaces\matlab" and run make.m .
This is required to activate the qpOASIS interface for solving the optimization problem.


send comments and questions to
#### arbabiha@gmail.com

H Arbabi

April 2018
