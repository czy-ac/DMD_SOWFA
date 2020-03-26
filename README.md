# Dynamic Mode Decomposition (DMD) for Simulator fOr Wind Farm Applications (SWOFA)

This project is inserted in the Master Thesis Dissertation entitled "Economic Model Predictive Control: an Extended Dynamic Mode Decomposition Approach", where Dynamic Mode Decomposition is used with subspace identification methods and the ideias behind the Koopman Operator to derive a tractable model to perform Economic Model Predictive Control.
The Master Thesis is currently being developed, so weekly commitments are expected.

The higher level goal of this code is to perform Dynamic Mode Decomposition in any set of data retrieved from situations in SOWFA, thus allowing one to derive models for wind turbine control. As simulations might differ in terms of control strategy (pitch control, yaw control), the user should thus take the following into consideration: (all changes should only be necessary - ideally - in the main script - MPC_wakesteering):


### Data
The program requires 2 different types of data

(1) The

(2)

### Directories

### Control Strategy

### Pre Processing (preprocessdmd.m)

### Detrending States

### Dynamic Mode Decomposition 



Comments and questions:
#### nassir.cassamo@tecnico.ulisboa.pt or n.rodriguescassamo@tudelft.nl

Nassir Rodrigues Cassamo

MSc. Mechanical Engineering

TU Delft, Delft Center for Systems and Control (DCSC), Delft, The Netherlands

Instituto Superior Técnico, Lisbon, Portugal 


March 2020
