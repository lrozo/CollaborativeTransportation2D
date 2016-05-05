# CollaborativeTransportation2D
The MATLAB code shows a simple collaborative transportation task carried out on a planar workspace. The models, algorithms and results given in this code are linked to the T-RO paper "Learning Physical Robot Collaborative Behaviors from Human Demonstrations"

### Compatibility

	- DEPENDENCIES:
		* CVX   Matlab toolbox for convex optimization.(http://cvxr.com/cvx/)

	The code has been tested with Matlab. Note that the CVX toolbox is currently not working with GNU Octave.
	However, this will be compatible with GNU Octave soon (http://cvxr.com/news/2014/03/cvx-on-octave/). 

### Usage

	The user can control the execution of different	parts of the code by setting the flags for learning, 
	BIC-based model selection and stiffness estimation. 
	
	- Flags
		needsModel:	Used for learning or loading a TP-GMM 
		BICselection:  	Used for carrying out a BIC-based model selection
		stiffEstimate:	USed for estimating/loading stiffness matrices (You need to have the CVX toolbox)

### Reference  
	
	L. Rozo, S. Calinon, D. Caldwell, P. Jimenez, and C. Torras. Learning Physical Collaborative Robot 
	Behaviors from Human Demonstrations. Transactions on Robotics. 2016

### Description

	This code shows a simple human-robot cooperative transportation task in a planar scenario. The task consists 
	of transporting an object from an initial position to a target location. This code shows:
 
	1. Computation of virtual attractor trajectories from the dynamics observed during the demonstrations of the 
	   collaborative task.
	2. TP-GMM learning of the collaborative task.
	3. BIC-based model selection.
	4. Stiffness estimation built on a convex optimization.
	5. Reproduction using GMR with adaption to new configurations of the task parameters.

	* TASK PARAMETERS
	Three different task parameters are defined in this experiment, namely: 
		- Initial position and orientation of the object
		- Target position and orientation of the object
		- Irrelevant frame that randonmly varies its position and orientation

### Authors

	Leonel Rozo and Sylvain Calinon
	http://programming-by-demonstration.org/
		
	This source code is given for free! In exchange, we would be grateful if you cite the following reference in any academic publication that uses this code or part of it:

	@article{RozoTRO2016,
		author	="Rozo, L. and Calinon, S. and Caldwell, D. and Jimenez, P. and Torras, C.",
		title	="Learning Physical Collaborative Robot Behaviors from Human Demonstrations",
		journal	="IEEE Trans. on Robotics",
		year	="2016"
	}

