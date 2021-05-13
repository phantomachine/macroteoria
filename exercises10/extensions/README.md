Unemployment over the business cycle and DMP with shocks
=========================================================

Empirical Case Study
--------------------
Here we replicate some famous dataset describing the phenomenon of job unemployment, vacancies, job finding rate, market tightness over a sample of the US Business cycle.

File: 

* ``Unemployment Regularities.ipynb``



DMP's Causal Mechanism and Quantitative Theory
-----------------------------------------------
Then we will work towards synthesizing the empirical regularities using the Nobel-prize winning model of Diamond, Mortensen and Pissarides. First, we'll start simple, by studying the DMP model's deterministic behavior using a local, first-order perturbation solution method. Then, we'll come back full-circle and employ a global, nonlinear solution method called *time interation* (TI) which is similar to what we studied earlier by successively getting better approximations of a nonlinear decision function characterized by a system of equilibrium/Euler conditions.

A local (deterministic) version of DMP to illustrate basic dynamics of DMP
--------------------------------------------------------------------------

This first part replicates the discussion in Miao's textbook 18.1.2 (transitional dynamics). We'll study the deterministic, dynamic behavior of the DMP model.

File: 

* `Mortensen-Pissarides S&M Model.ipynb`

A Non-linear global solution of DMP
-------------------------------------

Now we do the second part, using a TI method. Petrosky-Nadeau and Zhang demonstrated that in the DMP model the old-school practice of log-linearizing the model equilibrium mapping to get local approximant/solutions may be inaccurate. (Unlike in typical neoclassical and NK models, models with search market frictions may generate considerably greater nonlinearity due to externalities in decentralized market trades.)

This example computes a version of the Diamond-Mortensen-Pissarides labor search and matching model.
This version replicates Petrosky-Nadeau and Zhang's QE example of Hagedorn and Manovskii.

There is an ocassionally binding non-negativity constraint on vacancies :math:`V_{t}`.

Expect that this constraint tends to bind in states with low employment and low productivity.

Here we illustrate a time-iteration solution method.

This finds the fixed point (policy function) satisfying the system of nonlinear FOCs.

We solve this example using a version of a finite-element method on local hierarchical sparse grids.

Files:

* ``dmp.py`` is the class file containing model equilibrium setup and tools

* ``Business-cycle Search and Matching.ipynb`` Example Jupyter Notebook

Things to do:

* Alternative version with Markov chain shocks instead of AR(1)

* Speedups comparisons using alternatives:

	* NUMBA
	
	* OpenMPI

Dependencies:

* TASMANIAN. See:

	* [website](https://tasmanian.ornl.gov/) 
	* [website for Python interface](https://pypi.org/project/Tasmanian/)
	* Install using ``pip install Tasmanian --user``

* STATSMODEL. See [website](https://www.statsmodels.org/)

(c) 2020++ T. Kam (tcy.kam@gmail.com)
