# Ising2D

### Summary
This code, consisting of the macros Metropolis.C and Wolff.C, was created to simulate the 2D Ising model. "Metropolis.C" is an implementation of the familiar Metropolis algorithm, while "Wolff.C" uses the Wolff cluster-flip algorithm:
- Pick a random lattice site. Remember the value of the spin at that site, and add it to the "cluster".
- Check all adjacent sites. If they have the same spin as the cluster, add them to the cluster with probability p = 1 - exp(-2J/(kT)). (Here, J = 1 is the coupling strength, k = 1 is the Boltzmann constant)
- Repeat the previous step recursively for each spin that has just been added to the cluster.
- Flip the whole cluster.

### Usage
The two macros are littered with a number of functions. For Metropolis.C, the important ones are:
```sh
void init(UInt_t size = initial_size) // initialize the lattice, size is the side-length of the (square) lattice
void setT(Double_t T = Tc) // set the temperature, default is critical temperature
void setB(Double_t B = 0 ) // set the magnetic field, default is zero
void run() // runs the animation (to be used after init() has been called)
void stop() // stops the animation
void rg(UInt_t power = 1) // perform real-space renormalization by blocking spins into clusters of 4, using "majority rules". Perform it "power" # of times.
```
Wolff.C has the above functions. In addition:
```sh
void stepMany(TString id = "0") // run for size^2 steps, save the resulting configuration in a TCanvas. The TString "id" is used in the root file name
```
### Dependencies
[![N|Solid](https://d35c7d8c.web.cern.ch/sites/d35c7d8c.web.cern.ch/files/website-banner-allnew-croped_3.png)](https://root.cern.ch)

These macros makes use of ROOT, the data analysis framework developed at CERN (see https://root.cern.ch for information on installing ROOT).
The code was run in the ROOT interpreter. For example, one may run the following:
```sh
$ root -l
$ root [0] .L Wolff.C+
$ root [1] init(100)
$ root [2] setT(0.99 * Tc)
$ root [3] run()
$ root [4] rg()
$ root [5] stop()
```

Jan Offermann
