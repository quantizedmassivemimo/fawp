# Simulator for Finite-Alphabet Wiener Filter Precoding (FAWP)
(c) 2019 Oscar Castañeda and Christoph Studer
e-mail: caoscar@ethz.ch & studer@ethz.ch

More information about our research can be found at [https://iip.ethz.ch] and [https://sites.google.com/site/durisi].

### Important information

If you are using the simulator (or parts of it) for a publication, then you *must* cite our paper:

Oscar Castañeda, Sven Jacobsson, Giuseppe Durisi, Tom Goldstein, and Christoph Studer, "Finite-Alphabet Wiener Filter Precoding for mmWave Massive MU-MIMO Systems," Asilomar Conference on Signals, Systems, and Computers, November 2019, pp. 178-183

and clearly mention this in your paper.

### How to start a simulation:

Simply run

```sh
fa_precoder_sim
```

which starts a simulation of a 256 BS antenna, 16 user, 16-QAM massive MIMO system using the four finite-alphabet Wiener filter 
(FAWP) precoding methods proposed in our paper, Pre-FAWP-WF, Post-FAWP-WF, Pre-FAWP-FBS and Post-FAWP-FBS, all using the 1-bit alphabet. As a baseline, the simulator also includes the conventional Wiener filter (WF) precoder, whose entries are not quantized. The simulation runs in an i.i.d. Rayleigh-fading channel.

The simulator runs with predefined parameters. You can specify your own system and simulation parameters by passing your own "par"-structure (see the simulator for an example). Note that we use default parameters for the considered system configuration; if you want to run the simulation with a different system configuration, you will need to tune the parameters of the Pre-FAWP-FBS and Post-FAWP-FBS algorithms to obtain the best performance.

We highly recommend you to execute the code step-by-step (using MATLAB's debug mode) in order to get a detailed understanding of the simulator.

### Version history
* Version 0.1: ofcastaneda - initial version for GitHub release
