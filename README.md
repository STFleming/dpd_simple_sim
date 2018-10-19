# Simple Dissipative Particle Dynamics (DPD) Simulator
#### 0.0.1

<p align="center">
    <img src="dpd_example.gif"/>
</p>

This is a simple DPD simulator based on the description in chapter 8 of [this article](http://www.cse.scitech.ac.uk/ccg/software/DL_MESO/MANUAL/USRMAN.pdf). 
It contains a C++ API for describing, beads (particles), the toroidal simulation universe, and the pairwise interaction forces between beads. It also has a javascript web-based interface (see above) that can be used to view the simulation in real-time, replay the simulation faster from the start, or to playback previous simulations.

### examples

To run one of the example in the `examples/` directory navigate to the appropriate directory and type:
```bash
make
make run 
```

Then using a web-browser on your local machine navigate to `http://localhost:3000`, where you should see a live rendering of the particles (similar the the gif at the top of this page).

### install requirements

```bash
apt-get install libjsoncpp-dev
apt-get install npm
apt-get install nodejs
```

### changelog

* 0.0.1 - initial release
