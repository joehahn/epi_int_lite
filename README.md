## epi_int_lite

by Joe Hahn,<br />
jmh.datasciences@gmail.com,<br />
10 January 2018<br />
git branch=lindblad_resonance


### Intro:

This is the epi_int_lite N-body integrator (which is the successor to my earlier epi_int N-body
integrator) that is used simulate the dynamical evolution of gravitating rings.


### Installation:


1 Install Anaconda python. Browse this url and find the Anaconda python installer that is appropriate to your OS:

    https://repo.anaconda.com/archive

I'm using a MacOSX laptop , and I prefer to code in python 2.7 (and NOT 3.X), with v5.3.0 being
the latest version when I first started developing epi_int_lite, and the following
uses wget to download Anaconda2-5.3.0-MacOSX-x86_64.sh to my laptop:

    wget https://repo.anaconda.com/archive/Anaconda2-5.3.0-MacOSX-x86_64.sh
    chmod +x ./Anaconda2-5.3.0-MacOSX-x86_64.sh

2 the following installs v5.3.0 Anaconda python into the anaconda2 folder in my top directory:

    #rm -rf ~/anaconda2
    ./Anaconda2-5.3.0-MacOSX-x86_64.sh -p ~/anaconda2 -b
    rm Miniconda2-4.3.21-MacOSX-x86_64.sh

3 Set the PYTHON_PATH bash variable to point at python, and check your work:

    PYTHON_PATH=~/anaconda2/bin
    $PYTHON_PATH/python --version

4 enable jupyter_dashboards:

    $PYTHON_PATH/conda install -c conda-forge -y jupyter_dashboards
    $PYTHON_PATH/jupyter nbextension enable jupyter_dashboards --py --sys-prefix
    
5 Clone this repo:

    git clone https://github.com/joehahn/epi_int_lite.git
    cd epi_int_lite


### Run test simulations:

The tests folder contains additional folders, each of which test some aspect of epi_int_lite.
The test detailed here uses epi_int_lite to monitor the libration of six self-gravitating ringlets
that differ in total mass. To simulate one of those 6 ringlets:

    cd tests/libration/ring_mass_1.5e-10
    $PYTHON_PATH/python ./epi_int_lite.py

Then navigate to another folder and simulate a somewhat higher-mass ringlet:

    cd ../ring_mass_5.0e-10
    ./epi_int_lite.py

Then repeat the above for the remaining ring_mass_* folders, this generates all the output
the is needed by the Jupyter notebook that will compare simulation output to theory. To see that
comparison, navigate to the tests/libration folder and start Jupyter:

    cd tests/libration
    jupyter notebook
    
and click libration_frequency.ipynb to load that notebook and then click >> to refresh the output.
That notebook loads the output generated by the six epi_int_lite simulations, computes each
ringlet's libration period, and then compares the simulated ringlets' libration periods
to the predicted by Borderies Goldriech & Tremaine (BGT). Note that the simulated periods
agree with the BGT for higher-mass ringlets, as expected, but that simulated results diverge
from theory for low-mass rings, but that is OK too since a close examination of the BGT
paper will show that their assumptions are violated by the low-mass ringlets.

Many of the test simulations include an animation script called animate.py,
and execute the following to view those animations:

    cd tests/libration/ring_mass_1.5e-8
    ./animate.py

Also check out the remaining tests, their notebooks, and their animations.

When attempting to run those tests, keep in mind that epi_int_lite is still under development
and is evolving. So if any of my code tweaks alter epi_int_lite's expected input or output, then either
./epi_int_lite.py will fail, or a notebook will fail to read its output. This just means that the
particular test that you are looking at is out of sync with the source, and that that particular
test's input.py or notebook needs to be tweaked. If so, send me an email
and I'll eventually get that test debugged and this repo updated.


### To do:

1 finish this README

2 bump this 2D model up to 3D

3 do a buncha scientifically relevant runs

4 write frickn paper

