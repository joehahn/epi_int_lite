## epi_int_lite

by Joe Hahn,<br />
jmh.datasciences@gmail.com,<br />
10 January 2018<br />
git branch=master


### Intro:

This is the epi_int_lite N-body integrator (which is the successor to my earlier epi_int N-body
integrator) that is used simulate the dynamical evolution of gravitating rings.


### Installation:

1 Clone this repo:

    git clone https://github.com/joehahn/epi_int_lite.git
    cd epi_int_lite
    
2 install Anaconda python v4.3.21:

    wget https://repo.continuum.io/miniconda/Miniconda2-4.3.21-MacOSX-x86_64.sh
    chmod +x ./Miniconda2-4.3.21-MacOSX-x86_64.sh
    rm -rf /Users/joe/miniconda2
    ./Miniconda2-4.3.21-MacOSX-x86_64.sh -b -p /Users/joe/miniconda2
    rm Miniconda2-4.3.21-MacOSX-x86_64.sh

3 On my Mac laptop I edit ~/.bash_profile to let my PATH know where Anacoda python is

    export PATH="/Users/joe/miniconda2/bin:$PATH"
    echo $(conda --version)

4 install additional python libraries:

    conda install -y ipython
    #conda install -y scipy
    conda install -y pandas
    conda install -y matplotlib
    conda install -y seaborn
    conda install -y jupyter
    #conda install -y jupyter_dashboards -c conda-forge


### Run test simulations:

1 the tests folder contains additional folders, each of which test some aspect of epi_int_lite,
and the following tests the ring particles unperturbed motion by confirming that particle orbits
precess at the expected rate when orbiting an oblate planet. To run that test simulation:

    cd tests/rates
    ./epi_int_lite.py

then start jupyter notebook:

    jupyter notebook
    
and navigate to check_rates.ipynb and click Kernel > Run All to refresh the output. Consult
notebook comments for further details.


### To do:

1 finish this README

2 fix libration test

3 use interpolation to compute nearest site on perturbing streamline

4 bump this 2D model up to 3D

5 use an anaconda environment to version-control the python libraries used here


