## epi_int_lite

by Joe Hahn,<br />
jmh.datasciences@gmail.com,<br />
10 January 2018<br />
git branch=master


### Intro:

This is the epi_int_lite N-body integrator (which is the successor to the epi_int N-body code)
that is used simulate the dynamical evolution of gravitating rings.


### Installation:

1 Install Anaconda python v4.3.21:

    wget https://repo.continuum.io/miniconda/Miniconda2-4.3.21-MacOSX-x86_64.sh
    chmod +x ./Miniconda2-4.3.21-MacOSX-x86_64.sh
    rm -rf /Users/joe/miniconda2
    ./Miniconda2-4.3.21-MacOSX-x86_64.sh -b -p /Users/joe/miniconda2
    rm Miniconda2-4.3.21-MacOSX-x86_64.sh

2 On my Mac laptop I edit ~/.bash_profile to let my PATH know where Anacoda python is

    export PATH="/Users/joe/miniconda2/bin:$PATH"
    echo $(conda --version)

3 Install these libraries:

    conda install -y ipython
    #conda install -y scipy
    conda install -y pandas
    conda install -y matplotlib
    conda install -y seaborn
    conda install -y jupyter
    #conda install -y jupyter_dashboards -c conda-forge

4 Clone this repo:

    git clone https://github.com/joehahn/epi_int_lite.git


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


