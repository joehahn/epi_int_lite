## epi_int_lite

by Joe Hahn,<br />
jmh.datasciences@gmail.com,<br />
10 January 2018<br />
git branch=master

###Intro:

This is epi_int_lite, which is the successor of the epi_int code that I had 
developed to simulate the dynamical evolution of gravitating rings, with this version
of epi_int using all open source libraries.

### Installation:

1 Install Anaconda python v4.3.21:

    wget https://repo.continuum.io/miniconda/Miniconda2-4.3.21-MacOSX-x86_64.sh
    chmod +x ./Miniconda2-4.3.21-MacOSX-x86_64.sh
    rm -rf /Users/joe/miniconda2
    ./Miniconda2-4.3.21-MacOSX-x86_64.sh -b -p /Users/joe/miniconda2
    rm Miniconda2-4.3.21-MacOSX-x86_64.sh

2 On my Mac laptop I edit ~/.bash_profile to let my PATH know where A

    export PATH="/Users/joe/miniconda2/bin:$PATH"
    echo $(conda --version)

and then install these libraries:

    conda install -y ipython
    conda install -y scipy
    conda install -y pandas
    conda install -y matplotlib
    conda install -y seaborn
    conda install -y jupyter
    conda install -y jupyter_dashboards -c conda-forge

3 Clone this repo:

    git clone https://github.com/joehahn/epi_int_lite.git
    
### Intro:


