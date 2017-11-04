## spiral-waves

by Joe Hahn,<br />
jmh.datasciences@gmail.com,<br />
3 August 2017<br />
git branch=master

###Intro:

blah blah blah...

### Technical Notes:

1 Install current conda v4.3.21:

    wget https://repo.continuum.io/miniconda/Miniconda2-4.3.21-MacOSX-x86_64.sh
    chmod +x ./Miniconda2-4.3.21-MacOSX-x86_64.sh
    rm -rf /Users/joe/miniconda2
    ./Miniconda2-4.3.21-MacOSX-x86_64.sh -b -p /Users/joe/miniconda2
    rm Miniconda2-4.3.21-MacOSX-x86_64.sh

add this to ~/.bash_profile

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


    
### Intro:

