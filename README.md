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


### GPU


ami=Deep Learning AMI with Source Code (CUDA 8, Amazon Linux)
region=Oregon
instance type=m3.xl
$0.27/hr
ssh -i private/datasci.pem ec2-user@34.211.111.203
cat /home/ec2-user/src/README.md
install cudamat:
    git clone https://github.com/cudamat/cudamat
    sudo /home/ec2-user/src/anaconda2/bin/python2.7 setup.py install
    PYTHONPATH=$PYTHONPATH:/home/ec2-user/cudamat
    
install gnumpy:
    /home/ec2-user/src/anaconda2/bin/pip install -i https://pypi.anaconda.org/pypi/simple gnumpy
ipython2


