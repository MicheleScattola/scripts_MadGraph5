# Home use of MadGraph5 and Delphes on local linux machine
This guide provides basic steps to install and run madgraph, Pythia8 and Delphes on a local machine.

Note that all commands are written for Debian Linux distributions (Ubuntu) but can be very easily adapted to all linux versions
(e.g. if your package manager is not apt...)

# 1. UV virtual environment

References for UV manuals are found at [UV docs](https://docs.astral.sh/uv/getting-started/)

### Install

Install uv via the official standalone installer:

    curl -LsSf https://astral.sh/uv/install.sh | sh

check status with:

    uv

### Initialization

either create your project with uv (this creates the directory for you)

    uv init project-name

or create the directory and initialize

    mkdir project-name
    cd project-name
    uv init

### Environment and dependencies

Install latest python version 

    uv python install

or choose your version (e.g. 3.12)

    uv python install 3.12

you can check the installed python version in the environment with

    cat .python-version

Now add your necessary project dependencies:
(this creates automatically the virtual environment in .venv/bin/activate )

    uv add six numpy requests matplotlib

Python scripts can then be run via

    uv run file.py

Note that with this syntax there is no need to source the environment, uv does it on its own.
If you do want to source it manually before working, use

    source .venv/bin/activate

# 2. MadGraph5

### Install

**Navigate to your desired directory** and download the latest stable version

    wget https://launchpad.net/mg5amcnlo/3.0/3.6.x/+download/MG5_aMC_v3.6.7.tar.gz

Extract the folder in your current directory

    tar -xzf MG5_aMC_v3.6.7.tar.gz


### Dependencies (Debian)
Make sure you have the necessary requirements for compiling C++ and Fortran. Note that this section is specific for Debian distributions, however, the commands can easily be changed to suit your own Linux distribution with a quick google search.

Check if you already have the compilers installed:

    which gfortran
    which g++
    
Instructions for gfortran are found at [fortran-lang.org](https://fortran-lang.org/learn/os_setup/install_gfortran/). **For Debian based distributions**:

    sudo apt install gfortran
    sudo apt install build-essential

Note that build-essential includes other packages (gcc, make, etc.) which are often useful/required. If you strictly want g++ you can install via

    sudo apt install g++

### Ready to use
Now you are ready to enter MadGraph5 shell by navigating

    cd MG5_aMC_v3_6_7

and launch the shell

    ./bin/mg5_aMC

or either write your scripts and launch them

    ./bin/mg5_aMC script.txt

### Install hepmc and Pythia

Enter the shell and download the necessary tools:

    ./bin/mg5_aMC
    install hepmc
    install hepmc3
    install pythia8

# 3. ROOT and Delphes
Delphes requires root to be locally installed onto your machine.
Different linux distributions need different pre-compiled versions, but any information can be found over at: [ROOT CERN install](https://root.cern/install/#download-a-pre-compiled-binary-distribution)

### Dependencies
Some dependencies are needed (see [ROOT dependencies](https://root.cern/install/dependencies/)). For Ubuntu:

    sudo apt install binutils cmake dpkg-dev g++ gcc libssl-dev git libx11-dev \
    libxext-dev libxft-dev libxpm-dev python3 libtbb-dev libvdt-dev libgif-dev

### Download and extract
Navigate to the desired folder of installation for root, then download the latest pre-compiled binary. The latest version is shown at [ROOT latest](https://root.cern/install/all_releases/); choose your version and download. In March 2026 for Ubuntu:

    wget https://root.cern/download/root_v6.36.10.Linux-ubuntu24.04-x86_64-gcc13.3.tar.gz

and then extract the folder:

    tar -xzvf root_v6.36.10.Linux-ubuntu24.04-x86_64-gcc13.3.tar.gz

Now root is install in the folder ``current-path/root``. You can make it accessible in your current enviroment via:

    source root/bin/thisroot.sh


### Add to bashrc
To make ROOT accessible to all users automatically add the following line to the ``~/.bashrc`` file. Open the file with your favorite editor (code, nano, vim, etc...)

    code ~/.bashrc

scroll to the bottom of the file and add the full path to where you installed root:

    source /path/to/root/bin/thisroot.sh

now exit the terminal and open a new session. You should now be able to access root anywhere by just typing

    root

### Install Delphes
Navigat to your MadGraph installation folder   

    cd /path/to/MG5_aMC_v3_6_7

and launch the shell

    ./bin/mg5_aMC

Then install Delphes

    install delphes

