# HEALPix installation
### 1. Install Dependencies (`cfitsio`)
To read and write FITS files, HEALPix strictly requires the **CFITSIO** library. Without it, the C++ compilation will fail or lack FITS capabilities.
* **Ubuntu/Debian**: `sudo apt-get install libcfitsio-dev`
* **macOS**: `brew install cfitsio`
* Ensure you also have standard build tools (`g++` or `clang++`, and `make`) - on Debian they are obtainable at `sudo apt install build-tools`

### 2. Download and Extract HEALPix
Download the latest version of HEALPix from [SourceForge](https://sourceforge.net/projects/healpix/).
Once downloaded, extract it and navigate into the folder:

    tar -xzpf Healpix_3.83.tar.gz
    cd Healpix_3.83

### 3. Compilation
Compile via `./configure`. It is helpful to use the automated flag, in the case of C++ :

    ./configure --auto=cxx

Be aware that depending on the installation folder it could be necessary to link the directories of *cfitsio*.
On a standard linux install they can simply be passed with

    FITSDIR=/usr/lib FITSINC=/usr/include ./configure --auto=cxx

Alternatively your installation directory can be found via `find /usr -name "fitsio.h"`

### 4. Make and test
Now compile the project with `make` and test the installation with `make test`. Afterwards the test files can be cleaned up with `make tidy`.