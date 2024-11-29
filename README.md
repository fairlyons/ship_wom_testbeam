# SBT Cell SImulation

## Introduction
A Geant4 simulation of liquid scintillator cells being developed for the Surround Background Tagger of the SHiP experiment.

### Branches

<dl>
  <dt><code>master</code></dt>
  <dd>Main development branch. Contains most up to date parameters of 1 SBT cell (Aluminum with 5 mm wall thickness).</dd>
  <dt><code>4cellsv11</code></dt>
  <dd>4 cell prototype geometry.</dd>
  <dt><code>other</code></dt>
  <dd>Old, including compatibility with Geant v10.7.</dd>
</dl>

## Build Instructions

1. Download the software
    ```bash
    git clone https://github.com/fairlyons/ship_wom_testbeam.git
    ```

2. Make sure to setup ROOT and Geant4

    Recommended: using `lxplus`

    ```bash
    source /cvmfs/sft.cern.ch/lcg/contrib/gcc/11.3.0/x86_64-el9-gcc11-opt/setup.sh
    source /cvmfs/geant4.cern.ch/geant4/11.2/x86_64-el9-gcc11-optdeb/CMake-setup.sh
    ```
    Local
    ```bash
    source /usr/local/[...]/geant4make.sh
    ```
4. Move to code directory
    ```bash
    cd ship_wom_testbeam
    ```

5. Create build directory
    ```bash
    mkdir build
    ```
    
6. Move to build directory
    ```bash
    cd build
    ```
7. Run cmake
    ```bash
    cmake ..
    ```
8. Run make
    ```bash
    make
    ```

## Run Instructions

If running locally: set up Geant4

```bash
source /usr/local/[...]/geant4make.sh
```

To run in interactive mode simply execute inside the build directory with no arguments

```bash
./OpNovice
```    

Now you can simulate some events, run with a macro file

```bash
./OpNovice -m run1.mac
```


## Dependencies

This code is compatible with Geant4 releases 4.11.2. Not compatible with Geant4 version 4.10.x
