# SBT Cell SImulation

## Introduction
A Geant4 simulation of liquid scintillator cells being developed for the Surround Background Tagger of the SHiP experiment.

### Branches

<dl>
  <dt><code>master</code></dt>
  <dd>Main development branch, compatible with Geant4 v11.x. Contains most up to date parameters of 1 SBT cell (Aluminum with 5 mm wall thickness, matching the dimensions of the aluminum cells tested in November 2024 at the CERN PS testbeam. Reflectivity of walls can be varied, no coating layer).</dd>
  <dt><code>alumicern24v11</code></dt>
  <dd>Geometry of the aluminum cells tested in November 2024 at the CERN PS testbeam with 2 cm air gap along the left wall of the cell (looking at WOMs). Reflectivity of walls can be varied, no coating layer.</dd>
  <dt><code>steelcern24v11</code></dt>
  <dd>Geometry of the stainless steel cell tested in November 2024 at the CERN PS testbeam with 2 cm air gap along the left wall of the cell (looking at WOMs). Reflectivity of walls can be varied, no coating layer.</dd>
  <dt><code>reflectionhists</code></dt>
  <dd>Same as master (as of October 2025) with added histograms that output the number and type of reflections for photons which reach the WLS layer. Can be used to study behaviour of specular vs diffuse reflections.</dd>
  <dt><code>4cellsv11</code></dt>
  <dd>Geometry of 4 cell prototype tested in March/April 2024 at the CERN PS testbeam (10 mm steel walls). Contains BaSO4 coating layer with reflectivity according to data sheet, which can be scaled.</dd>
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

    Not recommended: local
    ```bash
    source /usr/local/[...]/geant4make.sh
    ```
3. Move to code directory
    ```bash
    cd ship_wom_testbeam
    ```

4. Create build directory
    ```bash
    mkdir build
    ```
    
5. Move to build directory
    ```bash
    cd build
    ```
6. Run cmake
    ```bash
    cmake ..
    ```
7. Run make
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

This code is compatible with Geant4 release 4.11.2. Not compatible with Geant4 version 4.10.x (unless using old branch).
