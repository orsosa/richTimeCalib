# RICH Time Calibration package
This package contains a set of tools to performs the time calibration of the CLAS12 RICH detector.

## Requirement
* Hipo4 c++ installed and environments variables ```HIPO4LIB``` and  ```HIPO4INC``` set to the hipo4 library and headers respectively
* Root version >= 5.34

## Getting ready
```
git clone https://github.com/orsosa/richTimeCalib
cd richTimeCalib
make
make install
setenv RICHTIMECALIB `pwd`
```
add the enviroment variable ```RICHTIMECALIB``` to your login script.
## Use
```
cd $RICHTIMECALIB/script
./send_richTimeCalib.py <path-to-calib-files or list-of-calib-files> <run-number>
```
## Structure
```
├── include
│   ├── Clas12Banks4.h
│   ├── Clas12Detectors.h
│   ├── Clas12Pid.h
│   ├── process_opt.h
│   ├── RichGeometry.h
│   ├── RichHW.h
│   ├── RichTimeCorr.h
│   └── RichTT.h
├── macros
│   ├── richTimeOffsets.C
│   └── richTimeWalks.C
├── Makefile
├── maps
│   ├── SspRich_mapCHANNEL2PIXEL.txt
│   └── SspRich_mapFIBER2PMT_sortbyPMT.txt
├── README.md
├── script
│   ├── richTimeCheck.C
│   ├── run_richTiming.sh
│   ├── run_TC.sh
│   ├── run_TO_corr.sh
│   ├── run_TW_corr.sh
│   └── send_richTimeCalib.py
└── src
    ├── richTimeCheck.cxx
    ├── richTimeOffsets.cxx
    ├── richTimeWalks.cxx
    ├── richTiming_app.cxx
    ├── richTiming.cxx
    └── richTiming_expert.cxx
```
After compilation two folders are created ```bin/``` and ```.obj/```
## FAQ
