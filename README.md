# fm128_radar_knmi
Tool to create WRFDA fm128_radar ascii based files from KNMI data

## Installation

fm128_radar_knmi is installable via pip:
```
pip install git+https://github.com/ERA-URBAN/fm128_radar_knmi
```
fm128_radar depends on the following packages:
```
fm128-radar
asteval
ConfigArgParse
dateutils
h5py
netCDF4
numpy
python-dateutil
```

## Usage
```
usage: fm128_radar_knmi [-h] [-c MY_CONFIG] --DOWNLOAD_DIR DOWNLOAD_DIR
                        --OUTPUT_DIR OUTPUT_DIR [--TMP_DIR TMP_DIR]
                        [--remove_intermediate REMOVE_INTERMEDIATE]
                        datetime

Download and convert KNMI radar data to WRFDA fm128_radar ascii Args that
start with '--' (eg. --DOWNLOAD_DIR) can also be set in a config file
($INSTALL_DIR/etc/fm128_radar_knmi/radar.config or specified via -c).
Config file syntax allows: key=value, flag=true, stuff=[a,b,c] (for details,
see syntax at https://goo.gl/R74nmi). If an arg is specified in more than one
place, then commandline values override environment variables which override
config file values which override defaults.

positional arguments:
  datetime              Datetime string "%Y%m%d%H%M"

optional arguments:
  -h, --help            show this help message and exit
  -c MY_CONFIG, --my-config MY_CONFIG
                        config file path
  --DOWNLOAD_DIR DOWNLOAD_DIR
                        Download directory tar files [env var: DOWNLOAD_DIR]
  --OUTPUT_DIR OUTPUT_DIR
                        Output directory of ob.radar file [env var: OUTPUT_DIR]
  --TMP_DIR TMP_DIR     Directory where intermediate files are saved, defaults
                        to DOWNLOAD_DIR [env var: TMP_DIR]
  --remove_intermediate REMOVE_INTERMEDIATE
                        Remove intermediate files

```
