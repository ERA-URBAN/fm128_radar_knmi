'''
description:    Download KNMI radar hdf5 files (inside a tar archive) from ftp
license:        APACHE 2.0
author:         Ronald van Haren, NLeSC (r.vanharen@esciencecenter.nl)
                Natalie Theeuwes
'''
import datetime
import numpy
import h5py
import datetime
from netCDF4 import Dataset
from netCDF4 import date2num
import time
import asteval
import time


def initialize_nans(shape):
  '''
  Initialize an array with nans
  '''
  nans = numpy.empty(shape)
  nans.fill(numpy.nan)
  return nans


class convert_to_netcdf:
  def __init__(self, filename, dt, outfile='radar.nc'):
    self.filename = filename
    self.outfile = outfile
    self.dt = dt
    self.define_constants()
    self.calculate_reflectivity()
    self.create_netcdf()

  def define_constants(self):
    self.LAT_bilt = 52.10168
    self.LON_bilt = 5.17834
    self.r_earth = 6378137.
    self.angles = [0.3, 0.4, 0.8, 1.1, 2.0, 3.0, 4.5, 6.0, 8.0, 10.0, 12.0,
                   15.0, 20.0,  25.0 ]
    self.scans = ['scan1', 'scan2', 'scan3', 'scan4', 'scan5', 'scan6',
                  'scan7', 'scan8', 'scan9', 'scan10', 'scan11', 'scan12',
                  'scan13', 'scan14']

  def calculate_reflectivity(self):
    '''
    Calculate reflectivity for each point
    '''
    self.ZZ_sc1 = initialize_nans([1, len(self.angles), 360, 240])
    self.degr = numpy.arange(0,360, 1)
    self.r = numpy.arange(0, 240, 1)*1000  # convert from KM to M
    dx = (numpy.meshgrid(self.r, numpy.cos(numpy.deg2rad(self.degr)))[0] *
          numpy.meshgrid(self.r, numpy.cos(numpy.deg2rad(self.degr)))[1])
    self.lon = (self.LON_bilt + (180/numpy.pi)*(dx/self.r_earth)/numpy.cos(
      self.LAT_bilt*numpy.pi/180))
    dy = (numpy.meshgrid(self.r, numpy.sin(numpy.deg2rad(self.degr)))[0] *
          numpy.meshgrid(self.r, numpy.sin(numpy.deg2rad(self.degr)))[1])
    self.lat = self.LAT_bilt + (180/numpy.pi)*(dy/self.r_earth)
    meshgr = (numpy.meshgrid(numpy.ones(len(self.degr)),
                        numpy.sin(numpy.deg2rad(self.angles)), self.r))
    self.z = meshgr[0] * meshgr[1] * meshgr[2]
    try:
      h5file = h5py.File(self.filename,'r')
    except Exception:
      # TODO: add meaningfull exception handling
      import pdb; pdb.set_trace()
    for x in range(len(self.scans)):
      scan1 = h5file.get(self.scans[x])
      PV = numpy.array(scan1.get('scan_Z_data'))
      cal_sc1 = scan1.get('calibration')
      str = numpy.array_str(cal_sc1.attrs.get(
        'calibration_Z_formulas')).strip("['']").split('=')
      # evaluate expression in a "safe" manner
      Z_sc1 = asteval.Interpreter(symtable={"PV": PV}).eval(str[1])
      # set negative values equal to netcdf fill_value
      Z_sc1[Z_sc1<0] = -999
      for i in  range(0, len(self.r)):
        for j in self.degr:
          self.ZZ_sc1[0, x, j, i] = Z_sc1[j,i]

  def create_netcdf(self):
    '''
    Write netcdf output file
    '''
    ntime = 1; nangles = len(self.scans); ndegr = len(self.degr); ndist = len(self.r)
    # open output file
    ncfile = Dataset(self.outfile, 'w')
    # create dimensions
    ncfile.createDimension('time', ntime)
    ncfile.createDimension('angles', nangles)
    ncfile.createDimension('degrees', ndegr)
    ncfile.createDimension('distance', ndist)
    # create variables
    data = ncfile.createVariable('reflectivity','f4',
                                 ('time','angles','degrees','distance'),
                                 zlib=True, fill_value=-999)
    data1 = ncfile.createVariable('latitude', 'f4', ('degrees','distance'),
                                  zlib=True)
    data2 = ncfile.createVariable('longitude', 'f4', ('degrees','distance'),
                                  zlib=True)
    data3 = ncfile.createVariable('altitude', 'f4',
                                  ('angles','degrees','distance'),
                                  zlib=True)
    timevar = ncfile.createVariable('time', 'f4', ('time',),
                                    zlib=True)
    # time axis UTC
    dt = date2num(self.dt, calendar='gregorian',
                  units='minutes since 2010-01-01 00:00:00')
    # define attributes
    timevar.units = 'minutes since 2010-01-01 00:00:00'
    timevar.calendar = 'gregorian'
    timevar.standard_name = 'time'
    timevar.long_name = 'time in UTC'
    data3.units = 'meters'
    data.units = 'dBZ'
    data1.units = 'degrees_east'
    data1.standard_name = 'Longitude'
    data2.units = 'degrees_north'
    data2.standard_name = 'Latitude'
    # write data
    data[:] = self.ZZ_sc1
    data1[:] = self.lat
    data2[:] = self.lon
    data3[:] = self.z
    timevar[:] = dt
    # Add global attributes
    ncfile.description = 'KNMI radar data converted with fm128_radar_knmi'
    ncfile.history = 'Created ' + time.ctime(time.time())
    ncfile.radar_name = 'deBilt'
    ncfile.radar_longitude = '5.17834'
    ncfile.radar_latitude = '52.10168'
    ncfile.radar_height = '44'
    ncfile.radar_height_units = 'meters'
    # close output file
    ncfile.close()

if __name__=="__main__":
  dt = datetime.datetime(2014,1,1)
  convert_to_netcdf('RAD_NL60_VOL_NA_201401010000.h5', dt)
