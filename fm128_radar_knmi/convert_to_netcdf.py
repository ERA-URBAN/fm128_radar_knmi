'''
description:    Download KNMI radar hdf5 files (inside a tar archive) from ftp
license:        APACHE 2.0
author:         Ronald van Haren, NLeSC (r.vanharen@esciencecenter.nl)
                Natalie Theeuwes
'''
from wradlib import georef
import datetime
import numpy
import datetime
from netCDF4 import Dataset
from netCDF4 import date2num
import time
import asteval
import time
import h5py


def initialize_nans(shape):
  '''
  Initialize an array with nans
  '''
  nans = numpy.empty(shape)
  nans.fill(numpy.nan)
  return nans


class convert_to_netcdf:
  def __init__(self, filename, dt, pdry, outfile='radar.nc'):
    self.filename = filename
    self.outfile = outfile
    self.dt = dt
    self.pdry = pdry
    self.define_constants()
    self.calculate_reflectivity()
    self.create_netcdf()

  def define_constants(self):
    self.LAT_bilt = 52.10168
    self.LON_bilt = 5.17834
    self.height_bilt = 44.
    self.r_earth = 6378137.
    self.angles = [0.3, 0.4, 0.8, 1.1, 2.0, 3.0, 4.5, 6.0, 8.0, 10.0, 12.0,
                   15.0, 20.0,  25.0 ]
    self.scans = ['scan1', 'scan2', 'scan3', 'scan4', 'scan5', 'scan6',
                  'scan7', 'scan8', 'scan9', 'scan10', 'scan11', 'scan12',
                  'scan13', 'scan14']
    self.r_points = 240
    self.az_points = 360

  def calculate_reflectivity(self):
    '''
    Calculate reflectivity for each point
    '''
    self.ZZ_sc1 = initialize_nans([1, len(self.angles), self.az_points,
                                   self.r_points])
    self.degr = numpy.arange(0, self.az_points, 1)
    # distance for small angles is 1km per point
    # ['scan1', 'scan2', 'scan3', 'scan4', 'scan5']:
    r = numpy.arange(1, self.r_points+1, 1)*1000  # convert from KM to M
    lon_s, lat_s, z_s = self.calculate_lon_lat_z(self.angles[0:5], r)
    # distance for large angles is 0.5km per point (scan6 and higher)
    r = numpy.arange(1, self.r_points+1, 1)*500  # convert from KM to M
    lon_l, lat_l, z_l = self.calculate_lon_lat_z(self.angles[5:], r)
    # combine
    self.lon = numpy.vstack((lon_s, lon_l))
    self.lat = numpy.vstack((lat_s, lat_l))
    self.z = numpy.vstack((z_s, z_l))
    try:
      h5file = h5py.File(self.filename,'r')
    except Exception:
      # TODO: add meaningfull exception handling
      raise
    for x in range(len(self.scans)):
      scan1 = h5file.get(self.scans[x])
      PV = numpy.array(scan1.get('scan_Z_data'))
      cal_sc1 = scan1.get('calibration')
      str = numpy.array_str(cal_sc1.attrs.get(
        'calibration_Z_formulas')).strip("['']").split('=')
      # evaluate expression in a "safe" manner
      Z_sc1 = asteval.Interpreter(symtable={"PV": PV}).eval(str[1])
      # create random mask with approx self.pdry% of data points
      if self.pdry:
        if self.pdry<=100:
          max_int = numpy.int(numpy.round((1/(self.pdry/100.))))
        else:
          max_int = 1
        #mask = numpy.random.randint(0, max_int,
        #                            size=Z_sc1.shape).astype(numpy.bool)
        # Set only points in mask that are <0 to NaN (corresponds to dry cases)
        #Z_sc1[(mask) & (Z_sc1<0)] = -999
      for i in  range(0, len(r)):
        for j in range(0, len(self.degr)):
          self.ZZ_sc1[0, x, j, i] = Z_sc1[j,i]

  def calculate_lon_lat_z(self, angles, r):
    '''
    calculate lon,lat,alt
    '''
    meshgr = (numpy.meshgrid(self.degr,  angles, r))
    # calculate lon/lat/altitude using wradlib.georef
    # which takes into account refraction
    lla = georef.polar2lonlatalt_n(meshgr[2].reshape(-1),
                                   meshgr[0].reshape(-1), meshgr[1].reshape(-1),
                                   (self.LON_bilt, self.LAT_bilt,
                                    self.height_bilt))
    # return to correct shape
    lon = lla[0].reshape((len(angles), len(self.degr), len(r)))
    lat = lla[1].reshape((len(angles), len(self.degr), len(r)))
    z = lla[2].reshape((len(angles), len(self.degr), len(r)))
    return lon, lat, z

  def create_netcdf(self):
    '''
    Write netcdf output file
    '''
    ntime = 1; nangles = len(self.scans); ndegr = len(self.degr);
    ndist = self.r_points
    # open output file
    ncfile = Dataset(self.outfile, 'w')
    # create dimensions
    ncfile.createDimension('time', ntime)
    ncfile.createDimension('angle', nangles-1)
    ncfile.createDimension('degree', ndegr)
    ncfile.createDimension('distance', ndist)
    # create variables
    data = ncfile.createVariable('reflectivity','f4',
                                 ('time','angle','degree','distance'),
                                 zlib=True, fill_value=-999)
    data1 = ncfile.createVariable('latitude', 'f4', ('angle', 'degree','distance'),
                                  zlib=True)
    data2 = ncfile.createVariable('longitude', 'f4', ('angle', 'degree','distance'),
                                  zlib=True)
    data3 = ncfile.createVariable('altitude', 'f4',
                                  ('angle','degree','distance'), zlib=True)
    data4 = ncfile.createVariable('angle', 'f4', ('angle',), zlib=True)
    data5 = ncfile.createVariable('degree', 'f4', ('degree',), zlib=True)
    timevar = ncfile.createVariable('time', 'f4', ('time',), zlib=True)
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
    data1.units = 'degree_east'
    data1.standard_name = 'longitude'
    data2.units = 'degree_north'
    data2.standard_name = 'latitude'
    data4.unites = 'degree'
    data4.description = 'angle with respect to the horizontal'
    data5.unites = 'degree'
    data5.description = 'angle with respect to viewing direction'
    # write data
    data[:] = self.ZZ_sc1[0,1:,:]
    data1[:] = self.lat[1:,:]
    data2[:] = self.lon[1:,:]
    data3[:] = self.z[1:]
    data4[:] = self.angles[1:]
    data5[:] = self.degr
    timevar[:] = dt
    # Add global attributes
    ncfile.description = 'KNMI radar data converted with fm128_radar_knmi'
    ncfile.history = 'Created ' + time.ctime(time.time())
    ncfile.radar_name = 'deBilt'
    ncfile.radar_longitude = '5.17834'
    ncfile.radar_latitude = '52.10168'
    ncfile.radar_height = self.height_bilt
    ncfile.radar_height_units = 'meters'
    # close output file
    ncfile.close()

if __name__=="__main__":
  dt = datetime.datetime(2014,1,1)
  convert_to_netcdf('RAD_NL60_VOL_NA_201401010000.h5', dt)
