'''
description:    Download KNMI radar hdf5 files (inside a tar archive) from ftp
license:        APACHE 2.0
author:         Ronald van Haren, NLeSC (r.vanharen@esciencecenter.nl)
'''

from netCDF4 import Dataset
from netCDF4 import num2date
import datetime
import numpy
from fm128_radar.write_fm128_radar import *
from scipy.interpolate import griddata

class convert_to_ascii:
  def __init__(self, netcdf_file, outfile, interpolate=False):
    self.outputfile = outfile
    self.check_file_exists(netcdf_file, 'r')
    self.read_netcdf(netcdf_file)
    if interpolate:
      self.interpolate()
      self.write_ascii(single=True)
    else:
      self.write_ascii(single=False)

  def check_file_exists(self, filepath, mode):
    ''' Check if a file exists and is accessible. '''
    try:
      f = open(filepath, mode)
      f.close()
    except IOError as e:
      raise IOError('File ' + filepath + ' is not accessible')

  def read_netcdf(self, netcdffile):
    '''
    Read the data from the netcdffile
    The netcdf file should contain at least reflectivity data
    '''
    ncfile = Dataset(netcdffile, 'r')
    # reflectivity
    try:
      self.rf = ncfile.variables['reflectivity'][0,:]
    except KeyError:
      raise KeyError('netCDF file ' + netcdffile +
                     ' does not contain any reflectivity data')
    try:
      self.rf_qc = ncfile.variables['reflectivity_qc'][0,:]
    except KeyError:
      self.rf_qc = numpy.zeros(numpy.shape(self.rf))
    try:
      self.rf_err = ncfile.variables['reflectivity_err'][0,:]
    except KeyError:
      self.rf_err = 0.1 * self.rf
      self.rf_err[self.rf_err<0] = 0
    # radial velocity
    try:
      self.rv = ncfile.variables['radial_velocity'][0,:]
    except KeyError:
      self.rv = -88 * numpy.ones(numpy.shape(self.rf))
    try:
      self.rv_qc = ncfile.variables['radial_velocity_qc'][0,:]
    except KeyError:
      self.rv_qc = -88 * numpy.ones(numpy.shape(self.rf))
    try:
      self.rv_err = ncfile.variables['radial_velocity_err'][0,:]
    except KeyError:
      self.rv_err = -888888 * numpy.ones(numpy.shape(self.rf))
    # lon/lat/altitude/time
    self.latitude = ncfile.variables['latitude'][:]
    self.longitude = ncfile.variables['longitude'][:]
    self.altitude = ncfile.variables['altitude'][:]
    time_in = ncfile.variables['time']
    # extract radar information from global attributes
    self.lon0 = float(ncfile.radar_longitude)
    self.lat0 = float(ncfile.radar_latitude)
    self.elv0 = float(ncfile.radar_height)
    self.radar_name = ncfile.radar_name
    # convert integer time to correct string
    self.time = num2date(time_in[0], calendar=time_in.calendar,
                         units=time_in.units)

  def interpolate(self):
    '''
    interpolate to regular grid using nearest neighbor interpolation
    '''
    # define source x,y,z
    x=self.longitude
    y=self.latitude
    z=self.altitude
    # source grid
    src=numpy.vstack((x.reshape(-1), y.reshape(-1), z.reshape(-1)))
    # source values
    vals = self.rf.reshape(-1)
    # original shape
    orig_shape =  numpy.shape(self.rf)
    # define target coordinates
    xtrg = numpy.tile(self.longitude[0,:].reshape(-1), 13)
    ytrg = numpy.tile(self.latitude[0,:].reshape(-1), 13)
    ztrg = self.altitude[:,:].reshape(-1)
    trg = numpy.vstack((xtrg, ytrg, ztrg))
    # interpolate
    self.rf = griddata(src.T, vals, trg.T, method='nearest', rescale=True
                       ).reshape(orig_shape)
    # use the single lon/lat value for the output we interpolated to
    self.longitude = self.longitude[0]
    self.latitude = self.latitude[0]
    # set the error to 10%
    self.rf_err = 0.1 * self.rf
    self.rf_err[self.rf_err<0] = 0

  def write_ascii(self, single=False):
    '''
    Write data to fm128_radar ascii file
    '''
    write_fm128_radar(self.radar_name, self.lat0, self.lon0, self.elv0,
                      self.time, self.latitude, self.longitude, self.altitude,
                      self.rf, self.rf_qc, self.rf_err,
                      self.rv, self.rv_qc, self.rv_err, outfile=self.outputfile,
                      single=single)

if __name__ == "__main__":
  read_netcdf_file()
