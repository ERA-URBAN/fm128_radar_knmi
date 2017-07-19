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
from geopy.distance import vincenty

class convert_to_ascii:
  def __init__(self, netcdf_file, outfile, pdry, interpolate=False):
    self.outputfile = outfile
    self.pdry = pdry
    self.check_file_exists(netcdf_file, 'r')
    self.read_netcdf(netcdf_file)
    if interpolate:
      self.interpolate(interpolate)
      self.set_mask_interpolate()
      self.write_ascii(single=True)
      #self.plot_refl()
    else:
      self.set_mask()
      self.write_ascii(single=False)

  def plot_refl(self):
    from mpl_toolkits.basemap import Basemap, cm
    import matplotlib.pyplot as plt
    # plot
    x1= numpy.min(self.longitude)
    x2 = numpy.max(self.longitude)
    y1 = numpy.min(self.latitude)
    y2 = numpy.max(self.latitude)

    for alt in range(0,len(self.altitude)):
      fig=plt.figure()
      map = Basemap(resolution='i',projection='merc',
                    llcrnrlat=y1,urcrnrlat=y2,llcrnrlon=x1,urcrnrlon=x2,
                    lat_ts=(x1+x2)/2)
      map.drawcoastlines(linewidth=0.25)
      map.drawcountries(linewidth=0.25)
      map.drawmeridians(numpy.arange(0,360,30))
      map.drawparallels(numpy.arange(-90,90,30))
      x,y = map(self.longitude, self.latitude)
      cs = map.pcolormesh(x,y, self.rf[alt,:],cmap=cm.s3pcpn_l, vmin=0, vmax=50)
      cbar = map.colorbar(cs, spacing='proportional',location='right',pad="5%")
      cbar.set_label('dBz')
      plt.savefig('rf_' + str(alt) + '.png', bbox_inches='tight')
      plt.close(fig)
      fig=plt.figure()
      map = Basemap(resolution='i',projection='merc',
                    llcrnrlat=y1,urcrnrlat=y2,llcrnrlon=x1,urcrnrlon=x2,
                    lat_ts=(x1+x2)/2)
      map.drawcoastlines(linewidth=0.25)
      map.drawcountries(linewidth=0.25)
      map.drawmeridians(numpy.arange(0,360,30))
      map.drawparallels(numpy.arange(-90,90,30))
      x,y = map(self.longitude, self.latitude)
      cs = map.pcolormesh(x,y, self.rv[alt,:],cmap=cm.s3pcpn_l, vmin=-20, vmax=20)
      cbar = map.colorbar(cs, spacing='proportional',location='right',pad="5%")
      cbar.set_label('m/s')
      plt.savefig('rv_' + str(alt) + '.png', bbox_inches='tight')
      plt.close(fig)
      fig=plt.figure()
      map = Basemap(resolution='i',projection='merc',
                    llcrnrlat=y1,urcrnrlat=y2,llcrnrlon=x1,urcrnrlon=x2,
                    lat_ts=(x1+x2)/2)
      map.drawcoastlines(linewidth=0.25)
      map.drawcountries(linewidth=0.25)
      map.drawmeridians(numpy.arange(0,360,30))
      map.drawparallels(numpy.arange(-90,90,30))
      x,y = map(self.longitude, self.latitude)
      cs = map.pcolormesh(x,y, self.altitude[alt,:],cmap=cm.s3pcpn_l, vmin=0, vmax=10000)
      cbar = map.colorbar(cs, spacing='proportional',location='right',pad="5%")
      cbar.set_label('m')
      plt.savefig('height_' + str(alt) + '.png', bbox_inches='tight')
      plt.close(fig)

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
      self.rf_err = 2.0 * numpy.ones(numpy.shape(self.rf))
    # radial velocity
    try:
      self.rv = ncfile.variables['radial_velocity'][0,:]
    except KeyError:
      self.rv = -88 * numpy.ones(numpy.shape(self.rf))
    try:
      self.rv_qc = ncfile.variables['radial_velocity_qc'][0,:]
    except KeyError:
      self.rv_qc = numpy.zeros(numpy.shape(self.rf))
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
    ncfile.close()

  def get_interpolate_grid(self, filename):
    '''
    Get grid to interpolate to from netcdf file containing XLAT/XLONG
    If XLAT/XLONG is (time, south_north_stag, west_east), use timestep 0
    '''
    ncfile = Dataset(filename, 'r')
    XLAT = ncfile.variables['XLAT'][0,:]
    if len(numpy.shape(XLAT))==1:
      XLAT = ncfile.variables['XLAT'][:]
    XLONG = ncfile.variables['XLONG'][0,:]
    if len(numpy.shape(XLONG))==1:
      XLONG = ncfile.variables['XLONG'][:]
    return XLAT, XLONG

  def interpolate(self, filename):
    '''
    interpolate to regular grid using nearest neighbor interpolation
    '''
    # get target grid from filename
    xlat, xlong = self.get_interpolate_grid(filename)
    # hardcode position of de Bilt for now
    LAT_bilt = 52.10168
    LON_bilt = 5.17834
    height_bilt = 44.
    # calculate horizontal distance for each xlat, xlong to de bilt
    # define source x,y,z
    x=self.longitude
    y=self.latitude
    z=self.altitude
    # source grid
    src=numpy.vstack((x.reshape(-1), y.reshape(-1), z.reshape(-1)))
    # source values, set last measurement in range to nan so we don't
    # extrapolate when using nearest neighbor
    self.rf[:,:,-1] = numpy.nan  # reflectivity
    vals_rf = self.rf.reshape(-1)
    self.rv[:,:,-1] = numpy.nan  # radial velocity
    vals_rv = self.rv.reshape(-1)
    self.rv_err[:,:,-1] = numpy.nan  # variance radial velocity
    vals_rv_err = self.rv_err.reshape(-1)
    # original shape
    orig_shape =  numpy.shape(self.rf)
    # define target coordinates
    xtrg = numpy.tile(xlong.reshape(-1), 13)
    ytrg = numpy.tile(xlat.reshape(-1), 13)
    # calculate horizontal distance for each gridpoint from base
    base_point = (LAT_bilt, LON_bilt)
    xlat_target = xlat.reshape(-1)
    xlon_target = xlong.reshape(-1)
    hor_dist = [vincenty(base_point, (xlat_target[x], xlon_target[x])).m for
                x in range(0, len(xlat_target))]
    ke = 1.3333333333  # adjustment factor to account for refractivity
    re = 6370040.0  # radius Earth
    angles = [0.4, 0.8, 1.1, 2.0, 3.0, 4.5, 6.0, 8.0, 10.0, 12.0,
              15.0, 20.0,  25.0 ]
    # calculate altitude for each beam angle at each grid point
    self.altitude = numpy.array([(numpy.sqrt((hor_dist/numpy.cos(numpy.deg2rad(theta)))**2
                                 +(ke*re)**2 + 2*(hor_dist/numpy.cos(
                                   numpy.deg2rad(theta)))*ke*re*
                                 numpy.sin(numpy.deg2rad(theta))
                                 ) - ke*re).reshape(numpy.shape(xlat)) for
                                 theta in angles]) + height_bilt
    ztrg = self.altitude[:,:].reshape(-1)
    trg = numpy.vstack((xtrg, ytrg, ztrg))
    # interpolate (using nearest neighbor)
    # scipy.interpolate.griddata is faster than wradlib.ipol
    self.rf = griddata(src.T, vals_rf, trg.T, method='nearest', rescale=True
                       ).reshape(numpy.shape(self.altitude))
    self.rv = griddata(src.T, vals_rv, trg.T, method='nearest', rescale=True
                       ).reshape(numpy.shape(self.altitude))
    self.rv_err = griddata(src.T, vals_rv_err, trg.T, method='nearest', rescale=True
                       ).reshape(numpy.shape(self.altitude))
    # use the single lon/lat value for the output we interpolated to
    self.longitude = xlong
    self.latitude = xlat
    # set the error to 2.0dBz
    self.rf_err = 2.0 * numpy.ones(numpy.shape(self.rf))
    self.rf_qc = numpy.zeros(numpy.shape(self.rf))
    self.rv_qc = numpy.zeros(numpy.shape(self.rf))

  def set_mask(self):
    '''
    set mask
    '''
    # create random mask with approx self.pdry% of data points
    # mask for altitudes above 10 km
    if self.pdry:
      if (self.pdry<=100) and (self.pdry>0):
        # correction factor on self.pdry to take into account dry points
        # max should not be larger than 100
        self.pdry = min(100,
                        len(self.rf.flatten())/float(
                          numpy.sum(self.rf<0)) * self.pdry)
        max_int = numpy.int(numpy.round((1/(self.pdry/100.))))
      else:
        max_int = 1
    else:
      max_int = 1
    # mask dry points
    mask = numpy.random.randint(0, max_int,
                                size=self.rf.shape).astype(numpy.bool)
    self.rf[self.rf<-30] = -30
    self.rv_err[self.rv<=-24] = -888888
    self.rv_qc[self.rv<=-24] = -88
    self.rv[self.rv<=-24] = -888888
    if not self.pdry==0:
      mask_dry = (mask) & (self.rf<0)
    else:
      mask_dry = (self.rf<0)
    # mask altitude
    mask_altitude = self.altitude > 6000
    # mask NaN
    mask_nan = numpy.isnan(self.rf)
    # combine masks and apply to reflectivity
    self.rf = numpy.ma.masked_where((mask_dry | mask_altitude | mask_nan),
                                    self.rf)
    self.rv = numpy.ma.masked_where((mask_dry | mask_altitude | mask_nan),
                                    self.rv)

  def set_mask_interpolate(self):
    '''
    set mask
    '''
    # create random mask with approx self.pdry% of data points
    # mask for altitudes above 10 km
    if self.pdry:
      if (self.pdry<=100) and (self.pdry>0):
        # correction factor on self.pdry to take into account dry points
        # max should not be larger than 100
        self.pdry = min(100,
                        len(self.rf[0,:].flatten())/float(
                          numpy.sum(self.rf[0,:]<0)) * self.pdry)
        max_int = numpy.int(numpy.round((1/(self.pdry/100.))))
      else:
        max_int = 1
    else:
      max_int =1
    # mask dry points
    mask_xy = numpy.random.randint(0, max_int,size=self.rf[0,:].shape
                                   ).astype(numpy.bool)
    mask = numpy.vstack([mask_xy[numpy.newaxis,:]] * numpy.shape(self.rf)[0])
    self.rf[self.rf<-30] = -30
    self.rv_err[self.rv<=-24] = -888888
    self.rv_qc[self.rv<=-24] = -88
    self.rv[self.rv<=-24] = -888888
    if not self.pdry==0:
      mask_dry = (mask) & (self.rf<=7)
    else:
      mask_dry = (self.rf<=7)
    # mask altitude
    mask_altitude = self.altitude > 6000
    # mask NaN
    mask_nan = numpy.isnan(self.rf)
    # combine masks and apply to reflectivity
    self.rf = numpy.ma.masked_where((mask_dry | mask_altitude | mask_nan),
                                    self.rf)
    self.rv = numpy.ma.masked_where((mask_dry | mask_altitude | mask_nan),
                                    self.rv)

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
