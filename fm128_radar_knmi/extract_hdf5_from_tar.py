'''
description:    Download KNMI radar hdf5 files (inside a tar archive) from ftp
license:        APACHE 2.0
author:         Ronald van Haren, NLeSC (r.vanharen@esciencecenter.nl)
'''

#from fm128_radar_knmi import download_radar_data
import datetime
import tarfile


class extract_hdf5_from_tar:
  def __init__(self, tar_archive, dt):
    '''
    Extract hd5 file from tar archive
    '''
    self.tar_archive = tar_archive
    self.dt = dt
    self.check_tarfile()
    self.define_filename()
    self.get_single_file()

  def check_tarfile(self):
    '''
    Check if the tar_archive exists and is a valid tar file
    '''
    try:
      if not tarfile.is_tarfile(self.tar_archive):
        # file exists but is not a valid tar file
        raise tarfile.ReadError(self.tar_archive + ' is not a tar archive')
    except FileNotFoundError:
        # file does not exist
        raise tarfile.ReadError(self.tar_archive + ' is not found')

  def define_filename(self):
    '''
    Define hdf5 filename to extract
    '''
    basename = 'RAD_NL60_VOL_NA_'
    ext = '.h5'  # file extension
    year = str(self.dt.year).zfill(4)
    month = str(self.dt.month).zfill(2)
    day = str(self.dt.day).zfill(2)
    hour = str(self.dt.hour).zfill(2)
    minute = str(self.dt.minute).zfill(2)
    self.filename = basename + year + month + day + hour + minute + ext

  def get_single_file(self):
    '''
    Extract a single file from a tar archive
    '''
    # open tar file
    tf = tarfile.open(self.tar_archive)
    tf_names = tf.getnames()  # list of files in tar archive
    # check if file is in the tar archive
    if self.filename in tf_names:
      tf.extract(self.filename)  # extract single file from tar archive
    tf.close()  # close tar file

if __name__=="__main__":
  dt = datetime.datetime(2014,1,1)
  extract_hdf5_from_tar('radar.tar', dt)
