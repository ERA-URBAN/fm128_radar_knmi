'''
description:    Download KNMI radar hdf5 files (inside a tar archive) from ftp
license:        APACHE 2.0
author:         Ronald van Haren, NLeSC (r.vanharen@esciencecenter.nl)
'''

from ftplib import FTP
import os
import datetime
import tarfile

class download_radar_data:
  def __init__(self, dt, outputdir):
    '''
    Download radar data from KNMI ftp
    Each tar file contains one day of hdf5 files with a frequency of 15 minutes
      - dt: datetime object of day to download the tar archive
            for
      - outputdir: directory where outputfile should be saved
    '''
    self.dt = dt  # datetime object
    self.outputdir = outputdir
    self.connect_to_ftp()
    self.change_to_download_directory()
    self.define_filename()
    self.define_outputfile()
    if not self.check_file_exists():
      self.download_file()

  def connect_to_ftp(self):
    '''
    Connect to KNMI ftp server
    '''
    # connect to host, default port
    self.ftp = FTP('data.knmi.nl')
    self.ftp.login()  # user anonymous, passwd anonymous@    

  def change_to_download_directory(self):
    '''
    Change into the correct download directory on the ftp server
    '''
    self.year = str(self.dt.year).zfill(4)
    self.month = str(self.dt.month).zfill(2)
    self.day = str(self.dt.day).zfill(2)
    download_dir = os.path.join('download', 'radar_tar_volume_debilt', '1.0',
                                '0001', self.year, self.month, self.day)
    self.ftp.cwd(download_dir)

  def define_filename(self):
    '''
    Define the filename to download
    '''
    basename = "RAD60_OPER_O___TARVOL__L2__"
    ext = '.tar'
    day1 = datetime.datetime.strftime(self.dt, "%Y%m%d")
    day2 = datetime.datetime.strftime(self.dt +datetime.timedelta(days=1),
                                      "%Y%m%d")
    self.filename = (basename + day1 + 'T000000_' + day2 + 'T000000' + '_0001'
                     + ext)

  def download_file(self):
    '''
    Download the tar file from ftp server
    '''
    # Open output file for writing
    self.file = open(self.outputfile, 'wb')
    # retrieve file
    self.ftp.retrbinary('RETR %s' % self.filename, self.file.write)
    # close the output file
    self.file.close()

  def define_outputfile(self):
    '''
    Define name and location of the output file
    '''
    self.outputfile = os.path.join(self.outputdir, self.filename)

  def check_file_exists(self):
    '''
    Check if outputfile exists and is a tar file
    We don't want to redownload the tar file from ftp if the
    file is already there
    '''
    try:
      if not tarfile.is_tarfile(self.outputfile):
        # file exists but is not a valid tar file
        os.remove(self.outputfile, dir_fd=None)
        # remove output file
        return False
      else:
        return True
      except OSError as e:
          if e.errno == errno.ENOENT:
            return False
          else:
              raise


if __name__=="__main__":
  download_radar_data(datetime.datetime(2014,1,1))
