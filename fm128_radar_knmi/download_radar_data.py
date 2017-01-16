'''
description:    Download KNMI radar hdf5 files (inside a tar archive) from ftp
license:        APACHE 2.0
author:         Ronald van Haren, NLeSC (r.vanharen@esciencecenter.nl)
'''

from ftplib import FTP
import os
import datetime


class download_radar_data:
  def __init__(self, dt, outputname='radar.tar'):
    '''
    Download radar data from KNMI ftp
    Each tar file contains one day of hdf5 files with a frequency of 15 minutes
    '''
    self.dt = dt  # datetime object
    self.outputname = outputname
    self.connect_to_ftp()
    self.change_to_download_directory()
    self.define_filename()
    self.open_output_file()
    self.download_file()
    self.close_output_file()

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
    self.ftp.retrbinary('RETR %s' % self.filename, self.file.write)

  def open_output_file(self):
    '''
    open output file for writing
    '''
    self.file = open(self.outputname, 'wb')

  def close_output_file(self):
    '''
    close the output file
    '''
    self.file.close()


if __name__=="__main__":
  download_radar_data(datetime.datetime(2014,1,1))
