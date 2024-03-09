# https://stackoverflow.com/questions/18721762/matplotlib-polar-plot-is-not-plotting-where-it-should
# import numpy as np
# import matplotlib.pyplot as plt

from __future__ import print_function

import gc # garbage collector
import matplotlib
# matplotlib.use('Agg')

from string import Template
import os
import sys
import errno
import matplotlib.pyplot as plt

t = Template("$EDA2TV_PATH/coinc")
new_import_path=t.substitute(os.environ)
sys.path.append( new_import_path )
print("Added import path %s" % new_import_path)
import eda2_aavs2_concidence 
import numpy as np
import re

# options :
from optparse import OptionParser,OptionGroup

def read_text_file( filename , verbose=0, ncols=2 ) :
   az_list = []
   elev_list = []

   if os.path.exists(filename) and os.stat(filename).st_size > 0 :
      file=open(filename,'r')
      data=file.readlines()
      for line in data :
         if line[0] != "#" :
#            words = line.split(' ')
            words = re.split( '\s+' , line )

            if verbose > 0 : 
               print("DEBUG : line = %s -> |%s|%s|" % (line,words[0+0],words[1+0]))

            if ncols >= 2 :
               az = float(words[0+0])
               elev = float(words[1+0])

            az_list.append(az)
            elev_list.append(elev)

      file.close()
   else :
      print("WARNING : empty or non-existing file %s" % (filename))

   print("READ %d values from file %s" % (len(az_list),filename))
   

   return (az_list,elev_list)


def read_list_file( listfile ) :
   file=open(listfile,'r')

   # reads the entire file into a list of strings variable data :
   data=file.readlines()
   files=[]

   for line in data : 
      words = line.split(' ')

      if line[0] == '#' :
         continue

      if line[0] != "#" :
        filename=words[0+0]

      files.append(filename.strip())
      
   return (files)   

def generate_plot(Az,El,options,imagefile="out.jpg"):
    outfile=options.outdir + "/" + imagefile
    
    if not os.path.exists(outfile) :
       fig = plt.figure( figsize=(10, 10) ,  dpi = options.dpi )
#    fig = plt.figure( figsize=(100, 30), dpi = options.dpi, tight_layout=True)
    # ax = fig.add_subplot(111, aspect='equal')
       ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
       ax.set_theta_zero_location('N')
       ax.set_theta_direction(-1)
       ax.plot(np.deg2rad(Az),90-El,'*')
       ax.set_yticks(range(0, 90+10, 30))                   # Define the yticks
       yLabel = ['90', '60','30','0']
       ax.set_yticklabels(yLabel)
       ax.set_xticks(np.arange(0, np.pi*2, np.pi/2))        # Define the xticks
       xLabel = ['N' , 'E','S','W']
       ax.set_xticklabels(xLabel)
       # plt.text(-0.9,0.5,imagefile, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes,fontsize=100)
       plt.title(imagefile)
    #    plt.show()
#    plt.savefig( options.outdir + "/" + imagefile, format = options.image_format, dpi = fig.dpi)
       plt.savefig( options.outdir + "/" + imagefile , format = options.image_format, dpi = options.dpi )
    else :
       print("File %s already exists -> skipped" % (outfile))

def mkdir_p(path):
   try:
      os.makedirs(path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST:
         pass

def parse_options(idx=0):
   usage="Usage: %prog [options]\n"
   usage+='\tPlotting candidates\n'
   parser = OptionParser(usage=usage,version=1.00)

   parser.add_option('--in_format','--input_format','--informat','--inputformat',dest="input_format",default="events", help="Input format [default %default]. events or sattest or azel (2 columns)")
   parser.add_option('--all',action="store_true",dest="all",default=False, help="Create all images (even without candidates) [default %s]")
   parser.add_option('--dpi',dest="dpi",default=100, help="Image DPI [default %default]",type="int") # 25
   parser.add_option('--out_format','--out_type','--ext','--out_image_type',dest="image_format",default="jpg", help="Output image type [default %default]")
   parser.add_option('--outdir','--out_dir','--output_dir','--dir',dest="outdir",default="images/",help="Output directory [default %default]")
   parser.add_option('--imagefile','--image_file','--out_file','--out_image',dest="image_file",default="satellites.png",help="Output image file name [default %default]")
   # parser.add_option('--list_high_freq',dest="list_high_freq",default="cand_list_aavs2",help="High frequency candidates list [default %default]")
   # parser.add_option('--list_low_freq',dest="list_low_freq",default="cand_list_eda2",help="Low frequency candidates list [default %default]")

   (options, args) = parser.parse_args(sys.argv[idx:])
   mkdir_p( options.outdir )

   return (options, args)


if __name__ == '__main__':
   (options, args) = parse_options()

   listfile="fits_list_I_diff"
#   filename="chan_204_20220914T004101_I_diff_cand.txt"
   if len(sys.argv) > 1:   
      listfile = sys.argv[1]



   print("######################################################")
   print("PARAMETERS :")  
   print("######################################################")
   print("List file    = %s" % (listfile))
   print("Input format = %s" % (options.input_format))
   print("Outdir       = %s" % (options.outdir))      
   print("######################################################")

   print("Read filenames from list file %s:" % (listfile))
   fits_file_list = read_list_file( listfile )
   for fitsfile in fits_file_list :
      gc.collect() # explicit invokation of garbage collector to clean-up memory (otherwise uses memory like crazy!)
                   # see : https://stackoverflow.com/questions/1316767/how-can-i-explicitly-free-memory-in-python
      candfile=fitsfile.replace(".fits","_cand.txt")
      candidates=[]
      imagefile=options.image_file
      if options.input_format == "sattest" :
         imagefile=fitsfile.replace(".txt","_cand.jpg")         
         print("%s -> %s" % (fitsfile,candfile))
         candidates=eda2_aavs2_concidence.read_sattest_file( candfile, 0.00, 0.00 )
      elif options.input_format == "azel" :
         (az_list,elev_list) = read_text_file( candfile )
      else : 
         imagefile=fitsfile.replace(".fits","_cand.jpg")
         print("%s -> %s" % (fitsfile,candfile))
         candidates=eda2_aavs2_concidence.read_candidate_file( candfile )
      
      Az=[]
      El=[]
      if (candidates is not None and len(candidates) > 0) or options.all or (len(az_list)>0 and len(elev_list)>0) :
         for cand in candidates :
            Az.append(cand.azim_deg)
            El.append(cand.elev_deg)
            
         for i in range(0,len(az_list)):
            Az.append(az_list[i])
            El.append(elev_list[i])
                        
         Az= np.array(Az)
         El= np.array(El)

         generate_plot(Az,El,options=options,imagefile=imagefile)   
         candidates=None         
         gc.collect()
      else :
         print("\tNo candidates in %s -> ignored" % (candfile))
#      sys.exit(0)
         
      
#   sys.exit(0)   

   # FITSNAME X Y FLUX[Jy] SNR ThreshInSigma RMS RA[deg] DEC[deg] AZIM[deg] ELEV[deg] UXTIME IMAGE_MIN IMAGE_MAX
   # chan_204_20220914T004101_I_diff.fits 28 117 117.47 74.97 5.00 1.57 171.061511 3.552535 66.890791 31.935454 1663116061.00 -12.5540 117.4658
   # Az=[90, 180, 270,]
   # El=[20, 30, 10]
#   Az=[66.890791]
#   El=[31.935454]
#   Az= np.array(Az)
#   El= np.array(El)
#   generate_plot(Az,El,options)

