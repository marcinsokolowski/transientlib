# transientlib
   
   Based on library for finding optical transients in "Pi of the Sky" data

# building :

  mkdir build
  cd build
  cmake ..
  make
  sudo make install

# Example usage after building 

  Get TwoLineElement file (for example from https://github.com/marcinsokolowski/Two_Line_Elements ), let's assume its name is 20230325.tle

  In order to generate a list of objects visible from the MRO at a particular UNIX TIME = 1663111525 execute the following command:

  sattest 1663111525 -all -qth=config/mro.qth -tle=20230325.tle -print_header -mwa

  in order to narrow down the search to objects within 5 degrees from a specific (RA,DEC) [deg] position, execute:

  sattest 1663111525 -all -qth=config/mro.qth -tle=20230325.tle -print_header -mwa  -ra=0.59 -dec=169.481659 -radius=5


  in order to see other options, execute :

  sattest -h
  
  

# Acknowledgements

   If you have been using this software or the related web interface for your research. Please, cite the paper:
   
    Advances in Astronomy, 2010, article id. 463496 (  https://ui.adsabs.harvard.edu/abs/2010AdAst2010E..54S/abstract )

    or 

    PhD thesis : Investigation of astrophysical phenomena in short time scales with "Pi of the Sky" apparatus
      ( PhD Thesis, 2008 : https://ui.adsabs.harvard.edu/abs/2008PhDT.......442S/abstract )
