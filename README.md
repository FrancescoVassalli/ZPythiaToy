# ZPythiaToy
A Toy Pythia Model for Z event multiplicity.
The histors folder is for making plots using the data generated from the generators folder. 
In the generators folder the main executable is the ZPythiaGenerator. This takes up to 5 arguments. 
  1: pT Hat Min
  2: Number of events 
  3: output filename - must end in .root
  4: a string arugment to control some processes of pythia -read the code
Note that the data may not be output in a way that files will link together on your system. 
  You may need to change some of the string litterals in the histors so that they find the output from the generators on your file system.
The ZStop generator is the same as ZPythiaGenerator but it prints the event anytime there is a track such that dphi<pi/2 and pt>40.
The minBias only does minbias pythia as defined by the CMS configuration:
  https://github.com/CmsHI/genproductions/blob/master/python/HI/minbias/Pythia8_MinBias_pp_5020GeV_cfi.py
  The arguments it takes are also different
 
