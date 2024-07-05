#If you want to replicate results from Section 6.2 in the manuscript, follow these steps:

#First, open CPCM.txt file and run the file (it contains the CPCM function)
#Second, download bQCD-master.zip (https://github.com/tagas/bQCD) and loci-master.zip (from https://github.com/AlexImmer/loci) 
#Third, in order to run the codes coded in python, open a python modul, set directory to the extracted loci-master.zip and run the following chuck of code:

from causa.datasets import AN, LS, MNU, SIMG, ANs, CausalDataset, Tuebingen, SIM, LSs
from causa.heci import HECI
from causa.loci import loci

#Fourth, set directory to the extracted bQCD-master.zip file and run Data_generators.R and Baseline_methods_in_R.R 
#finally, open Simulations 2.txt and run it


#Note: this is not the most elegant code. The code in python for the methods is coded more elegantly in https://github.com/AlexImmer/loci and in R in https://github.com/tagas/bQCD
