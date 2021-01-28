# STEP2_QUENCHING

In order to perform DM for CHLA you should 
1. STEP 1 : Estimate the DARK OFFSET of the CHLA 

- https://github.com/catsch/dark_offset_chla on one core 
  or 
- https://github.com/qjutard/dark_offset_chla on several cores

2. STEP 2 : Correct the profile from the DARK and estimate the quenching correction (Xing 2018, Terrats 2020, Sackmann 2008) 
You need : 

  => A working directory containing Argo B files and C files (Be careful the program will write a new CHLA_ADJUSTED in your B files !!!) 
  => WMO.DARK : a file containing the OFFSET of the CHLA (the CHLA will be corrected as CHLA - OFFSET ) of the WMO float of interest
  => coriolis_CHLA.list is the list of all the WMO floats you want to correct (for every entry, you need a WMO.DARK file)
  
  lance_STEP2 is my launching script :
      -It defines the pathway to the working directory
      -It reads the coriolis_CHLA.list to know the WMO 
      -It reads the working directory to build the list of the Bfiles : "liste_all_B" 
      
  MLD_calc.R is calculating the MLD
  IPAR_15_DEPTH.R is estimating the depth at which the PAR is equal to 15 (Limit under which there is no longer quenching) 
  RunningFilter.R is used to filter (mean, median) the data 
  Quenching_correction.R is the main program (you will need oce, stringR, ncdf4 library)
  
  The outputs are files ready for STEP3 : estimation of the SLOPE 
