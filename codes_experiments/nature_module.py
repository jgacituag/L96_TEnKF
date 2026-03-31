# -*- coding: utf-8 -*-
import numpy as np
from scipy import stats
import time
import os
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, 'fortran', 'model'))
sys.path.insert(0, os.path.join(ROOT, 'fortran', 'da'))

from model   import lorenzn as model          #Import the model (fortran routines)
from obsope  import common_obs as hoperator      #Import the observation operator (fortran routines)



def nature_run( conf ) :

   if 'RandomSeed' in conf.GeneralConf.keys() :
     print('Warning: Random seed is being fixed')
     np.random.seed( conf.GeneralConf['RandomSeed'] )

   #=================================================================
   # LOAD CONFIGURATION : 
   #=================================================================

   GeneralConf=conf.GeneralConf
   ModelConf  =conf.ModelConf
   ObsConf    =conf.ObsConf
   NatureConf =conf.NatureConf

   #=================================================================
   # INITIALIZATION : 
   #=================================================================
   #Get the number of parameters
   NCoef=ModelConf['NCoef']
   #Get the number of ensembles
   NEns=NatureConf['NEns']
   #Get the size of the large-scale state
   Nx=ModelConf['nx']
   #Get the size of the small-scale state
   NxSS=ModelConf['nxss']

   #Initialize model configuration, parameters and state variables.
   if not ModelConf['EnableSRF']    :
     XSigma=0.0
     XPhi=1.0
   else                             :
     XSigma=ModelConf['XSigma']
     XPhi  =ModelConf['XPhi']

   if not ModelConf['EnablePRF']    :
     CSigma=np.zeros(NCoef)
     CPhi=1.0
   else                     :
     CSigma=ModelConf['CSigma']
     CPhi  =ModelConf['CPhi']


   if not ModelConf['FSpaceDependent'] :
     FSpaceAmplitude=np.zeros(NCoef) 
   else                   :
     FSpaceAmplitude=ModelConf['FSpaceAmplitude']

   FSpaceFreq=ModelConf['FSpaceFreq']

   CRF0=np.zeros([NEns,NCoef])
   RF0 =np.zeros([Nx,NEns])

   X0=np.zeros((Nx,NEns))
   XSS0=np.zeros((NxSS,NEns))
   C0=np.zeros((Nx,NEns,NCoef))

   #Generate a random initial conditions and initialize deterministic parameters
   for ie in range(0,NEns)  :

      X0[:,ie]=ModelConf['Coef'][0]/2 + np.random.normal( size=Nx )

      for ic in range(0,NCoef) :
         C0[:,ie,ic]=ModelConf['Coef'][ic] + FSpaceAmplitude[ic]*np.cos( FSpaceFreq[ic]*2*np.pi*np.arange(0,Nx)/Nx )

   #=================================================================
   # RUN SPIN UP : 
   #=================================================================

   #Do spinup
   print('Doing Spinup')
   start = time.time()

 
   nt=int( NatureConf['SPLength'] / ModelConf['dt'] )    #Number of time steps to be performed.
   ntout=int( 2 )                       #Output only the last time.

   #Runge Kuta 4 integration of the LorenzN equations
   [XSU , XSSSU , DFSU , RFSU , SSFSU , CRFSU , CSU ]=model.tinteg_rk4( nens=1 , nt=nt ,  ntout=ntout                              , 
                                   x0=X0     , xss0=XSS0   , rf0=RF0     , phi=XPhi     , sigma=XSigma                          , 
                                   c0=C0     , crf0=CRF0   , cphi=CPhi   , csigma=CSigma, param=ModelConf['TwoScaleParameters'] ,
                                   nx=Nx     , ncoef=NCoef , dt=ModelConf['dt']         , dtss=ModelConf['dtss'] )
   
   print('Spinup up took', time.time()-start, 'seconds.')

   #=================================================================
   # RUN NATURE : 
   #=================================================================

   #Run nature
   print('Doing Nature Run')
   start = time.time()

   X0=XSU[:,:,-1]                  #Start large scale variables from the last time of the spin up run.
   XSS0=XSSSU[:,:,-1]              #Start small scale variables from the last time of the sipin up run.  
   CRF0=CRFSU[:,:,-1]              #Spin up for the random forcing for the parameters   
   RF0=RFSU[:,:,-1]                #Spin up for the random forcing for the state variables


   nt=int( NatureConf['Length'] / ModelConf['dt'] )                     #Number of time steps to be performed.
   ntout=int( nt / ObsConf['Freq'] )  + 1                               #Output the state every ObsFreq time steps.


   #Runge Kuta 4 integration of the LorenzN equations.
   [XNature , XSSNature , DFNature , RFNature , SSFNature , CRFNature , CNature ]=model.tinteg_rk4( nens=1 , nt=nt ,  ntout=ntout                      ,
                                                       x0=X0     , xss0=XSS0   , rf0=RF0     , phi=XPhi     , sigma=XSigma                          ,
                                                       c0=C0     , crf0=CRF0   , cphi=CPhi   , csigma=CSigma, param=ModelConf['TwoScaleParameters'] ,
                                                       nx=Nx     , ncoef=NCoef , dt=ModelConf['dt'] , dtss=ModelConf['dtss'] )


   print('Nature run took', time.time()-start, 'seconds.')

   #=================================================================
   # GENERATE OBSERVATIONS : 
   #=================================================================

   #Apply h operator and transform from model space to observation spaece.
   print('Generating Observations')
   start = time.time()

   #Get the total number of observations
   NObs = hoperator.get_obs_number( ntype=ObsConf['NetworkType'] , nx=Nx , nt=ntout ,
                                      space_density=ObsConf['SpaceDensity']  ,
                                      time_density =ObsConf['TimeDensity'] )

   #Get the space time location of the observations based on the network type and 
   #observation density parameters.
   ObsLoc = hoperator.get_obs_location( ntype=ObsConf['NetworkType'] , nx=Nx , nt=ntout , no=NObs ,
                                      space_density=ObsConf['SpaceDensity']  , 
                                      time_density =ObsConf['TimeDensity'] )

   #Assume that all observations are of the same type.
   ObsType=np.ones( np.shape(ObsLoc)[0] )*ObsConf['Type']

   #Set the time coordinate corresponding to the model output.
   TLoc=np.arange(1,ntout+1)

   #Get the observed value (without observation error)
   [YObs , YObsMask]=hoperator.model_to_obs( nx=Nx   , no=NObs   , nt=ntout , nens=1    ,
                                      obsloc=ObsLoc , x=XNature , obstype=ObsType    ,
                                      obsval=np.zeros(NObs) , obserr=np.ones(NObs)   ,
                                      xloc=ModelConf['XLoc']    , tloc= TLoc         ,
                                      gross_check_factor = 1.0e9 , low_dbz_per_thresh = 1.1 )

   #Get the time reference in number of time step since nature run starts.
   ObsLoc[:,1]=( ObsLoc[:,1] - 1 )*ObsConf['Freq']

   #Add a Gaussian random noise to simulate observation errors
   ObsError=np.ones( np.shape(YObs) )*ObsConf['Error']
   ObsBias =np.ones( np.shape(YObs) )*ObsConf['Bias']


   YObs = hoperator.add_obs_error(no=NObs ,  nens=1  ,  obs=YObs  ,  obs_error=ObsError  ,
                               obs_bias=ObsBias , otype = ObsConf['Type'] ) 

   print('Observations took', time.time()-start, 'seconds.')

   #=================================================================
   # SAVE THE OUTPUT : 
   #=================================================================

   if NatureConf['RunSave']   :
    
      FNature = DFNature + RFNature + SSFNature  #Total Nature forcing.
    
      filename=GeneralConf['DataPath'] + '/' + GeneralConf['NatureFileName']
      print('Saving the output to ' +  filename  )
      start = time.time()
    
      if not os.path.exists( GeneralConf['DataPath'] + '/' )  : 
         os.makedirs(  GeneralConf['DataPath'] + '/'  )

      #Save Nature run output
      np.savez_compressed( filename ,   XNature=XNature       , FNature=FNature
                      ,   CNature=CNature       
                      ,   YObs=YObs , NObs=NObs , ObsLoc=ObsLoc   
                      ,   ObsType=ObsType       , ObsError=ObsError 
                      ,   ModelConf=ModelConf   , NatureConf=NatureConf 
                      ,   ObsConf=ObsConf       , GeneralConf=GeneralConf 
                      ,   XSSNature=XSSNature )
   
   
      #Print XNature and XSSNature as a CSV
      fileout=GeneralConf['DataPath'] + '/XNature.csv' 
      np.savetxt(fileout, np.transpose( np.squeeze( XNature ) ), fmt="%6.2f", delimiter=",")
      #np.squeeze(XNature).tofile(fileout,sep=',',format='%6.2f' , newline='\n' )

      #Print XNature and XSSNature as a CSV
      fileout=GeneralConf['DataPath'] + '/XSSNature.csv' 
      np.savetxt(fileout, np.transpose( np.squeeze( XSSNature ) ), fmt="%6.2f", delimiter=",")

      print('Saving took ', time.time()-start, 'seconds.')