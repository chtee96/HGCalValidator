#!/bin/bash

#Close by photons - No second particle will be generated. 
# ----- GSD -----
if [ $1 == "GSD" ]
then
  echo \#CE_E_Front_120um
  echo \#-----
  echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 320.99 --zMax 321.01 --rMin 54.99 --rMax 55.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_E_Front_120um 
  python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 320.99 --zMax 321.01 --rMin 54.99 --rMax 55.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_E_Front_120um 

  echo \#CE_E_Front_200um
  echo \#-----
  echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 320.99 --zMax 321.01 --rMin 89.99 --rMax 90.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_E_Front_200um 
  python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 320.99 --zMax 321.01 --rMin 89.99 --rMax 90.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_E_Front_200um 

  #echo \#CE_E_Front_300um - Old position not run anymore
  #echo \#---------
  #echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 320.99 --zMax 321.01 --rMin 134.99 --rMax 135.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_E_Front_300um 
  #python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 320.99 --zMax 321.01 --rMin 134.99 --rMax 135.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_E_Front_300um 

  echo \#CE_E_Front_300um_Var1
  echo \#---------
  echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 320.99 --zMax 321.01 --rMin 124.99 --rMax 125.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_E_Front_300um_Var1 
  python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 320.99 --zMax 321.01 --rMin 124.99 --rMax 125.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_E_Front_300um_Var1 

  echo \#CE_H_Coarse_300um
  echo \#-------
  echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 431.99 --zMax 432.01 --rMin 79.99 --rMax 80.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_300um 
  python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 431.99 --zMax 432.01 --rMin 79.99 --rMax 80.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_300um 

  echo \#CE_H_Coarse_Scint
  echo \#-------
  echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 431.99 --zMax 432.01 --rMin 179.99 --rMax 180.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint 
  python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 431.99 --zMax 432.01 --rMin 179.99 --rMax 180.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint 

  echo \#CE_H_Fine_120um
  echo \#-------
  echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 362.519 --zMax 362.521 --rMin 49.99 --rMax 50.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Fine_120um 
  python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 362.519 --zMax 362.521 --rMin 49.99 --rMax 50.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Fine_120um

  echo \#CE_H_Fine_200um
  echo \#--------
  echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 362.519 --zMax 362.521 --rMin 89.99 --rMax 90.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Fine_200um 
  python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 362.519 --zMax 362.521 --rMin 89.99 --rMax 90.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Fine_200um

  #echo \#CE_H_Fine_300um - Old position not run anymore
  #echo \#--------
  #echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 362.519 --zMax 362.521 --rMin 134.99 --rMax 135.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Fine_300um 
  #python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 362.519 --zMax 362.521 --rMin 134.99 --rMax 135.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Fine_300um

  echo \#CE_H_Fine_300um_Var1
  echo \#--------
  echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 362.519 --zMax 362.521 --rMin 124.99 --rMax 125.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Fine_300um_Var1 
  python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 362.519 --zMax 362.521 --rMin 124.99 --rMax 125.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Fine_300um_Var1 

  echo \#CE_H_Fine_Scint
  echo \#-------
  echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 408.99 --zMax 409.01 --rMin 179.99 --rMax 180.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Fine_Scint 
  python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 408.99 --zMax 409.01 --rMin 179.99 --rMax 180.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Fine_Scint

  # Some variations of the positions not run anymore
  # Maybe useful for the future though
  #echo \#CE_H_Fine_Scint_Var1
  #echo \#-------
  #echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 405.99 --zMax 406.01 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Fine_Scint_Var1
  #python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 405.99 --zMax 406.01 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Fine_Scint_Var1

  #echo \#CE_H_Fine_Scint_Var2
  #echo \#-------
  #echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 408.99 --zMax 409.01 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Fine_Scint_Var2
  #python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 408.99 --zMax 409.01 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Fine_Scint_Var2

  #echo \#CE_H_Coarse_Scint_Var1
  #echo \#-------
  #echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 426.99 --zMax 427.01 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint_Var1
  #python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 426.99 --zMax 427.01 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint_Var1

  #echo \#CE_H_Coarse_Scint_Var2
  #echo \#-------
  #echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 432.99 --zMax 433.01 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint_Var2
  #python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 432.99 --zMax 433.01 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint_Var2

  #echo \#CE_H_Coarse_Scint_4285
  #echo \#-------
  #echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 428.49 --zMax 428.51 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint_4285
  #python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 428.49 --zMax 428.51 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint_4285

  #echo \#CE_H_Coarse_Scint_4295
  #echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 429.49 --zMax 429.51 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint_4295
  #python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 429.49 --zMax 429.51 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint_4295

  #echo \#CE_H_Coarse_Scint_4305
  #echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 430.49 --zMax 430.51 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint_4305
  #python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 430.49 --zMax 430.51 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint_4305

  #echo \#CE_H_Coarse_Scint_4315
  #echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 431.49 --zMax 431.51 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint_4315
  #python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 431.49 --zMax 431.51 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint_4315

  #echo \#CE_H_Coarse_Scint_4325
  #echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 432.49 --zMax 432.51 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint_4325
  #python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 432.49 --zMax 432.51 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint_4325

  #echo \#CE_H_Coarse_Scint_4335
  #echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 433.49 --zMax 433.51 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint_4335
  #python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 433.49 --zMax 433.51 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint_4335

  #echo \#CE_H_Coarse_Scint_4345
  #echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 434.49 --zMax 434.51 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint_4345
  #python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 434.49 --zMax 434.51 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint_4345

  #echo \#CE_H_Coarse_Scint_4354
  #echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 435.39 --zMax 435.41 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint_4354
  #python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 435.39 --zMax 435.41 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint_4354

  #echo \#CE_H_Coarse_Scint_4364
  #echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 436.39 --zMax 436.41 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint_4364
  #python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --zMin 436.39 --zMax 436.41 --rMin 169.99 --rMax 170.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_H_Coarse_Scint_4364

fi

