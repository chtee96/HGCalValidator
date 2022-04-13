#!/bin/bash

thedate="20220327"

#Close by photons - No second particle will be generated. 

#--------------------------------------------------------------------------
# GSD 
#--------------------------------------------------------------------------

if [ $1 == "GSD" ]
then
  echo \#CE_E_Front_120um
  echo \#-----
  echo python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 1000 --queue tomorrow --partID 22 --nPart 1 --zMin 320.99 --zMax 321.01 --rMin 54.99 --rMax 55.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_E_Front_120um 
  python production_withdeltas.py --datTier GSD --nevts 10000 --evtsperjob 1000 --queue tomorrow --partID 22 --nPart 1 --zMin 320.99 --zMax 321.01 --rMin 54.99 --rMax 55.01 --gunMode closeby --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --tag CE_E_Front_120um 

fi

#--------------------------------------------------------------------------
# RECO 
#--------------------------------------------------------------------------
if [ $1 == "RECO" ]
then

  samples=(CE_E_Front_120um)
  #samples=(CE_E_Front_120um CE_E_Front_200um CE_E_Front_300um_Var1 CE_H_Coarse_300um CE_H_Coarse_Scint CE_H_Fine_120um CE_H_Fine_200um CE_H_Fine_300um_Var1 CE_H_Fine_Scint)
  #older samples just for historical reason
  #samples = (CE_E_Front_300um CE_H_Fine_300um CE_H_Fine_Scint_Var1 CE_H_Fine_Scint_Var2 CE_H_Coarse_Scint_Var1 CE_H_Coarse_Scint_Var2 CE_H_Coarse_Scint_4285 CE_H_Coarse_Scint_4295 CE_H_Coarse_Scint_4305 CE_H_Coarse_Scint_4315 CE_H_Coarse_Scint_4325 CE_H_Coarse_Scint_4335 CE_H_Coarse_Scint_4345 CE_H_Coarse_Scint_4354 CE_H_Coarse_Scint_4364)
  for i in "${samples[@]}"
  do
      echo "----- $i ------"
      echo python production_withdeltas.py --datTier RECO --evtsperjob 1000 --queue tomorrow --partID 22 --nPart 1 --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --gunMode closeby --date ${thedate} --tag $i
      python production_withdeltas.py --datTier RECO --evtsperjob 1000 --queue tomorrow --partID 22 --nPart 1 --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --gunMode closeby --date ${thedate} --tag $i
  done
    
fi

#--------------------------------------------------------------------------
# NTUP 
#--------------------------------------------------------------------------
if [ $1 == "NTUP" ]
then

  samples=(CE_E_Front_120um CE_E_Front_200um CE_E_Front_300um_Var1 CE_H_Coarse_300um CE_H_Coarse_Scint CE_H_Fine_120um CE_H_Fine_200um CE_H_Fine_300um_Var1 CE_H_Fine_Scint)
  #older samples just for historical reason
  #samples = (CE_E_Front_300um CE_H_Fine_300um CE_H_Fine_Scint_Var1 CE_H_Fine_Scint_Var2 CE_H_Coarse_Scint_Var1 CE_H_Coarse_Scint_Var2 CE_H_Coarse_Scint_4285 CE_H_Coarse_Scint_4295 CE_H_Coarse_Scint_4305 CE_H_Coarse_Scint_4315 CE_H_Coarse_Scint_4325 CE_H_Coarse_Scint_4335 CE_H_Coarse_Scint_4345 CE_H_Coarse_Scint_4354 CE_H_Coarse_Scint_4364)
  for i in "${samples[@]}"
  do
      echo "----- $i ------"
      echo python production_withdeltas.py --datTier NTUP --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --gunMode closeby --date ${thedate} --tag $i
      python production_withdeltas.py --datTier NTUP --evtsperjob 100 --queue tomorrow --partID 22 --nPart 1 --eosArea /eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies --gunMode closeby --date ${thedate} --tag $i
  done

fi





