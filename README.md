# HGCalValidator 

This package is devoted to the production and analysis of objects relevant to HGCAL. It is called 
HGCalValidator, since it monitors the majority of the objects we follow in the HGCalValidator 
module in CMSSW, however the architecture is different because here we want to produce trees. 

There are many things that are borrowed from the official unmaintained HGCAL tools, namely
[reco-prodtools](https://github.com/CMS-HGCAL/reco-prodtools), [reco-ntuples](https://github.com/CMS-HGCAL/reco-ntuples) and 
[ntuple-tools](https://github.com/CMS-HGCAL/ntuple-tools). These packages contain many things as a lot of users 
worked on them in the past. Here we try to limit ourselves to the absolutely necessary material. On the other hand, 
here we go far beyond as far as objects are concerned, reaching even associators. 

Due to the large amount of instructions to run this package, we decided to create the documentation in the 
wiki pages of this package ([link](https://github.com/apsallid/HGCalValidator/wiki)).

