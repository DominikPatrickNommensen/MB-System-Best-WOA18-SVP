#!/bin/bash


#create a raw datalist
mbm_makedatalist -S.all -Odatalist_raw.mb-1

#create a platform file for Kongsberg MBES
mbmakeplatform --swath=datalist_raw.mb-1 --verbose --output=Meteor_EM122.plf

#mbpreprocess the raw datalist and platform file
mbpreprocess --input=datalist_raw.mb-1 --verbose --platform-file=Meteor_EM122.plf

#create a mb59 datalist
mbm_makedatalist -S.mb59 -P-V

#create a datalist for future processed files
mbdatalist -Z

#extrac ray info
mbbackangle -A1 -A2 -Q -V

#set the Amplitude file and the side scan file in the parameterfiles of all items of the datalist
mbset -PAMPCORRFILE:datalist.mb-1_tot.aga

mbset -PSSCORRFILE:datalist.mb-1_tot.sga
