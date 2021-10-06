#!/bin/bash
########################################################################
####  script for parallization of the bHPE
########################################################################
# Loop over parallelisations
for parallel_count in {11..20}
do
    matlab -nodesktop -nosplash -r "Simulate_LKCestim3D_isotropic( 10 , 'bHPE', 'isotropic','num2str($parallel_count)' ); exit" &
done