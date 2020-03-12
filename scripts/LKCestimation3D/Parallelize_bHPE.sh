#!/bin/bash
########################################################################
####  script for parallization of the bHPE
########################################################################
# Loop over parallelisations
for parallel_count in {1..10}
do
    matlab -nodesktop -nosplash -r "Simulate_LKCestim3D_isotropic( 10 , methods, "isotropic",'num2str($parallel_count)' ); exit" &
done