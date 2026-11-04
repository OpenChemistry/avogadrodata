#!/bin/sh
# h2o GFN2-xTB single point at the fixed geometry (xtb 6.7.1); writes molden.input
xtb ../geometry.xyz --gfn 2 --molden > h2o.out
mv molden.input h2o.molden
