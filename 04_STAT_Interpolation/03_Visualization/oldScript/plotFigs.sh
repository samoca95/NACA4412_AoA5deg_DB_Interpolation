#!/bin/bash

bash cpData.sh

python beta.py

python cf.py 

python R.py 

python Reth_.py

python u_Vel.py --x 0.4
python u_Vel.py --x 0.75 

python v_Vel.py --x 0.4
python v_Vel.py --x 0.75 

python uu_Vel.py --x 0.4 
python uu_Vel.py --x 0.75 

python uv_Vel.py --x 0.4 
python uv_Vel.py --x 0.75 