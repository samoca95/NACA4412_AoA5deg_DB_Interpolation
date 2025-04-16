#!/bin/bash


python 01_Plotcf.py
python 02_EnergyGain.py
python 03_mean_profile.py --x 0.4
python 03_mean_profile.py --x 0.75
python 04_Reynolds_Stress.py --x 0.4
python 04_Reynolds_Stress.py --x 0.75
python 05_uu_Vel.py --x 0.75
python 06_turbBudget.py --x 0.4
python 06_turbBudget.py --x 0.75