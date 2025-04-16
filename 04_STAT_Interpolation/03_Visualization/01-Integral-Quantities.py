"""
Visualisation of the profiles 
@yuningw
"""
import matplotlib.pyplot as plt
import struct
import numpy as np
import pandas as pd
from   scipy import io as sio
from   scipy.interpolate import interp1d
from   scipy.integrate import quad
from  lib.plot import *  
from  lib.configs import *
import argparse
import copy
parser = argparse.ArgumentParser()
parser.add_argument('--x',default=0.75,type=float)
parser.add_argument('--s',default="SS",type=str)
args = parser.parse_args()
plt_setUp()

def calculate_relative(p,g):

  return (p - g )/g 

AOA = 5 
Rec = 200
fldr='../outputs_data/' 
sides = ['SS',"PS"]

#######
# del data['ppo']
#######


#######################################
# OVER SUCTION/PRESSURE SIDE 
#######################################
for caseName in data.keys():
    name = data[caseName]['fileName']
    fname= name_file(fldr,name,rec=Rec,naca=4412)
    d = sio.loadmat(fname)
    d = d[f'Re{Rec}k'][0][0]
    data[caseName]['data_SS'] = d['top'][0][0]['ref']
    data[caseName]['data_PS'] = d['bot'][0][0]['ref']
    data[caseName]['data_EX'] ={}
    print(f"[IO] DATA: {fname}")

#------------------------------------------
# Integral quantities : Cf 
#------------------------------------------
control_region_cfg = {
                    "xmin":0.40,
                    "xmax":0.41,
                    'color':cc.gray,
                    "alpha":0.5}

### ZOOM IN THE T.E For CF on S.S
var = 'Cf'
# fig,axs = plt.subplots(**single_fig_cfg)
fig,axss = plt.subplots(1,2,figsize=(15,8))
axs = axss[0]
# axins = inset_axes(axs, 
#                   width="100%", 
#                   height="100%",
#                   bbox_to_anchor=(  0.68,  0.57,   0.3,   0.4),
#                   bbox_transform=axs.transAxes,
#                           )
x_start = 0.1,
x_end  = 0.99
x_c_zoom = 0.35
x_c_end = 0.45
control_region_cfg2 = {
                    "xmin":x_c_zoom,
                    "xmax":0.86,
                    'color':cc.grays,
                    "alpha":0.5}

var_Name = var_name_dict[var]['name']
legend_list=[]
for case_name in data.keys():
  style_dict = copy.deepcopy(data[case_name]['style'])
  style_dict['marker'] = None
  style_dict['linestyle'] = "-"
  style_dict['lw'] = 3.0
  fig,axs = plot_integral_quantities(data[case_name]['data_SS'],fig,axs,
                                    x_start,x_end,
                                    var,var_Name,style_dict,interval=3)
  axs.set(**var_name_dict[var]['axs'])
  legend_list.append(data[case_name]['label'])

axs.yaxis.set_major_formatter(formatter2)
axs.grid(**grid_setup)
axs.axvspan(**control_region_cfg)
axs.axhline(0,**support_line1)

var_Name = var_name_dict[var]['name']
legend_list=[]
for case_name in data.keys():
  style_dict = copy.deepcopy(data[case_name]['style'])
  style_dict['marker'] = None
  style_dict['linestyle'] = "--"
  style_dict['lw'] = 3.0
  fig,axs = plot_integral_quantities(data[case_name]['data_PS'],fig,axs,
                                    x_start,x_end,
                                    var,var_Name,style_dict,
                                    interval=3,
                                    )
  axs.set(**var_name_dict[var]['axs'])
  legend_list.append(data[case_name]['label'])
axs.grid(**grid_setup)
# axs.axvspan(**control_region_cfg)
axs.yaxis.set_major_formatter(formatter2)


axs = axss[1]
var = 'R'
var_Name=r"$R \ [\%]$"

for case_name in data.keys():
  if 'reference' in case_name:
    ref_data = copy.deepcopy(data[case_name]['data_SS']['Cf'])

for case_name in data.keys():
  cre_data = copy.deepcopy(data[case_name]['data_SS']['Cf'])
  dff_data = calculate_relative(cre_data,ref_data)

  data[case_name]['data_EX']['R'] = dff_data
  data[case_name]['data_EX']['xa']  = data[case_name]['data_SS']['xa']


for case_name in data.keys():
  style_dict = copy.deepcopy(data[case_name]['style'])
  style_dict['marker'] = None
  style_dict['linestyle'] = "-"
  style_dict['lw'] = 3.0

  fig,axs = plot_integral_quantities(data[case_name]['data_EX'],fig,axs,
                                    x_start,x_end,
                                    var,var_Name,style_dict,interval=3)
  axs.set(**var_name_dict[var]['axs'])
  legend_list.append(data[case_name]['label'])

axs.yaxis.set_major_formatter(formatter2)
axs.grid(**grid_setup)
axs.axvspan(**control_region_cfg)
axs.axhline(0,**support_line1)

var_Name = var_name_dict[var]['name']
legend_list=[]
for case_name in data.keys():
  style_dict = copy.deepcopy(data[case_name]['style'])
  style_dict['marker'] = None
  style_dict['linestyle'] = "--"
  style_dict['lw'] = 3.0
  fig,axs = plot_integral_quantities(data[case_name]['data_EX'],fig,axs,
                                    x_start,x_end,
                                    var,var_Name,style_dict,
                                    interval=3,
                                    )
  axs.set(**var_name_dict[var]['axs'])
  legend_list.append(data[case_name]['label'])
axs.grid(**grid_setup)
fig.subplots_adjust(wspace=0.4)
# axs.axvspan(**control_region_cfg)
axs.yaxis.set_major_formatter(formatter2)





fig.savefig(f'Figs_{Rec}k/02-BL-DEVELOP/cf_cp_BothSides.jpg',
              **figkw
              )



#------------------------------------------
# Integral quantities : Others on both sides
#------------------------------------------
plt.rc("xtick",labelsize = 30)
plt.rc("ytick",labelsize = 30)
plt.rc("font",size = 35)



VarList =[
          'beta',
          'Reth',
          "Ret",
          "H12",
          ]

fig, axss = plt.subplots(**quadra_fig_22)
axss = axss.flatten()
AlphaList = [['(a)',"(b)","(c)","(d)"],["(e)","(f)","(g)","(h)"]]
for jl, side in enumerate(sides): 
  # fig,axss = plt.subplots(**quadra_fig_larger)
  for il, var in enumerate(VarList):
    axs=axss[il]
    # fig,axs = plt.subplots(**single_fig_cfg)
    x_c = 0.16
    var_Name = var_name_dict[var]
    legend_list=[]

    for case_name in data.keys():
      d = data[case_name][f'data_{side}']
      if 'ref'  in case_name:
        x_end = 0.86
      else:
        x_end = 0.92
      
      label = data[case_name]['label']
      style_dict = data[case_name]['style']
      style_dict['lw']  = 3.5
      style_dict['marker']  = None
      if side == 'SS':
        style_dict['linestyle']='-'
      elif side == 'PS':
        style_dict['linestyle']='--'
            
      fig,axs = plot_integral_quantities(d,
                                        fig,axs,x_c,x_end,
                                        var,var_Name,style_dict,
                                        interval=4
                                        )
      

      legend_list.append(data[case_name]['label'])
    axs.set(**var_name_dict[var]['axs'])
    axs.grid(**grid_setup)
    axs.axvspan(**control_region_cfg)
    axs.set_title(AlphaList[0][il],**title_setup)
    # axs.yaxis.set_major_formatter(formatter3)
fig.subplots_adjust(**{"hspace":0.4,"wspace":0.25})
fig.savefig(f'Figs_{Rec}k/02-BL-DEVELOP/BL_{il+1}VARS.jpg',
                **figkw
                )

