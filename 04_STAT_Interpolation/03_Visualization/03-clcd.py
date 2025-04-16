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


AOA = 11 
Rec = 200
fldr='../outputs_data/' 
sides = ['SS',"PS"]

#######
# del data['ppo']
#######

############## Functions #####################
def rel_modif(g,p,positive=True):
  """
  Calculate the relative change of Cl,Cd
  g : reference value 
  p : controlled value
  """
  if positive:
    return (100*(p-g)/g) 
  else:
    return (100*(g-p)/g)


def get_Aero(d:dict):
    """Get Aerodynamic Characteristics"""
    d = d[0,0]
    char_list = ['Cl',"Cd","Cd_p","Cd_f","LD"]
    data={var:d[var][0][0] for var in char_list}
    return data
###############################################


###############
# Make a table for dtaa 
###############
df = {
    "Name":[],
    "Cl":[],
    "dCl":[],
    "Cd_p":[],
    "dCd_p":[],
    "Cd_f":[],
    "dCd_f":[],
    "Cd":[],
    "dCd":[],
    "LD":[],
    "dLD":[],
      }
char_list = ['Cl',"Cd","Cd_p","Cd_f","LD"]

#######################################
# OVER SUCTION/PRESSURE SIDE 
#######################################
for caseName in data.keys():
    name = data[caseName]['fileName']
    fname= name_file(fldr,name,rec=Rec)
    d = sio.loadmat(fname)
    d = d[f'Re{Rec}k'][0][0]
    data[caseName]['data'] = get_Aero(d['top'][0][0]['ref'])

    df['Name'].append(caseName)
    if 'ref' in caseName:
      ref_d = {char:data[caseName]['data'][char] for char in char_list }
    
    for char in char_list:
      df[char].append(data[caseName]['data'][char])
      df["d"+char].append(rel_modif(
                                    g=ref_d[char],
                                    p=data[caseName]['data'][char]))
pd.DataFrame(df).to_csv('CLCd.csv',float_format="%.5f")
pd.DataFrame(df).to_latex('CLCd.tex',float_format="%.5f")





######################
# Fig 1 Bar for Cd
########################
fig, axss = plt.subplots(1,2,figsize = (14,6))
axs = axss[0]
caselist = [ data[k]['label']  for k in data.keys() ]
print(caselist)
bottom = np.zeros(shape=(len(caselist)))
b1 = axs.bar(caselist,df['Cd_f'],bottom=bottom,color = cc.grays,label=r'$C_{d,f}$')
bottom += df['Cd_f']
b2 = axs.bar(caselist,df['Cd_p'],bottom=bottom,color = cc.gray,label=r'$C_{d,p}$')
axs.axhline(df["Cd_f"][0],linestyle='-.',color = cc.red)
axs.axhline(df["Cd"][0]     ,linestyle='-.',color = cc.red)
axs.set_ylabel(r'$C_d = C_{d,f} + C_{d,p}$',fontsize = 20)
axs.set_xlabel(r'CASE',fontsize = 20)
axs.legend(bbox_to_anchor=(1.0,0.5,0.0,0.5))
# axs.legend(loc='upper left',ncol=2)
axs.set_title("(a)",**title_setup)


################
# Fig2 Scatter for Cl/Cd
################

legend_list = []
# fig, axs = plt.subplots(1,1,figsize=(6,6))
axs = axss[1]
for il, case in enumerate(data.keys()):
    
    style_dict = data[case]['style']
    labelName  = data[case]['label']
    
    if 'ref' in case:
      style_dict['marker']='o'
    
    axs.plot(df['Cd'][il],df['Cl'][il],
              linestyle='none',
              marker=style_dict['marker'],
              c=style_dict['c'],
              markersize=25)
    legend_list.append(Line2D([0],[0],
                            linestyle='none',
                            marker=style_dict['marker'],
                            c=style_dict['c'],
                            markersize=12,
                            label=labelName
                                  ))
# axs.xaxis.set_major_formatter(formatter2)
axs.set(**{
          'xlabel':r"$C_d$",
          # 'xlim':[0.051,0.056],
          "ylabel":r"$C_l$",
          'ylim':[0.88,0.92],
          })
axs.axhline(df["Cl"][0]     ,linestyle='-.',color = cc.grays)
axs.axvline(df["Cd"][0]     ,linestyle='-.',color = cc.grays)
axs.grid(**grid_setup)
axs.legend(
  handles = legend_list,
  loc ='upper left',
  ncol = len(legend_list)//2,
  # bbox_to_anchor=(1.,0.5,0.0,0.5),
  prop={'size':15}
  )
axs.set_title("(b)",**title_setup)
# axs.grid(which='major')
fig.tight_layout()
fig.savefig(f'Figs_{Rec}k/01-CTRL-EFFECT/cl_VS_Cd.jpg',**figkw)
