"""
Visualize the Turbulent Budget of opposition control 
"""


import  os
import  numpy       as np 
import  scipy.io    as sio 
import  matplotlib.pyplot as plt 
from    utils.plot  import colorplate as cc 
from    utils       import plt_rc_setup
from matplotlib.ticker import LogLocator, NullFormatter

import argparse

parser  = argparse.ArgumentParser()
parser.add_argument("--x","-x",default=0.4,type=float,help="The location on x/c")
args    = parser.parse_args()
base_dir = '01_data/'

case_list   = [  
                "public",
                "reference",

                # 'oppo',
                # "body",
                ]


# Drawing the plot 

draw_dict_public = {
        "Pkp"   :{"lw":2.5,      "c":cc.blue,        "linestyle":"-.",    "label":"Production"},
        "Dkp"   :{"lw":2.5,      "c":cc.red,         "linestyle":"-.",    "label":"Dissapation"},
        "Tkp"   :{"lw":2.5,      "c":cc.deepgreen,   "linestyle":"-.",    "label":"Turb.Transport"},
        "VDkp"  :{"lw":2.5,      "c":cc.brown,       "linestyle":"-.",    "label":"Visc.Diffusion"},
        "VPGkp" :{"lw":2.5,      "c":cc.black,       "linestyle":"-.",    "label":"Vel-PG Tensor"},
        "Ckp"   :{"lw":2.5,      "c":cc.pink,        "linestyle":"-.",    "label":"Convection"},
            }

draw_dict_self = {
        "Pkp"   :{"lw":2.5,      "c":cc.blue,    "linestyle":"-",    "label":"Production"},
        "Dkp"   :{"lw":2.5,      "c":cc.red,   "linestyle":"-",    "label":"Dissapation"},
        "Tkp"   :{"lw":2.5,      "c":cc.deepgreen,   "linestyle":"-",    "label":"Turb.Transport"},
        "VDkp"  :{"lw":2.5,      "c":cc.brown,   "linestyle":"-",    "label":"Visc.Diffusion"},
        "VPGkp" :{"lw":2.5,      "c":cc.black,   "linestyle":"-",    "label":"Vel-PG Tensor"},
        "Ckp"   :{"lw":2.5,      "c":cc.pink,   "linestyle":"-",    "label":"Convection"},
            }

legends = ["Production", 'Dissapation', 'Turb.Transport', 'Visc.Diffusion', "Vel-PG Tensor"]

# Setup for output 
figkw  ={'bbox_inches':"tight", "dpi":300}



def read_Budget_Vinuesa(fileName):
    """
    Reading Cf from the public database  
    """
    d        = sio.loadmat(fileName)
    
    data_dict = {}
    
    data_dict["xa"]     = np.array([i for i in d['top2n']["xa"][0]]).flatten()
    
    data_dict["yn"]     = np.array([i for i in d['top2n']["yn"][0]])
    data_dict["nu"]     = np.array([i for i in d['top2n']["nu"][0]]).flatten()
    data_dict["ut"]     = np.array([i for i in d['top2n']["ut"][0]]).flatten()
    
    data_dict["Pk"]     = np.array([i for i in d['top2n']["Pk"][0]])
    data_dict["Dk"]     = np.array([i for i in d['top2n']["Dk"][0]])
    data_dict["Tk"]     = np.array([i for i in d['top2n']["Tk"][0]])
    data_dict["VDk"]    = np.array([i for i in d['top2n']["VDk"][0]])
    data_dict["VPGk"]   = np.array([i for i in d['top2n']["VPGk"][0]])
    
    ID     = np.where(np.round(data_dict["xa"],2) == args.x)[0]
    ID_ = 0
    print(ID[ID_], data_dict['xa'][ID[ID_]])
    data_dict["yn"]     = data_dict["yn"][ID[ID_],:,:].flatten()
    data_dict["nu"]     = data_dict["nu"][ID[ID_]]
    data_dict["ut"]     = data_dict["ut"][ID[ID_]]
    data_dict["yp"]     = data_dict["yn"]*data_dict["ut"]/data_dict["nu"]
    id_y                = np.where((data_dict["yp"]>1.0))
    data_dict['yp']      = data_dict["yp"][id_y]
    
    # Inner-Scaling 
    data_dict["Pkp"]     = (data_dict['Pk'][ID[ID_],:,:].flatten()*data_dict['nu']/data_dict['ut']**4/2)[id_y]
    data_dict["Dkp"]     = (data_dict['Dk'][ID[ID_],:,:].flatten()*data_dict['nu']/data_dict['ut']**4/2)[id_y]
    data_dict["Tkp"]     = (data_dict['Tk'][ID[ID_],:,:].flatten()*data_dict['nu']/data_dict['ut']**4/2)[id_y]
    data_dict["VDkp"]    = (data_dict['VDk'][ID[ID_],:,:].flatten()*data_dict['nu']/data_dict['ut']**4/2)[id_y]
    data_dict["VPGkp"]   = (data_dict['VPGk'][ID[ID_],:,:].flatten()*data_dict['nu']/data_dict['ut']**4/2)[id_y]
    
    # # Inner-Scaling 
    # data_dict["Pkp"]     = data_dict['Pk'][id_y]
    # data_dict["Dkp"]     = data_dict['Dk'][id_y]
    # data_dict["Tkp"]     = data_dict['Tk'][id_y]
    # data_dict["VDkp"]    = data_dict['VDk'][id_y]
    # data_dict["VPGkp"]   = data_dict['VPGk'][id_y]
    
    del data_dict['nu'], data_dict['ut'], data_dict['xa'],data_dict['yn']
    del data_dict['Pk'], data_dict['Dk'], data_dict['Tk'],data_dict['VDk'],data_dict['VPGk']
    return data_dict 


def readBudget(fileName):
    data     = sio.loadmat(fileName)
    print(f"The {fileName} loaded, has keys:\n {data.keys()}")
    data     = data['Re200k'][0][0] 
    top      = data['top'][0][0]["ref"]
    xaa      = np.concatenate([ top[0,i]["xa"] for i in range(top.shape[1]) ]).flatten()
    
    ID       = np.where(np.round(xaa,2) == args.x)
    
    print(ID)
    ID       = ID[0][0]
    print(f"XAA= {xaa[ID]}")
    top_data = top[0,ID] 
    yp = top_data['yp'].flatten()
    id_y    =  np.where((yp >= 1.0)) 
    
    data_dict = {}
    data_dict['yp']      = yp[id_y]
    data_dict["Pkp"]     = np.reshape(top_data['Pkp'],(-1,))[id_y]
    data_dict["Dkp"]     = np.reshape(top_data['Dkp'],(-1,))[id_y]
    data_dict["Tkp"]     = np.reshape(top_data['Tkp'],(-1,))[id_y]
    data_dict["VDkp"]    = np.reshape(top_data['VDkp'],(-1,))[id_y]
    data_dict["VPGkp"]   = np.reshape(top_data['VPGkp'],(-1,))[id_y]
    # data_dict["Ckp"]     = np.reshape(top_data['Ckp'],(-1,))[id_y]
    
    return data_dict

fig, axs = plt.subplots(1,1, figsize =(7,6))
for ind, case in enumerate(case_list):
    fileName = base_dir + case + ".mat"
    print(f"INFO: Reading: {fileName}")
    locmin = LogLocator(base=1.0,subs=np.arange(1,10)*0.1, numticks=100)
    axs.xaxis.set_minor_locator(locmin)
    axs.xaxis.set_minor_formatter(NullFormatter())
    
    
    if case == 'public':
        data_dict = read_Budget_Vinuesa(fileName)
        draw_dict = draw_dict_public
    else:
        data_dict = readBudget(fileName)
        draw_dict = draw_dict_self
    
    for key in data_dict:
        if key != 'yp':
            axs.semilogx(data_dict['yp'],data_dict[key],**draw_dict[key])
    
    axs.set_xscale("log")
    axs.set_xlabel(r"$y_{n}^{+}$",fontsize=25)
    axs.set_ylabel("TKE Budget")
    axs.grid(which="both",c=cc.gray)
    axs.set_title(f"x/c = {args.x:.2f}", fontsize = 25)
    axs.set_xlim(data_dict['yp'].min()*0.9,data_dict['yp'].max()*1.1)
fig.tight_layout()

from matplotlib.lines import Line2D

custom_lines = [Line2D([0], [0],linestyle='-.', lw=1.5,color=cc.black),
                Line2D([0], [0],linestyle="-" , lw=1.5,color=cc.black),
                ]

axs.legend(legends,loc = 'lower right')

# axs.legend(custom_lines, ['Vinuesa et al. 2018', "Self Study"],loc='upper right')
fig.savefig(f'Figs_200k/Budget_Comparison_{int(args.x*100)}.jpg',**figkw)
