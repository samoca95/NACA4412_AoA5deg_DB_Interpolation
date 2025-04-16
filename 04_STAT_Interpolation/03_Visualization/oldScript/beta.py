"""
A script to load matlab file and plot the Calusure coefficient beta
"""
import  os
import  numpy       as np 
import  scipy.io    as sio 
import  matplotlib.pyplot as plt 
from    utils.plot  import colorplate as cc 
from    utils       import plt_rc_setup
from matplotlib.ticker import LogLocator, NullFormatter
import matplotlib.ticker as ticker

import argparse

parser  = argparse.ArgumentParser()
parser.add_argument("--x","-x",default=0.5,type=float,help="The location on x/c")
args    = parser.parse_args()
wing = '4412'
## Physical 

x_start  = 0.2457; x_end = 0.8619

base_dir    = "01_data/"

base_dir    = "01_data/"
case_list   = [  
                # "public",
                "reference",
                'oppo',
                "blw_ss_01",
                # "body",
                ]


style_list  =   [
                # {"lw":3,    "c":cc.grays,   "linestyle":"-."},
                {"lw":3,      "c":cc.black,   "linestyle":"-"},
                {"lw":2.5,    "c":cc.red,    "linestyle":"-"},
                {"lw":2.5,    "c":cc.blue,      "linestyle":"-"},
                # {"lw":2.5,    "c":cc.yellow,   "linestyle":"-"},
]


Labels     =   [
                # "Vinuesa et al. 2018",
                "Benchmark",
                "Opposition Control",
                "Suction-side: Blowing "    +    r"$0.1\%U_{\infty}$",
                #  "Body Force ",
                ]

# Setup for output 
figkw  ={'bbox_inches':"tight", "dpi":300}

################################
# DataReader 
################################
def read_Glob_Cf_VinuesaData(fileName):
    """
    Reading Cf from the public database  
    """
    d        = sio.loadmat(fileName)
    Cf     = np.array(d['top2n']["Cf"][0]).flatten()
    xa     = np.array(d['top2n']["xa"][0]).flatten()
    return xa, Cf 


def read_beta(fileName,ifglob=False):
    """
    Read the Cf vs xa from the results
    """
    data     = sio.loadmat(fileName)
    print(f"The {fileName} loaded, has keys:\n {data.keys()}")
    data     = data['Re200k'][0][0] 
    top      = data['top'][0][0]["ref"]
    xaa      = np.concatenate([ top[0,i]["xa"] for i in range(top.shape[1]) ])
    xaa      = xaa.flatten()
    xaa_ind  = np.where((xaa >0.2) & (xaa <=0.9))
    xaa      = xaa[xaa_ind]
    
    Cf       = np.array([top[0,i]['beta'].squeeze() for i in range(top.shape[1])])
    Cf       = Cf[xaa_ind]
    print(len(Cf))
    return xaa, Cf 


################################
# Visualisation
################################
fig, axs    = plt.subplots(1,1, figsize =(6,4))
for ind, case in enumerate(case_list):

    fontdict = style_list[ind]
    fileName = base_dir + case + ".mat"
    
    xaa, betal = read_beta(fileName)    

    axs.semilogy(xaa, betal, **fontdict)

axs.grid()
axs.vlines(x_start,ymin=0,ymax=100,linestyles='dotted',colors=cc.grays)
axs.vlines(x_end,ymin=0,ymax=100,linestyles='dotted',colors=cc.grays)
axs.set_xlabel(r"$x/c$",fontsize=25)
axs.set_ylabel(r"$\beta$",fontsize=25)
axs.legend(Labels,loc='upper left')
# axs.set_title(r"Clauser pressure-gradient parameter")

fig.savefig(f"Figs_200k/Beta_naca4412.jpg", **figkw)
fig.savefig(f"Figs_200k/Beta_naca4412.pdf", **figkw)