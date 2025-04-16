"""
A script to load matlab file and plot the Velocity profile 
"""
import  os
import  numpy       as np 
import  scipy.io    as sio 
import  matplotlib.pyplot as plt 
from    utils.plot  import colorplate as cc 
from    utils       import plt_rc_setup
import matplotlib.ticker as ticker
import argparse

parser  = argparse.ArgumentParser()
parser.add_argument("--x","-x",default=0.5,type=float,help="The location on x/c")
args    = parser.parse_args()

################################
# Parameter Setup  
################################

## Physical 

x_start  = 0.2457; x_end = 0.861
# x_start  = 0.3; x_end = 0.6

base_dir    = "01_data/"
case_list   = [  
                # "public",
                "reference",
                "oppo_uni050",
                "oppo_uni125",
                'oppo',
                "oppo_2sides",
                "oppo_125",
                "oppo_175",
                "oppo_300",
                "oppo_500",
                "oppo_650",
                "blw_ss_01",
                # "body",
                ]


style_list  =   [
                # {"lw":3,    "c":cc.grays,   "linestyle":"-."},
                {"lw":3,      "c":cc.black,   "linestyle":"-"},
                {"lw":1.5,    "c":cc.purple,    "linestyle":"-",'marker':"o",'markersize':3.5},
                {"lw":1.5,    "c":cc.orange,    "linestyle":"-",'marker':"s",'markersize':3.5},
                {"lw":1.5,    "c":cc.red,    "linestyle":"-",'marker':"o",'markersize':3.5},
                {"lw":1.5,    "c":cc.yellow,    "linestyle":"-",'marker':"o",'markersize':3.5},
                {"lw":1.5,    "c":cc.brown,    "linestyle":"-",'marker':"v",'markersize':3.5},
                {"lw":1.5,    "c":cc.deepgreen,    "linestyle":"-",'marker':"v",'markersize':3.5},
                {"lw":1.5,    "c":cc.deepblue,    "linestyle":"-",'marker':"v",'markersize':3.5},
                {"lw":1.5,    "c":cc.lightblue,    "linestyle":"-",'marker':"v",'markersize':3.5},
                {"lw":1.5,    "c":cc.green,    "linestyle":"-",'marker':"v",'markersize':3.5},
                {"lw":2.5,    "c":cc.blue,      "linestyle":"-"},
                # {"lw":2.5,    "c":cc.yellow,   "linestyle":"-"},
]


Labels     =   [
                # "Vinuesa et al. 2018",
                "Benchmark",
                "Oppo Ctrl: "+ r"$y^+ = 5.0$",
                "Oppo Ctrl: "+ r"$y^+ = 12.5$",
                "Oppo Ctrl: " + r"$y^+ = 15.0$",
                "Oppo Ctrl: "+ r"$y^+ = 15.0$" + " at "+ r"S.S & P.S",
                "Oppo Ctrl: "+ r"$y^+ = 12.5$" + " if "+ r"x/c > 0.4",
                "Oppo Ctrl: "+ r"$y^+ = 17.5$" + " if "+ r"x/c > 0.4",
                "Oppo Ctrl: "+ r"$y^+ = 30.0$" + " if "+ r"x/c > 0.4",
                "Oppo Ctrl: "+ r"$y^+ = 50.0$" + " if "+ r"x/c > 0.4",
                "Oppo Ctrl: "+ r"$y^+ = 65.0$" + " if "+ r"x/c > 0.4",
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


def read_Cf(fileName,ifglob=False):
    """
    Read the Cf vs xa from the results
    """
    data     = sio.loadmat(fileName)
    print(f"The {fileName} loaded, has keys:\n {data.keys()}")
    data     = data['Re200k'][0][0] 
    top      = data['top'][0][0]["ref"]
    xaa      = np.concatenate([ top[0,i]["xa"] for i in range(top.shape[1]) ])
    xaa      = xaa.flatten()
    xaa_ind  = np.where((xaa >=0.15))
    # xaa_ind  = np.where((xaa >=0.24) & (xaa<=0.5))
    xaa      = xaa[xaa_ind]
    
    # Load Global Cf for local 
    if ifglob:
        Cf       = np.array([top[0,i]['Cf'].squeeze() for i in range(top.shape[1])])
    else:
        Cf       = np.array([top[0,i]['Cfinf'].squeeze() for i in range(top.shape[1])])
    
    Cf       = Cf[xaa_ind]
    print(len(Cf))
    return xaa, Cf 

def read_Utau(fileName,ifglob=False):
    """
    Read the Cf vs xa from the results
    """
    data     = sio.loadmat(fileName)
    print(f"The {fileName} loaded, has keys:\n {data.keys()}")
    data     = data['Re200k'][0][0] 
    top      = data['top'][0][0]["ref"]
    xaa      = np.concatenate([ top[0,i]["xa"] for i in range(top.shape[1]) ])
    xaa      = xaa.flatten()
    xaa_ind  = np.where((xaa >=0.15))
    # xaa_ind  = np.where((xaa >=0.24) & (xaa<=0.5))
    xaa      = xaa[xaa_ind]
    
    # Load Global Cf for local 
    if ifglob:
        Cf       = np.array([top[0,i]['ut'].squeeze() for i in range(top.shape[1])])
    else:
        Cf       = np.array([top[0,i]['ut'].squeeze() for i in range(top.shape[1])])
    
    Cf       = Cf[xaa_ind]
    print(len(Cf))
    return xaa, Cf 

################################
# Visualisation
################################
fig, axs    = plt.subplots(1,1, figsize =(12,6))
for ind, case in enumerate(case_list):

    fontdict = style_list[ind]
    fileName = base_dir + case + ".mat"

    xaa, Cf = read_Cf(fileName)    
    axs.plot(xaa, Cf, **fontdict)

axs.grid()
y_labels = axs.get_yticks()

axs.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
axs.set_ylim(0.0,0.015)
# axs.set_xlim(x_start,x_end)
axs.vlines(x_start,ymin=0,ymax=0.015,linestyles='dotted',colors=cc.grays)
axs.vlines(x_end,ymin=0,ymax=0.015,linestyles='dotted',colors=cc.grays)
axs.set_xlabel(r"$x/c$",fontsize=25)
axs.set_ylabel(r"$C_f$",fontsize=25)
axs.legend(Labels,loc='upper right',ncol=1,prop={'size':10})

fig.savefig(f"Figs_200k/Cf_oppo_naca4412_{x_start:.1f}_{x_end:.1f}.jpg", **figkw)
fig.savefig(f"Figs_200k/Cf_oppo_naca4412_{x_start:.1f}_{x_end:.1f}.pdf", **figkw)

fig, axs    = plt.subplots(1,1, figsize =(12,6))

indices = [0,1,2,3,4,5,10]

labels_ = [Labels[l] for l in indices]

for ind in (indices):
    case = case_list[ind]
    fontdict = style_list[ind]
    fileName = base_dir + case + ".mat"

    xaa, Utau = read_Utau(fileName)    
    if case == 'reference': ut = Utau
    Cf = (1 - (Utau/ut)**2 )*100
    axs.plot(xaa, Cf, **fontdict)

axs.grid()
y_labels = axs.get_yticks()

axs.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
axs.set_ylim(-15.0,100)
# axs.set_xlim(x_start,x_end)
axs.axvline(x_start,ls='dotted',c=cc.grays)
axs.axvline(x_end,ls='dotted',c=cc.grays)
axs.set_xlabel(r"$x/c$",fontsize=25)
axs.set_ylabel(r"$R$",fontsize=25)
axs.legend(labels_,loc='upper left',ncol=1,prop={'size':10})

fig.savefig(f"Figs_200k/R_oppo_naca4412_{x_start:.1f}_{x_end:.1f}.jpg", **figkw)
fig.savefig(f"Figs_200k/R_oppo_naca4412_{x_start:.1f}_{x_end:.1f}.pdf", **figkw)