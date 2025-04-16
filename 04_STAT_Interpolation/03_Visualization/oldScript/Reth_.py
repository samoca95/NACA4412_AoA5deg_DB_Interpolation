"""
A script to load matlab file and plot the Velocity profile 
"""
import  os
import  numpy       as np 
import  scipy.io    as sio 
from scipy.interpolate import interp1d
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

x_start  = 0.2457; x_end = 0.8619

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
                {"lw":3,      "c":cc.black,     "linestyle":"-",'marker':"o"},
                {"lw":2.5,    "c":cc.red,       "linestyle":"-",'marker':"o"},
                {"lw":2.5,    "c":cc.blue,      "linestyle":"-",'marker':"o"},
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


def read_ReTheta(fileName,if_Re = True):
    """
    Read the Cf vs xa from the results
    """
    data     = sio.loadmat(fileName)
    # print(f"The {fileName} loaded, has keys:\n {data.keys()}")
    data     = data['Re200k'][0][0] 
    top      = data['top'][0][0]["ref"]
    xaa      = np.concatenate([ top[0,i]["xa"] for i in range(top.shape[1]) ])
    xaa      = xaa.flatten()
    xaa_ind  = np.where((xaa >=0.15))
    xaa      = xaa[xaa_ind]
    
    if if_Re:
        Theta       = np.array([top[0,i]['Reth'].squeeze() for i in range(top.shape[1])])
    else: 
        Theta       = np.array([top[0,i]['theta'].squeeze() for i in range(top.shape[1])])

    Theta       = Theta[xaa_ind]
    # print(len(Theta))
    return xaa, Theta 



################################
# Visualisation
################################
Theta_min = 500 
Theta_max = 500
fig, axs  = plt.subplots(1,1, figsize =(8,4))


for ind, case in enumerate(case_list):

    fontdict = style_list[ind]
    fileName = base_dir + case + ".mat"

    xaa, Theta = read_ReTheta(fileName)    


    axs.plot(xaa, Theta, **fontdict)

    if Theta_max < Theta.max(): Theta_max = Theta.max()
    if Theta_min > Theta.min(): Theta_min = Theta.min()

axs.grid()
y_labels = axs.get_yticks()

axs.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
axs.set_ylim(Theta_min*0.9, Theta_max*1.1)
axs.vlines(x_start,ymin=Theta_min*0.9,ymax=Theta_max*1.1,linestyles='dotted',colors=cc.grays)
axs.vlines(x_end,ymin=Theta_min*0.9,ymax=Theta_max*1.1,linestyles='dotted',colors=cc.grays)
axs.set_xlabel(r"$x/c$",fontsize=25)
axs.set_ylabel(r"$Re_{\theta}$",fontsize=25)
axs.legend(Labels,loc='upper left',prop={'size':12})

fig.savefig(f"Figs_200k/Theta_X_naca4412.jpg", **figkw)
fig.savefig(f"Figs_200k/Theta_X_naca4412.pdf", **figkw)



# #### Figure 2 The displacement of momentum thickness 
# Theta_min = 0 
# Theta_max = 0
# fig, axs  = plt.subplots(1,1, figsize =(8,4))

# datadict = {}
# for ind, case in enumerate(case_list):

#     fontdict = style_list[ind]
#     fileName = base_dir + case + ".mat"

#     xaa, Theta = read_ReTheta(fileName,False)    

#     xaa_line = np.linspace(xaa.min(), xaa.max(), 70)
#     interp_Theta = interp1d(xaa,Theta)(xaa_line)

#     datadict[case] = interp_Theta
    
#     dTheta = np.zeros_like(xaa_line)
#     print(f"At case = {case}")
#     for i, theta0 in enumerate(datadict['reference']):
#         j         = np.argmin((theta0 - datadict[case][i:])**2)
#         dTheta[i] = (xaa_line[j] - xaa_line[i])
        
#         print(f"For X = {xaa_line[i]}, The minial distance = {dTheta[i]}, The nearest = {xaa_line[j]}")
    
#     # dTheta = datadict["reference"] - datadict[case]

#     axs.plot(xaa_line, dTheta, **fontdict)
#     if Theta_max < dTheta.max(): dTheta_max = Theta.max()
#     if Theta_min > dTheta.min(): dTheta_min = Theta.min()

# axs.grid()
# y_labels = axs.get_yticks()

# axs.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
# # axs.set_ylim(Theta_min*0.9, Theta_max*1.1)
# # axs.vlines(x_start,ymin=Theta_min*0.9,ymax=Theta_max*1.1,linestyles='dotted',colors=cc.grays)
# # axs.vlines(x_end,ymin=Theta_min*0.9,ymax=Theta_max*1.1,linestyles='dotted',colors=cc.grays)
# axs.set_xlabel(r"$x/c$",fontsize=25)
# axs.set_ylabel(r"$\Delta x$",fontsize=25)
# axs.legend(Labels,loc='upper left',prop={'size':12})

# fig.savefig(f"Figs_200k/dTheta_X_naca4412.jpg", **figkw)
# fig.savefig(f"Figs_200k/dTheta_X_naca4412.pdf", **figkw)
