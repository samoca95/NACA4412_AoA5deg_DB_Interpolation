"""
A script to load matlab file and plot the Velocity profile 
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
parser.add_argument("--x","-x",default=0.5,type=float,help="The location on x/c")
args    = parser.parse_args()


# os.system('bash "cpData.sh"')

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
def read_V_VinuesaData(fileName):
    d        = sio.loadmat(fileName)
    V      = np.array([i for i in d['top2n']["V"][0]])
    xa     = np.array([i for i in d['top2n']["xa"][0]]).flatten()
    yn     = np.array([i for i in d['top2n']["yn"][0]])
    nu     = np.array([i for i in d['top2n']["nu"][0]]).flatten()
    ut     = np.array([i for i in d['top2n']["ut"][0]]).flatten()

    print(V.shape,xa.shape,yn.shape,nu.shape,ut.shape)
    print(xa)
    ID     = np.where( np.round(xa,1) == args.x)[0]
    print(ID)
    print(xa[ID])
    yn     = yn[ID[3],:,:].flatten()
    V      = V[ID[3],:,:].flatten()
    nu     = nu[ID[3]]
    ut     = ut[ID[3]]
    yp     = yn*ut/nu
    V      = V/ut
    istort = np.where(yp<=400)
    yp     = yp[istort]
    V     = V[istort]
    return yp, V 


def read_V(fileName):
    data     = sio.loadmat(fileName)
    print(f"The {fileName} loaded, has keys:\n {data.keys()}")
    data     = data['Re200k'][0][0] 
    top      = data['top'][0][0]["ref"]
    xaa      = np.concatenate([ top[0,i]["xa"] for i in range(top.shape[1]) ])
    print(xaa.shape)

    ID       = np.where(np.round(xaa,2) == args.x)
    ID       = ID[0][0]
    print(f"XAA= {xaa[ID]}")

    top_data = top[0,ID] 
    print(f"Nu = {top_data['nu']}")
    nu       = top_data['nu'][0]
    utz      = top_data['ut'][0]
    print(f"The nu is = {nu}")

    V       =  top_data["V"]/utz
    print(len(V))
    yp      =  top_data['yp']
    id_y    =  np.where((yp >= 1)&(yp<=400))[0]
    id_y    = id_y[-1]

    yp,V = yp[:id_y],V[:id_y]
    return yp, V 

fig, axs = plt.subplots(1,1, figsize =(6,6))

for ind, case in enumerate(case_list):
    locmin = LogLocator(base=1.0,subs=np.arange(1,10)*0.1, numticks=100)    
    axs.xaxis.set_minor_locator(locmin)
    axs.xaxis.set_minor_formatter(NullFormatter())

    fontdict = style_list[ind]
    fileName = base_dir + case + ".mat"
    yp, V = read_V(fileName)
    axs.semilogx(yp, V, **fontdict)

axs.set_xscale("log")
axs.set_xlabel(r"$y_{n}^{+}$",fontsize=25)
axs.set_ylabel(r"$V_{n}^{+}$",fontsize=25)
axs.grid(which="both",c=cc.gray)
axs.set_title(f"x/c = {args.x:.2f}", fontsize = 25)
fig.tight_layout()

plt.legend(Labels,loc='upper left')
plt.savefig(f"Figs_200k/V_Vel_naca4412_{int(args.x*100)}.pdf", 
            **figkw)
plt.savefig(f"Figs_200k/V_Vel_naca4412_{int(args.x*100)}.jpg", 
            **figkw)