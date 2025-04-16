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
parser.add_argument("--x","-x",default=0.4,type=float,help="The location on x/c")
args    = parser.parse_args()


# os.system('bash "cpData.sh"')

base_dir    = "01_data/"
base_dir    = "01_data/"
case_list   = [  
                # "public",
                "reference",
                'oppo',
                "blw_ss_01",
                "synjet_xc06-08_amp01U_freq27e-1",
                "synjet_xc08-10_amp01U_freq27e-1",
                ]


style_list  =   [
                # {"lw":3,    "c":cc.grays,   "linestyle":"-."},
                {"lw":3,      "c":cc.black,   "linestyle":"-"},
                {"lw":1.5,    "c":cc.red,    "linestyle":"-",},
                {"lw":1.5,    "c":cc.blue,      "linestyle":"-"},
                {"lw":1.5,    "c":cc.yellow,    "linestyle":"-",},
                {"lw":1.5,    "c":cc.deepgreen,    "linestyle":"-"},
                {"lw":1.5,    "c":cc.brown,    "linestyle":"-",'marker':"v",},
                {"lw":1.5,    "c":cc.deepblue,    "linestyle":"-",'marker':"v",'markersize':3.5},
                {"lw":1.5,    "c":cc.lightblue,    "linestyle":"-",'marker':"v",'markersize':3.5},
                {"lw":1.5,    "c":cc.green,    "linestyle":"-",'marker':"v",'markersize':3.5},
                # {"lw":2.5,    "c":cc.yellow,   "linestyle":"-"},
]


Labels     =   [
                # "Vinuesa et al. 2018",
                "Benchmark",
                "Oppo Ctrl: " + r"$y^+ = 15.0$",
                "Suction-side: Blowing "    +    r"$0.1\%U_{\infty}$",
                "Periodic Control: "+ r"$x/c = 0.6 \sim 0.8$" + " " + r"$A_c = 0.1\%U_{\infty}$",
                "Periodic Control: "+ r"$x/c = 0.8 \sim 1.0$" + " " + r"$A_c = 0.1\%U_{\infty}$",
                #  "Body Force ",
                ]

# Setup for output 
figkw  ={'bbox_inches':"tight", "dpi":300}
# figkw  ={"dpi":300}


################################
# DataReader 
################################
def read_U_VinuesaData(fileName):
    """
    Reading Cf from the public database  
    """
    d        = sio.loadmat(fileName)
    U      = np.array([i for i in d['top2n']["U"][0]])
    xa     = np.array([i for i in d['top2n']["xa"][0]]).flatten()
    yn     = np.array([i for i in d['top2n']["yn"][0]])
    nu     = np.array([i for i in d['top2n']["nu"][0]]).flatten()
    ut     = np.array([i for i in d['top2n']["ut"][0]]).flatten()

    print(U.shape,xa.shape,yn.shape,nu.shape,ut.shape)
    print(xa)
    ID     = np.where( np.round(xa,1) == args.x)[0]
    print(ID)
    print(xa[ID])
    yn     = yn[ID[0],:,:].flatten()
    U      = U[ID[0],:,:].flatten()
    nu     = nu[ID[0]]
    ut     = ut[ID[0]]
    yp     = yn*ut/nu
    U      = U/ut

    isort  = np.where((yp>=0.5)&(yp<=400))
    yp     = yp[isort]
    U     = U[isort]

    return yp, U 



def read_U_yp(fileName,ifglob=False):
    """
    Read the U vs xa from the results
    """
    data     = sio.loadmat(fileName)
    print(f"The {fileName} loaded, has keys:\n {data.keys()}")
    data     = data['Re200k'][0][0] 
    top      = data['top'][0][0]["ref"]
    xaa      = np.concatenate([ top[0,i]["xa"] for i in range(top.shape[1]) ])
    
    ID       = np.where( xaa >= args.x)
    ID       = ID[0][0]
    
    top_data = top[0,ID] 
    nu       = top_data['nu'][0]
    utz      = top_data['ut'][0]
    Ue       = top_data['Uinf'][0]
    yp      =  top_data['yp']
    
    id_y    =  np.where((yp >= 1)&(yp<=400))[0]
    id_y    = id_y[-1]
    
    datdict = {}
    datdict['ut'] = utz
    datdict['Uinf'] = Ue
    datdict['yp'] = yp[:id_y]
    datdict['U']  = top_data["U"][:id_y]
    
    datdict['V']  = top_data["V"][:id_y]
    

    datdict['Uinner'] = datdict["U"]/datdict['ut']
    datdict['Uouter'] = datdict["U"]/datdict['Uinf']
    
    datdict['Vinner'] = datdict["V"]/datdict['ut']
    datdict['Vouter'] = datdict["V"]/datdict['Uinf']
    
    return datdict


#### Visualisation
fig1, axs1      = plt.subplots(1,1, figsize =(6,6))
fig2, axs2      = plt.subplots(1,1, figsize =(6,6))
fig3, axs3      = plt.subplots(1,1, figsize =(6,6))
fig4, axs4      = plt.subplots(1,1, figsize =(6,6))
for ind, case in enumerate(case_list):
    fontdict = style_list[ind]
    fileName = base_dir + case + ".mat"


    locmin = LogLocator(base=1.0,subs=np.arange(1,10)*0.1, numticks=100)
    axs1.xaxis.set_minor_locator(locmin)
    axs1.xaxis.set_minor_formatter(NullFormatter())

    if case == "public":
        yp, U = read_U_VinuesaData(fileName)
        axs1.semilogx( yp,U , **fontdict )
    else:
        datdict = read_U_yp(fileName)
        axs1.semilogx( datdict["yp"],datdict["Uinner"] , **fontdict )
        axs2.semilogx( datdict["yp"],datdict["Uouter"] , **fontdict )
        axs3.semilogx( datdict["yp"],datdict["Vinner"] , **fontdict )
        axs4.semilogx( datdict["yp"],datdict["Vouter"] , **fontdict )
        
    
axs1.set_xscale("log")
axs1.set_xlabel(r"$y_{n}^{+}$",fontsize=25)
axs1.set_ylabel(r"$U_{t}^{+}$",fontsize=25)
axs1.grid(which="both",c=cc.gray)
axs1.set_title(f"x/c = {args.x:.2f}", fontsize = 25)

axs2.set_xscale("log")
axs2.set_xlabel(r"$y_{n}^{+}$",fontsize=25)
axs2.set_ylabel(r"$U_{t}/U_e$",fontsize=25)
axs2.grid(which="both",c=cc.gray)
axs2.set_title(f"x/c = {args.x:.2f}", fontsize = 25)

axs3.set_xscale("log")
axs3.set_xlabel(r"$y_{n}^{+}$",fontsize=25)
axs3.set_ylabel(r"$V_{n}^{+}$",fontsize=25)
axs3.grid(which="both",c=cc.gray)
axs3.set_title(f"x/c = {args.x:.2f}", fontsize = 25)

axs4.set_xscale("log")
axs4.set_xlabel(r"$y_{n}^{+}$",fontsize=25)
axs4.set_ylabel(r"$V_{n}/U_e$",fontsize=25)
axs4.grid(which="both",c=cc.gray)
axs4.set_title(f"x/c = {args.x:.2f}", fontsize = 25)


fig1.tight_layout()
fig2.tight_layout()
axs1.legend(Labels,loc='upper left',prop={'size':6})
axs2.legend(Labels,loc='upper left',prop={'size':6})
axs3.legend(Labels,loc='upper left',prop={'size':6})
axs4.legend(Labels,loc='upper left',prop={'size':6})

fig1.savefig(f"Figs_Jet/U_Vel_naca4412_{int(args.x*100)}_inner.pdf", 
            **figkw)
fig1.savefig(f"Figs_Jet/U_Vel_naca4412_{int(args.x*100)}_inner.jpg", 
            **figkw)
fig2.savefig(f"Figs_Jet/U_Vel_naca4412_{int(args.x*100)}_outer.pdf", 
            **figkw)
fig2.savefig(f"Figs_Jet/U_Vel_naca4412_{int(args.x*100)}_outer.jpg", 
            **figkw)

fig3.savefig(f"Figs_Jet/V_Vel_naca4412_{int(args.x*100)}_inner.pdf", 
            **figkw)
fig3.savefig(f"Figs_Jet/V_Vel_naca4412_{int(args.x*100)}_inner.jpg", 
            **figkw)
fig4.savefig(f"Figs_Jet/V_Vel_naca4412_{int(args.x*100)}_outer.pdf", 
            **figkw)
fig4.savefig(f"Figs_Jet/V_Vel_naca4412_{int(args.x*100)}_outer.jpg", 
            **figkw)
