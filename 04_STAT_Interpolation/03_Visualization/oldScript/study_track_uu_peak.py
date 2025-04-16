"""
Check the location of the second peak appers in after x/c = 0.5 over suction side 
@yuningw
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
parser.add_argument("--x","-x",default=0.8,type=float,help="The location on x/c")
args    = parser.parse_args()

base_dir    = "01_data/"
case_list   = [  
                "public",
                ]


style_list  =   [
                {"lw":3,    "c":cc.grays,   "linestyle":"-."},
]


Labels     =   [
                "Vinuesa et al. 2018",
                ]

def read_uu_VinuesaData(fileName):
    """
    Reading Cf from the public database  
    """
    
    ############
    # Load Data
    ###########
    d        = sio.loadmat(fileName)
    # Velocity  fluctutaion 
    uu      = np.array([i for i in d['top2n']["uu"][0]]).squeeze()
    # X - Axis
    x     = np.array([i for i in d['top2n']["xa"][0]]).squeeze()
    # Wall-normal distribution 
    yn     = np.array([i for i in d['top2n']["yn"][0]]).squeeze()
    # Viscousity # Constant 
    nu     = np.array([i for i in d['top2n']["nu"][0]]).squeeze()
    # Friction velocity 
    ut     = np.array([i for i in d['top2n']["ut"][0]]).squeeze()
    # Print the shape of the array 
    print(x.shape, yn.shape, nu.shape ,uu.shape, nu.shape, ut.shape)
    
    #############
    # Inner-scaling 
    #############
    for il in range(yn.shape[0]):
        yn[il,:] = yn[il,:]*ut[il]/nu[il]
        uu[il,:] = uu[il,:]/ut[il]**2
    return x, yn, uu


fileName = base_dir + case_list[0] + ".mat"
xa, yp, uu_inn = read_uu_VinuesaData(fileName)

#########
## Visualize the inner-scaled uu at each x/c 
#########
# for il, x in enumerate(xa):
#     fig, axs = plt.subplots(1,1)
#     axs.semilogx(yp[il,:],uu_inn[il,:],"-",c=cc.blue,lw=2.5)
#     axs.set_title(f"x/c = {x:.3f}")
#     fig.savefig('uu_figs/Test_uu_' + '0'*(2-len(str(il)))+ str(il) + '.jpg') 
#     plt.clf()


##########
# Exploration:
# Find the second peak in outer regin of the uu after Maximum chamber  
##########

## Given: y+ >= 1e2, let's narrow the range and do another visualisation
y_maxuu = []
icount = 1
for il, x in enumerate(xa):
    # Only after maximum chamber 
    if x >=0.8: 
        fig, axs = plt.subplots(1,1)

        yp1 = yp[il,:]
        uu1 = uu_inn[il,:]

        isort = np.where((yp1>=30) & (yp1<=400))
        yp1 = yp1[isort]
        uu1 = uu1[isort]
        axs.plot(yp1,uu1,"-o",c=cc.red,lw=2.5)
        axs.set_title(f"x/c = {x:.3f}")
        fig.savefig('uu_figs/Sort_uu_' + '0'*(2-len(str(icount)))+ str(icount) + '.jpg') 
        plt.clf()
        plt.close(fig)
        icount +=1 

        idx = np.argmax(uu1)
        y_maxuu.append(yp1[idx])
    elif x>0.4 and x<=0.8:
        fig, axs = plt.subplots(1,1)
        yp1 = yp[il,:]
        uu1 = uu_inn[il,:]
        isort = np.where((yp1>=60) & (yp1<=200))
        yp1 = yp1[isort]
        uu1 = uu1[isort]
        axs.plot(yp1,uu1,"-o",c=cc.red,lw=2.5)
        axs.set_title(f"x/c = {x:.3f}")
        fig.savefig('uu_figs/Sort_uu_' + '0'*(2-len(str(icount)))+ str(icount) + '.jpg') 
        plt.clf()
        plt.close(fig)
        icount +=1 
        idx = np.argmax(uu1)
        y_maxuu.append(yp1[idx])
    else:
        idx = np.argmax(uu_inn[il,:])
        y_maxuu.append(yp[il,idx])
    


#########
## Visualize the inner-scaled uu at each x/c with notation of the location 
#########
for il, x in enumerate(xa):
    fig, axs = plt.subplots(1,1,figsize=(6,4))
    axs.semilogx(yp[il,:],uu_inn[il,:],"-",c=cc.blue,lw=2.5)
    
    uu_min, uu_max = 0.9 * uu_inn[il,:].min(), 1.1 * uu_inn[il,:].max()
    axs.axvline(y_maxuu[il],ymin=uu_min,ymax=uu_max,color=cc.cyan)
    axs.set_title(f"x/c={x:.3f}; " + r"$y^+$" + f"={y_maxuu[il]:.3f}")
    fig.savefig('uu_figs/Mark_uu_' + '0'*(2-len(str(il)))+ str(il) + '.jpg') 
    plt.clf()
    plt.close(fig)


#######
# Finaly visualize the max y+ as func of x/c 
fig, axs = plt.subplots(1,1,figsize = (8,3))
isort = np.where(xa>=0.4)
y_maxuu = np.array(y_maxuu)
axs.plot(xa[isort], y_maxuu[isort],"o-",lw = 2.5,c=cc.black)
axs.axhline(np.mean(y_maxuu[isort]),color=cc.grays,linestyle ='--')
axs.set_xlabel('x/c',fontdict={'size':20})
axs.set_ylabel(r"$y^+_{\rm max(\overline{u'u'}^+)}$",fontdict={'size':20})
axs.set_title("Mean value of "+r"$y^+_{\rm max(\overline{u'u'}^+)}$"+f" ={np.mean(y_maxuu[isort]):.1f}",fontdict={'size':20})
fig.savefig('Figs_200k/y+_maxuu.jpg',bbox_inches='tight',dpi=300)