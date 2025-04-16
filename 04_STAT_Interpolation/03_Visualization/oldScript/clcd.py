"""
Load and plot the Cl Cd 
@yuningw
"""
import  os
import pandas as pd 
import  numpy       as np 
import  scipy.io    as sio 
import  matplotlib.pyplot as plt 
from    utils.plot  import colorplate as cc 
from    utils       import plt_rc_setup
import matplotlib.ticker as ticker
import argparse
figkw  ={'bbox_inches':"tight", "dpi":300}


df  =pd.read_csv('Cl_Cd.csv')

######################
# Fig 1 Bar for Cd
########################
fig, axs = plt.subplots(1,1,figsize = (6,6))
caselist = [c for c in df['CASE']]
caselist = [c for c in df['CASE']]
print(caselist)


bottom = np.zeros(shape=(len(caselist)))
b1 = axs.bar(caselist,df['Cdf'],bottom=bottom,color = cc.grays,label='Cf')
bottom += df['Cdf']
b2 = axs.bar(caselist,df['Cdp'],bottom=bottom,color = cc.gray,label='Cp')
axs.axhline(df["Cdf"][0],color = cc.red)
axs.axhline(df["Cd"][0],color = cc.red)
axs.set_ylabel(r'$C_d = C_f + C_p$',fontsize = 20)
axs.set_xlabel(r'CASE',fontsize = 20)
axs.legend()
fig.tight_layout()
fig.savefig('Figs_200k/cl_cd_bar.jpg',**figkw)



################
# Fig2 Scatter for Cl/Cd
################


colorList = [cc.grays, cc.green, cc.deepgreen, cc.lightblue,
            cc.deepblue, cc.pink, cc.yellow, cc.darkred, cc.red]
fig, axs = plt.subplots(1,1,figsize=(6,4))
for il, case in enumerate(caselist):
    if case != 'I':
        axs.plot(df['Cd'][il],df['Cl'][il],'o',markersize = 11.5,c=colorList[il])
    else:
        axs.plot(df['Cd'][il],df['Cl'][il],'*',markersize = 12.5,c=colorList[il])

    axs.text(df['Cd'][il]-1e-4,df['Cl'][il]-2e-3,case,fontsize = 12, ha='right',va='bottom',color=colorList[il])
axs.set_ylabel(r"$C_l$",fontsize = 20)
axs.set_xlabel(r"$C_d$",fontsize = 20)
# axs.set_yticks(np.linspace(0.75,0.95,7))
# axs.set_xticks(np.linspace(0.02,0.0235,3))
axs.grid(which='major')
fig.tight_layout()
fig.savefig('Figs_200k/cl_VS_cd.jpg',**figkw)

