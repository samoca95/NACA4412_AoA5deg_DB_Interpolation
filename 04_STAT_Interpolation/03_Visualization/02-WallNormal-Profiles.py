"""
Visualisation of the turbulence statistics  
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

parser = argparse.ArgumentParser()
parser.add_argument('--x',default=0.75,type=float)
parser.add_argument('--s',default="SS",type=str)
args = parser.parse_args()
plt_setUp_Smaller()

sides = ['SS',"PS"]
AlphaList = [['(a)',"(b)","(c)","(d)"],["(e)","(f)","(g)","(h)"]]
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
    fname= name_file(fldr,name,rec=Rec)
    d = sio.loadmat(fname)
    d = d[f'Re{Rec}k'][0][0]
    data[caseName]['data_SS'] = d['top'][0][0]['ref']
    data[caseName]['data_PS'] = d['bot'][0][0]['ref']
    
    print(f"[IO] DATA: {fname}")


"""
Mean Velocity profiles at different streamwise location Inner and outer scale 
"""
def Visual_Mean_Vel():
    VarList =['U','V']
    AlphaList = [["(a)","(b)"],["(c)","(d)"],]
    scales=['inner','outer']
    for var in VarList:
        fig,axss = plt.subplots(**quadra_fig_22_large2)
        # axss = axss.flatten()
        for jl, scale in enumerate(scales):
            for il, side in enumerate(sides): 
                axs = axss[il,jl]
                x_c = args.x
                var_Name = var_name_dict[var+scale]
                legend_list=[]
                
                for case_name in data.keys():

                    fig,axs = plot_Vel(data[case_name][f'data_{side}'],
                                    fig,axs,x_c,var,var_Name,
                                    data[case_name]['style'],
                                    grid_setup,
                                    scale=scale,
                                    interval=10)
                    legend_list.append(data[case_name]['label'])
                
                axs.set(**var_name_dict[var+scale]['axs'])
                axs.set_title(AlphaList[il][jl] + " " + rf"$x/c={x_c}$"+", "+f"{side_text[side]}",**title_setup)
                axs.xaxis.set_minor_locator(locmin)
                axs.xaxis.set_major_locator(locmin)
                axs.xaxis.set_minor_formatter(NullFormatter())
                axs.yaxis.set_major_formatter(formatter2)
        fig.subplots_adjust(**{"hspace":0.4,"wspace":0.25})
        fig.savefig(f'Figs_{Rec}k/03-STATS/{var}_{int(x_c*100)}.jpg',
                        **figkw
                        )
        # plt.clf()
        # plt.close(fig)


"""
Reynolds Stresses profiles 
"""
def Visual_Reynolds_Stress():
    VarList =['uu',
            #   "vv",'ww','uv',
                ]
    scales=['inner','outer']
    AlphaList = [["(a)","(b)"],["(c)","(d)"],]
    side = 'SS'
    for var in VarList:
        fig,axss = plt.subplots(**double_fig_larger)
        for jl, scale in enumerate(scales):
            # for il, side in enumerate(sides): 
                axs = axss[jl]
                x_c = args.x
                var_Name = var_name_dict[var+scale]
                legend_list=[]
                
                for case_name in data.keys():

                    fig,axs = plot_ReynoldStress(data[case_name][f'data_{side}'],
                                    fig,axs,x_c,var,var_Name,
                                    data[case_name]['style'],
                                    grid_setup,
                                    scale=scale)
                    legend_list.append(data[case_name]['label'])
                
                axs.set(**var_name_dict[var+scale]['axs'])
                # axs.set_title(AlphaList[il][jl] + " " + rf"$x/c={x_c}$"+", "+f"{side_text[side]}",**title_setup)
                axs.set_title(AlphaList[0][jl] + " " + rf"$x/c={x_c}$"+", "+f"{side_text[side]}",**title_setup)
                axs.xaxis.set_minor_locator(locmin)
                axs.xaxis.set_major_locator(locmin)
                axs.xaxis.set_minor_formatter(NullFormatter())
                axs.yaxis.set_major_formatter(formatter2)
        fig.subplots_adjust(**{"hspace":0.4,"wspace":0.3})
        fig.savefig(f'Figs_{Rec}k/03-STATS/{var}_{int(x_c*100)}.jpg',
                        **figkw
                        )
        plt.clf()
        plt.close(fig)


def Visual_Reynolds_Stress_All():
    VarList =[
            'uu',
            "vv",'ww','uv']
    NameList =[
        r'$\overline{u_t^2}$',
        r'$\overline{v_n^2}$',r'$\overline{w^2}$',r'$\overline{u_tv_n}$',
            ]
    # VarList =["vv",'ww','uv']
    scales=['inner',
            'outer']
    AlphaList = [["(a)","(c)"],["(b)","(d)"],]
    colors = plt.get_cmap('plasma')
    colors = colors(np.linspace(0,1,1+len(VarList)))[:-1]
    
    legend_list=[]
    figs_cfg = {'nrows':1,'ncols':2,'figsize':(14,6)}
    fig,axss = plt.subplots(**figs_cfg)


    ## Inspection of axess 


    side = 'SS'
    for jl, scale in enumerate(scales):
        # for il, side in enumerate(sides): 
            axs = axss[jl]
            x_c = args.x
            
            if side =='PS':
                ## Inspection 
                axins = inset_axes(axs,
                            width="100%", 
                            height="100%",
                            bbox_to_anchor=(  0.70,  0.55,   0.3,   0.4),
                            bbox_transform=axs.transAxes,
                            )
                for case_name in data.keys():
                    stylelist = data[case_name]['style']
                    for kl, var in enumerate(VarList):
                        stylelist['c'] = colors[kl]

                        var_Name = var_name_dict[var+scale]
                        fig,axins = plot_ReynoldStress(data[case_name][f'data_{side}'],
                                        fig,axins,x_c,var,var_Name,
                                        stylelist,
                                        grid_setup,
                                        scale=scale,
                                        interval=10,)
                
                if scale=='inner':
                    axins.set(**{'xlim':[0,100],'xscale':'log',
                                # 'ylim':[-20,50]
                                })
                    # axins.set_xticks([10,50,100])
                elif scale=='outer':
                    axins.set(**{'xlim':[0,100],'xscale':'log',
                                #  'ylim':[-0.005,0.01]
                                })
                
                axins.xaxis.set_minor_locator(locmin)
                # axins.xaxis.set_major_locator(locmin)
                axins.xaxis.set_minor_formatter(NullFormatter())
                axins.yaxis.set_major_formatter(formatter2)
        

            for case_name in data.keys():
                stylelist = data[case_name]['style']
                for kl, var in enumerate(VarList):
                    stylelist['c'] = colors[kl]

                    var_Name = var_name_dict[var+scale]
                    fig,axs = plot_ReynoldStress(data[case_name][f'data_{side}'],
                                    fig,axs,x_c,var,var_Name,
                                    stylelist,
                                    grid_setup,
                                    scale=scale,
                                    interval=10,)
                    
            axs.set(**var_name_dict['uiuj'+scale]['axs'])
            axs.set_title(AlphaList[0][jl] + " " + rf"$x/c={x_c}$"+", "+f"{side_text[side]}",**title_setup)
            axs.xaxis.set_minor_locator(locmin)
            axs.xaxis.set_major_locator(locmin)
            axs.xaxis.set_minor_formatter(NullFormatter())
            axs.yaxis.set_major_formatter(formatter2)
    
    for kl,var in enumerate(VarList):
        legend_list.append( Rectangle(xy=(1, 0), width=3, height=3,
                            color=colors[kl],
                            label=NameList[kl]))
    for case_name in data.keys():
        st = data[case_name]['style']
        label = data[case_name]['label']

        legend_list.append(Line2D([0],[0],
                            color=cc.grays,
                            lw  = 1.5, 
                            ls  = st['linestyle'],
                            marker = st['marker'] if case_name !='Ref' else None,
                            markersize=6,
                            fillstyle=st['fillstyle'],

                            label=label,
                            )
                            )

    axss[-1].legend(handles=legend_list,ncol=1,
                        loc='upper left', 
                        bbox_to_anchor=(1, 0.9, 
                                        0.1,0.1), 
                        borderaxespad=0, 
                        frameon=False,)

    # axss[-1,-1].legend(handles=legend_list,ncol=1,
    #                     loc='upper left', 
    #                     bbox_to_anchor=(1, 0.9, 
    #                                     0.1,0.1), 
    #                     borderaxespad=0, 
    #                     frameon=False,)

    fig.subplots_adjust(**{"hspace":0.35,"wspace":0.25})

    # fig.savefig(f'Figs_{Rec}k/03-STATS/StressAll_{int(x_c*100)}.pdf',
    #                 **figkw
    #                 )
    fig.savefig(f'Figs_{Rec}k/03-STATS/StressAll_{int(x_c*100)}.jpg',
                    **figkw
                    )
    plt.clf()
    plt.close(fig)





if __name__ == "__main__":
    Visual_Mean_Vel()

    Visual_Reynolds_Stress()
    
    Visual_Reynolds_Stress_All()