"""
A dictionary of all the cases here 
The order should be sorted by the control range and intensity 
The periodic and uniform cases are separated here. 
@yuningw
"""

"""
Rules: 
Color  === CONTROL INTENSITY / Freq
MARKER === CONTROL REGION 
LINE   === CONTROL METHOD
"""

from lib.plot import cc

lw1     = 3.0 
lw2     = 2.0 
mksize1 = 7.5 
ls0     = ':'
ls1     = '--'
ls2     = '-'
ls3     = ':'
mktype1 = 'o'
mktype2 = 'x'
mktype21 = 'X'
mktype3 = 's'
mktype4 = 'D'
mktype5 = 'v'
mktype6 = '^'

cs  = [cc.deepgreen,cc.blue,cc.purple,cc.yellow,cc.red,cc.pink]
ms = ['o','D','v','s','^']




##################### Config for Final data  ##############
data=     {
        'reference':{
                    'fileName':'benchmark',
                    'style':{
                            'lw':lw1,
                            'c':cc.black,
                            'linestyle':"-",
                            'marker':None,
                            'markersize':mksize1,
                            'fillstyle':'none',
                        },
                        'label':'Ref',
                        'config':{
                        'Re':200,
                        'region':[(0.0,0.0),(0,0)],
                        'intensity':[0.0,0.0],
                        "freq":[0.0,0.0],
                        'side':[False,False], # (S.S,P.S), 1==Yes, 0==No
                },
                        "d99":0.0477,
                        },
        
      }
