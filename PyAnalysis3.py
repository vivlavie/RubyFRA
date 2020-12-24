#PyAnalysis3.py
#Analysis using only Dataframe

import numpy as np
import math
#Distance of radiation extent to 6.3 kW/m2
from scipy.interpolate import interp1d

#for Event analyses
#Process Ara
Area = 'ProcessArea'
# element_dump_filename = 'Bv06_dump_PyExdCrv6'
element_dump_filename = 'Bv08_c5_dump'


#Offloaidng ara analysis
# 'Flammable dispersion' is not relevant??? All vaporized??
# Area = 'ReelStation'
# cExlFilename='Bv06_O_c'
# element_dump_filename = 'Bv06_offloading_dump'

#Utlity analysis
# cExlFilename='Bv06_u_c'
# # element_dump_filename = 'Bv06_utility_dump_2' #_2 with the updated S01 part count
# element_dump_filename = 'Bv06_utility_dump' #_2 with the updated S01 part count

#Hull deck analysis
# Area = 'Hull'
# cExlFilename='Bv06_h_c'
# element_dump_filename = 'Bv06_hull_dump'

import dill
element_dump_filename_DimCube5 = element_dump_filename + '_DimCube5'
with open(element_dump_filename_DimCube5,'rb') as element_dump:
# with open(element_dump_filename,'rb') as element_dump:
    lEvent = dill.load(element_dump)

Modules = []
for e in lEvent:    
    if e.Module not in Modules:
        Modules.append(e.Module)


def jffit(m):
    jl_lowe = 2.8893*np.power(55.5*m,0.3728)
    if m>5:
        jf = -13.2+54.3*math.log10(m)
    elif m>0.1:
        jf= 3.736*m + 6.
    else:
        jf = 0.
    # print(m, jl_lowe,jf)
    return jl_lowe


def PoolFireDuration(e,pool_type):
    if (pool_type == 'Early') and (e.EarlyPoolFire != None):
        dd = e.EarlyPoolFire.Diameter
    elif (pool_type == 'Late') and (e.LatePoolFire != None):
        dd = e.LatePoolFire.Diameter
    else:
        print(e.Key)
        print("Not a pool fire event")

    #Ms = e.Discharge.Ms[0] - e.Discharge.Ms[-1]
    m = e.TVD[:,1]
    Ms = max(m) - min(m)
    PFD = Ms/(3.14*dd*dd/4 * 0.062)*3/2 #Pool Fire Duration
    return PFD


#effect distance as ratio to thet jet length as per 'kW/m2'
# w_rad = {5:1.75, 12.5:1.45, 37.5:1.2}
j_rad = [5, 12.5, 27.5]
r_dist = [1.75, 1.45, 1.2]

#Distance of radiation extent to 6.3 kW/m2
from scipy.interpolate import interp1d
f_r_dist = interp1d(j_rad,r_dist)
d_6_3 = float(f_r_dist(6.3))

#Pool fire radiation distance
pool_diameter = [0.5, 2.5, 5, 10, 25, 30, 50, 80,100] #the lowest pool diamter and its radiation distance is added by SHK refering to ther ation at 2.5m
distance_12  = [1.26, 6.3, 9.4, 11, 13.2, 15.6, 25.2, 40.2, 50]
distance_5 = [2.4, 12, 19, 28, 30, 32, 43, 64, 77]
#Distance from the center of pool i.e. radius
f12 = interp1d(pool_diameter,distance_12)
f5 = interp1d(pool_diameter,distance_5)




import pandas as pd

listLeak = []
#To validate the radiation distance model by CMPT: compare JetFire.D12 & D37 
#TVD = [t, mm, rr ]      
for e in lEvent:
    key,hole,weather=e.Key.split("\\")
    Tf1 = e.TVD[-1,0]
    Tf2 = e.Discharge.Duration
    if "-L" in key:
        Fluid = 'L'
    else:
        Fluid = 'G'

    listLeak.append([key, hole, weather,e.Module,e.Deck,Fluid,e.Frequency,Tf1,e.Discharge.LiquidMassFraction])
pdLeak = pd.DataFrame(listLeak,columns=['Key','Hole','Weather','Module','Deck','Fluid','Frequency','Duration','LiquidMassFraction'])
pdLeakGrouped = pdLeak.groupby(['Module','Deck','Fluid'])
ftuples = [('Frequency','sum'),('LiquidMassFraction','mean')]
pdLeakGrouped['Frequency','LiquidMassFraction'].agg(ftuples)

# pdLeakGrouped['Frequency'].sum().plot.bar()
# plt.show() #This gives a bar graph with (# modules) x (# Decks) x (# of Fluid types) 



import matplotlib.pyplot as plt

pdJFM = pdJF.groupby('Module')
# aveLength = pdJFM['LengthFWeighted'].sum()/pdJFM['Frequency'].sum() #manual way of weighted average
get_wavg = lambda g: np.average(g['Length'],weights=g['Frequency'])
pdJFM.apply(get_wavg).plot.bar()
plt.show()



#Module-wise frequency analysis
pdF = pd.DataFrame(Frequency,index=['SM','ME','ME','LM'])
pdF = pdF.transpose()
pdF.hist() #to see how much frequencies are distributed...
plt.show()
pdF = pd.DataFrame(Frequency)
pdF = pdF.append(Modules,ignore_index=True)
pdF = pdF.append(Deck,ignore_index=True)
pdF.index = ['SM','ME','ME','LM','Module','Deck']
pdF = pdF.transpose()
pdFM = pdF.groupby(Modules)

def top(df, n=5, column='SM'):
    return df.sort_values(by=column)[-n:]

pdF_noNaN = pdF.dropna()

pdF_Process = pdF[pdF['Deck'].str.startswith("A")]
pdF_Process = pdF_noNaN[(pdF_noNaN.Deck in ["A","B","C"]).any(1)]
pdF_Process = pdF_noNaN[(pdF_noNaN.Deck == "A").any(1)]
pdFM_noNaN.apply(top) #gives top 5 high frequency of 'SM' per module


#String analysis in Dataframe
pdF_L = pdF.index.str.contains("131-01")
pdF_L = pdF[pdF_L]
pdF_LM = pdF_L.groupby(Modules)
pdF_LM.SM.max()
pdF_LM.SM.sum()
pdF_LM.SM.describe()


DC = {}
for e in lEvent:            
    pv,hole,weather = e.Key.split("\\")
    safety = hole[3:]    
    hole = hole[:2]        
    if e.JetFire != None:
        rr = e.Discharge.ReleaseRate
        dd = e.Discharge.Duration
        DC[(e.Module,e.Deck,pv,hole,safety,weather)] = [e.Pressure, e.Temperature, e.Frequency, rr, e.Discharge.LiquidMassFraction,dd]        
#Analysis for Radiation distances from PHAST and CMPT formula
pDC = pd.DataFrame(DC)
pDC = pDC.T
pDC.columns=['Pressure','Temperature','Frequency','ReleaseRate','LiqFrac','Duration']
pDC.index.names = ['Module','Deck','IS','Hole','Safety','Weather']
#Averge pressure per module-deck
pDC.mean(level=['Module','Deck']).Pressure.plot.bar()
plt.show()
# c1 = pDC['Module'] == 'S03'
# c2 = pDC['Module'] == 'P03'

pDC['VapRelRate'] = pDC['ReleaseRate'] * (1-pDC['LiqFrac'])
c3 = pDC['VapRelRate'] > 40

pDC.loc[['S03','P03']][c3].sum()

#Frequency of havign flammable gas from S03 and P03 to GTG exhaust
pDC.loc[['S03','P03']][c3].sum()/6*0.5/2
# Frequency        0.000263

 pdDC.xs('S02',level=0)
 pdDC.loc[('S02','B',slice(None))]
 pdDC.xs('S02',level=0)


listDC = []
for e in lEvent:            
    pv,hole,weather = e.Key.split("\\")
    safety = hole[3:]    
    hole = hole[:2]        
    if e.JetFire != None:
        rr = e.Discharge.ReleaseRate
        dd = e.Discharge.Duration
        listDC.append([e.Module,e.Deck,pv,hole,safety,weather,e.Pressure, e.Temperature, e.Frequency, rr, e.Discharge.LiquidMassFraction,dd])
#Analysis for Radiation distances from PHAST and CMPT formula
pdDC = pd.DataFrame(listDC)
pdDC.columns=['Module','Deck','Key','Hole','Safety','Weather','Pressure','Temperature','Frequency','ReleaseRate','LiqFrac','Duration']


is_ME1 = pdDC['Hole'] == 'ME'

is_MS2 = pdJFc['Hole'] == 'SM'
is_ME2 = pdJFc['Hole'] == 'ME'
is_MA2 = pdJFc['Hole'] == 'MA'
is_LA2 = pdJFc['Hole'] == 'LA'
is_LM2 = pdJFc['Hole'] == 'LM'

pdJFc[is_ME2].groupby('Module').Length.mean()
is_SM2 = pdJFc['Hole'] == 'SM'
pdJFc[is_SM2].groupby('Module').Length.mean()

pdSM = pdJFc[is_SM2]
pdME = pdJFc[is_ME2]
pd.merge(pdSM,pdME)
pdSM = pdJFc[is_SM2].groupby('Module')
pdME = pdJFc[is_ME2].groupby('Module')
pdMA = pdJFc[is_MA2].groupby('Module')
pdLA = pdJFc[is_LA2].groupby('Module').Length
pdLM = pdJFc[is_LM2].groupby('Module').Length
# #JF method 1 using a dictionary
# JF = {}
JF2 = {}
for e in lEvent:            
    pv,hole,weather = e.Key.split("\\")
    safety = hole[3:]    
    hole = hole[:2]
    if e.JetFire != None:
        jfl = e.JetFire.Length
        sep = e.JetFire.SEP
        dur = e.Discharge.Duration
        # JF[e.Key] = [pv,hole,safety,weather,e.Module,e.Deck,e.JetFire.Frequency,jfl,e.JetFire.D37,1.2*jfl, e.JetFire.D12,1.45*jfl,e.JetFire.D04,1.75*jfl]        
        JF2[(e.Module,e.Deck,pv,hole,safety,weather)] = [e.JetFire.Frequency,jfl,sep,dur]        
#Analysis for Radiation distances from PHAST and CMPT formula
pJF = pd.DataFrame(JF2)
# pJF = pJF.replace('n/a',np.nan)
pJF = pJF.T
# pJF.columns=['IS','Hole','Safety','Weather','Module','Deck','JFF','JFL','D37p','D37c','D12p','D12c','D04p','D05c']
pJF.columns=['JFF','JFL','SEP','Duration']
pJF.index.names = ['Module','Deck','IS','Hole','Safety','Weather']
# #to check all are correctly read
# for i in pJFW.index:
#     if np.nan in i:
#         print(i)

plt.subplot(211)
pJF.SEP.max(level='Module').plot.bar()
plt.xlabel(None)
plt.xticks([])
plt.ylim([0,390])
# pJF.SEP.describe(level='Module').plot.bar()
plt.title('Surface Emissive Power [kW/m2]\nMaximum')
plt.subplot(212)
pJF.SEP.mean(level='Module').plot.bar()
plt.title('Average')
plt.ylim([0,390])
plt.xticks(rotation=45,ha='right')
plt.tight_layout()
plt.show()


pJF.groupby('Module').SEP.describe()


pJFW = pJF.unstack() #Jet Fire with top column for 'W'eather
# pJFW.loc['P05','A',:,'SM']
# pJFW.loc['P04','A',:,'SM'].JFL.mean(level=2)
# Weather        14.5D     2.9F     7.7D
# IS
# 024-01-02-L  12.0563  15.9398  12.5933
# 024-02-G-DA  14.3680  12.1194  13.8216
# 025-02-02-L  11.9469  15.7952  12.4791


# pJFLW = pJFW['JFL']*pJFW['JFF']#Jet length weighted by Frequency
# pJFLW.sum(level='IS')/pJFW.sum(level='IS').JFF #IS-wise weighted jet length
# pJL_Module = pJFLW.sum(level='Module')/pJFW.sum(level='Module').JFF #IS-wise weighted jet length
# pJL_Module.plot.bar()
# plt.show()

#Optoin 1 - Arithmetic average
# JetLengthModule = pJFW[('JFL','7.7D')].mean(level='Module') #IS-wise weighted jet length
# plttitle ='Average Jet Length per Module' 

#Optoin 2 - Frequency-weighted average over module
pJFLW77 = pJFW[('JFL','7.7D')]*pJFW[('JFF','7.7D')]
# JetLengthModule = pJFLW77.sum(level='Module')/pJFW[('JFF','7.7D')].sum(level='Module') #IS-wise weighted jet length
# plttitle ='Frequency-Weighted Jet Length averaged over Module' 

#Optoin 3 - Frequency-weighted average over Process
JetLengthModule = pJFLW77.sum(level='Module')/pJFW[('JFF','7.7D')].sum() #IS-wise weighted jet length
plttitle ='Frequency-Weighted Jet Length averaged over Process' 

JetLengthModule.plot.bar()
plt.xticks(rotation=45,ha='right')
plt.title(plttitle)
plt.tight_layout()
# plt.title('Frequency-Weighted Jet Length averaged over each Module')
plt.show()

#Duration analysis




#JT method 2 using automatic indexing, looks more efficient
listJF = []
#To validate the radiation distance model by CMPT: compare JetFire.D12 & D37 
for e in lEvent:
    if (e.JetFire != None):
        jfl = e.JetFire.Length        
        d37 = e.JetFire.D37
        d12 = e.JetFire.D12
        key,hole,weather=e.Key.split("\\")
        safety = hole[3:]    
        hole = hole[:2]        
        listJF.append([e.Module,e.Deck,key, hole, safety, weather,e.JetFire.Frequency,jfl,d12, 1.45*jfl,d37,1.2*jfl, e.JetFire.SEP, e.Discharge.Duration])
pdJF = pd.DataFrame(listJF,columns=['Module','Deck','Key','Hole','Safety','Weather','Frequency','Length','D12_PHAST','D12_CMPT','D37_PHAST','D37_CMPT','SEP','Duration'])
pdJFc = pdJF.drop(['D37_PHAST','D12_PHAST'],axis='columns')
pdJ = pdJFc.drop(['D12_CMPT','D37_CMPT'],axis='columns')
pdJm = pdJ.set_index(['Module','Deck','Key','Hole','Safety'])#Multi-index

# #Over the whole process area weighted average
# pdJ['LengthFWeighted'] = pdJ['Length']*pdJ['Frequency']
# aveLengthFPSO = pdJ['LengthFWeighted'].sum()/pdJ['Frequency'].sum()



#to convert into a multi-level index
# pdJFc_MD = pdJFc.set_index(['Module','Deck','Key','Hole','Weather'])
pdJFc_MD = pdJFc.set_index(['Module','Deck','Key'])
pdJFc_MD.Length.mean(level='Module').plot.bar()

pdJFc.groupby(['Module','Deck']).agg({'Frequency':'sum','Length':'mean'})

pdJFc.groupby(['Module','Deck']).Length.mean().plot.bar()
plt.title('Avg Jet Flame Length - Module & Deck')
pdJFc.groupby(['Module']).Length.mean().plot.bar()
plt.title('Avg Jet Flame Length - Module')
plt.xticks(rotation=45,ha='right')
plt.tight_layout()
plt.show()

pdJFc_MD.Length.mean(level='Deck').plot.bar()


is_s02 = pJF['Module'] == 'S02'
is_A  = pJF['Deck'] == 'A'
is_B  = pJF['Deck'] == 'B'
is_SM = pJF['Hole'] == 'SM'
is_ME = pJF['Hole'] == 'ME'
pJF[is_s02 & is_A]


pJFc[is_p05 & is_A & (is_SM | is_ME)].JFL.mean()



pdJF['LengthFWeighted'] = pdJF['Length']*pdJF['Frequency']
aveLengthFPSO = pdJF['LengthFWeighted'].sum()/pdJF['Frequency'].sum()

pdJF['DurationFWeighted'] = pdJF['Duration']*pdJF['Frequency']
aveDurationFPSO = pdJF['DurationFWeighted'].sum()/pdJF['Frequency'].sum()




# get_wavgp = lambda g: np.average(g['D37p'], weights = g['JFF'])
# pJF_w.apply(get_wavgp)
# get_wavgc = lambda g: np.average(g['D37c'], weights = g['JFF'])
# pJF_w.apply(get_wavgc)

# pJF_W = pJF.dropna().groupby('Weather')
#          PHAST     CMPT  Ratio
# 2.9F     60.4       66.1 9.4%
# 7.7D     59.6       64.1 7.6%
# 14.5D    68.4       73.1 6.9%
# pJF_M = pJF.dropna().groupby('Module')

# pJFc = pJF.drop(['D37p','D12p','D04p'],axis='columns')
# pJF_M = pJFc.groupby('Module')



pJFc = pJF.drop(['D37p','D12p','D04p'],axis='columns')
is_p05 = pJFc['Module'] == 'P05'
is_A  = pJFc['Deck'] == 'A'
is_B  = pJFc['Deck'] == 'B'
is_SM = pJFc['Hole'] == 'SM'
is_ME = pJFc['Hole'] == 'ME'
pJFc[is_p05 & is_A & (is_SM | is_ME)].JFL.mean()
#31.8m
is_P05 = pJFc['Module'] == 'P05'
pJFc[is_P05 & is_A].JFF.sum()

 pJFc.groupby('Module').agg({'JFF':'sum'})

#Sum of frequency whose initial jet length is longe than 15m from P05
is_longjet = pJFc['JFL'] > 15
pJFc[is_p05 & is_A & is_longjet].JFF.sum()
# 5.6E-6


#Frequency weighted average jet flame length
get_wavgp = lambda g: np.average(g['JFL'], weights = g['JFF'])
#Frequency-weighted flame length
pJF_M.apply(get_wavgp)
#Average frequency
pJF_M['JFF'].mean()

#
# pJF_M.agg({'JFF':'mean','JFL':'mean','JFL':get_wavgp})
pJF_M.agg({'JFF':['sum','mean'],'JFL':'mean'})
# Out[251]: 
#              JFF                      JFL       JFL
#              sum          mean       mean freq-weighted mean      
# Module                                    
# KOD     0.000004  3.728599e-08  11.030258  7.950744
# P02     0.000002  4.243934e-08  19.859478 19.256705
# P03     0.000052  2.052144e-07  23.400352 21.552296
# P04     0.000082  3.538242e-07  49.648065 50.235762
# P05     0.000058  3.211504e-07  50.614229 48.361007
# S02     0.000002  8.321024e-08  13.284389 11.980998
# S03     0.000053  3.273099e-07  31.139052 37.181575
# S04     0.000064  4.817382e-07  26.607697 31.494818
# S05     0.000126  3.976437e-07  42.091717 43.279454
pJF_M.apply(get_wavgp)

#            JFL
# Module   freq-weighted mean
# KOD      7.950744
# P02     19.256705
# P03     21.552296
# P04     50.235762
# P05     48.361007
# S02     11.980998
# S03     37.181575
# S04     31.494818
# S05     43.279454


pJF[['D37p','D37c']].dropna().mean()
pJF[['D12p','D12c']].dropna().mean()
# pJF[['D37p','D04c']].dropna().mean()

#Ratio weather dependency
pJF_w = pJF.dropna().groupby('Weather')
pJF_w.mean()

#Weighted dependency
pdJF = pJF.dropna()


pdJF['LengthFWeighted'] = pdJF['JFL']*pdJF['JFF']
aveLengthFPSO = pdJF['LengthFWeighted'].sum()/pdJF['Frequency'].sum()





pool_diameter = [0, 0.5, 2.5, 5, 10, 25, 30, 50, 80,100] #the lowest pool diamter and its radiation distance is added by SHK refering to ther ation at 2.5m
distance_12  = [0, 1.26, 6.3, 9.4, 11, 13.2, 15.6, 25.2, 40.2, 50]
distance_5 = [0, 2.4, 12, 19, 28, 30, 32, 43, 64, 77]
#Distance from the center of pool i.e. radius
f12 = interp1d(pool_diameter,distance_12)
f5 = interp1d(pool_diameter,distance_5)

EPF = {}
for e in lEvent:            
    pv,hole,weather = e.Key.split("\\")
    safety = hole[3:]    
    hole = hole[:2]            
    if e.EarlyPoolFire != None:
        epfd = e.EarlyPoolFire.Diameter
        # EPF[(e.Module,e.Deck,pv,hole,safety,weather)] = [e.EarlyPoolFire.Frequency,epfd,float(f12(epfd)),float(f5(epfd))]        
        EPF[(e.Module,e.Deck,pv,hole,safety,weather)] = [e.EarlyPoolFire.Frequency,epfd]        
        # EPF[(e.Module,e.Deck,e.Key)] = [pv,hole[:2],hole[3:],weather,e.EarlyPoolFire.Frequency,epfd,f12(epfd),f5(epfd)]        
        # EPF[e.Key] = [pv,hole[:2],hole[3:],weather,e.Module,e.Deck,e.EarlyPoolFire.Frequency,epfd,float(f12(epfd)),float(f5(epfd))]        
        # fmt_jf = "{:10.2e}\t{:8.1f}\t{:8.1f}\t{:8.1f}\t{:8.1f}\t{:8.1f}\t".format(e.JetFire.Frequency,e.JetFire.SEP,jfl,e.JetFire.D37,e.JetFire.D12,e.JetFire.D04)
#Analysis for Radiation distances from PHAST and CMPT formula
pEPF = pd.DataFrame(EPF)
pEPF = pEPF.T
# pEPF.columns=['EPFF','EPFD','D12','D05']
pEPF.columns=['EPFF','EPFD']
pEPF.index.names = ['Module','Deck','IS','Hole','Safety','Weather']


pEPFW = pEPF['EPFD']*pEPF['EPFF']
# pEPFDW = pEPFW.sum(level='IS')/pEPF.sum(level='IS').EPFF #IS-wise weighted jet length
#Option 2
pEPFDW = pEPFW.sum(level='Module')/pEPF.sum(level='Module').EPFF #IS-wise weighted jet length
plttitle = 'Frequency-Weighted Early Pool Fire Diameter per Module'

#Option 3
# pEPFDW = pEPFW.sum(level='Module')/pEPF.sum().EPFF #IS-wise weighted jet length
# plttitle = 'Frequency-Weighted Early Pool Fire Diameter averaged over Process'

pEPFDW.plot.bar()
plt.xticks(rotation=45,ha='right')
plt.title(plttitle)
plt.tight_layout()
plt.show()


pEPF.loc[pEPF.EPFF == 0].EPFD.mean(level='Module')


pEPF = pEPF.sort_index()
pEPF.loc['S04']
pEPF.loc['S04':'S05']
pEPF.loc['S04','B'].unstack()[[('EPFF','7.7D'),('EPFD','7.7D')]]



# pEPF.groupby('Module').agg({'EPFF':['sum','mean'],'EPFD':'mean'})
#                 EPFF                     EPFD
#                  sum          mean       mean
# Module
# KOD     1.062780e-05  3.321187e-07   4.358239
# P02     8.162569e-07  2.040642e-07   7.160602
# P03     3.527005e-05  1.469585e-06   2.599001
# P04     1.471484e-05  1.226236e-06   1.221014
# S02     1.199273e-05  4.441751e-07  12.279412
# S03     4.151587e-06  3.459656e-07   2.726317
# S04     3.655330e-04  5.538378e-06   6.110860
# S05     2.744747e-04  6.535111e-06   9.400246

pEPF.groupby(['Module','Deck']).agg({'EPFF':['sum','mean'],'EPFD':'mean'})
#                      EPFF                     EPFD
#                       sum          mean       mean
# Module Deck
# KOD    A     1.062780e-05  3.321187e-07   4.358239
# P02    B     8.162569e-07  2.040642e-07   7.160602
# P03    A     3.527005e-05  1.469585e-06   2.599001
# P04    B     1.471484e-05  1.226236e-06   1.221014
# S02    A     1.199273e-05  4.441751e-07  12.279412
# S03    A     4.151587e-06  3.459656e-07   2.726317
# S04    A     1.595012e-05  1.329177e-06   8.709768
#        B     3.495829e-04  6.473757e-06   5.533324 #S04, Pool fire more frequently at Deck B
# S05    A     2.744747e-04  6.535111e-06   9.400246 #S05, Pool fire at Deck A which is chosen as dimensioning



# pEPF.mean(level='Module')
pEPF.mean(level=['Module','Deck'])
pEPF.sum(level=['Module','Deck'])

#indexing using index
pEPF_s02 = pEPF.loc["S02"]
pEPF_s02[pEPF_s02.D12 > 10].EPFF.sum()
# Out[52]: 1.8903532164063199e-06

#with multilevel index
pEPF.loc[('S02','A')]


pEPF[pEPF['D12']>10].EPFF.sum()

pEPF[pEPF['EPFD']>10].EPFF.sum(level=['Module','Deck'])


listEPF = []
for e in lEvent:
    if (e.EarlyPoolFire != None):        
        key,hole,weather=e.Key.split("\\")
        listEPF.append([key, hole, weather,e.Module,e.Deck, e.EarlyPoolFire.Frequency,e.EarlyPoolFire.Diameter])

pdEPF = pd.DataFrame(listEPF,columns=['Key','Hole','Weather','Module','Deck','Frequency','Diameter'])
pdEPFGrouped = pdEPF.groupby(['Module','Deck'])
# pdEPFGrouped['Frequency'].sum()
pdEPFGrouped['Frequency'].sum().plot.bar(logy=True)
plt.title('Early pool fire - Frequency per Module-Deck [#/year]')
plt.tight_layout()
plt.show()

is_s05 = pdEPF["Module"] == "S05"
pdEPF[is_s05].groupby('Key').sum()
is_NZF = pdEPF["Frequency"] > 0.
pdEPF[is_s05 & is_NZF][['Key','Hole','Frequency','Diameter']]


pdEPF[is_s05][['Key','Hole','Frequency','Diameter']]


is_s04 = pdEPF["Module"] == "S04"
pdEPF[is_s04].groupby('Key').sum()
is_ZF = pdEPF["Frequency"] == 0.
pdEPF[is_s04 & is_NZF][['Key','Hole','Frequency','Diameter']]


# for m in ['P05','P04','P03','P02','S05','S04','S03','S02']:
#     print(m,EquipCoamingArea[m]*DrainRateModuleVol[m]/DrainRateModuleArea[m])




#Pool fire scenarios where pool is not possible due to drain
for e in lEvent:
    if (e.EarlyPoolFire != None) and (e.EarlyPoolFire.Frequency == 0):
        print(e.Key, e.EarlyPoolFire.Frequency, e.EarlyPoolFire.Diameter)


for e in lEvent:
    if (e.EarlyPoolFire != None) and (e.EarlyPoolFire.Frequency == 0) :
        print("{:30s} {:8.1f} {:10.2e} {:s}".format(e.Key, e.EarlyPoolFire.Diameter,e.EarlyPoolFire.Frequency, ' '.join(["{:8.1f}".format(v) for v in e.EarlyPoolFire.Ts[:3]])))


#Late pool fire
listLPF = []
for e in lEvent:
    if (e.EarlyPoolFire != None):        
        key,hole,weather=e.Key.split("\\")
        listLPF.append([key, hole, weather,e.Module,e.Deck, e.LatePoolFire.Frequency,e.LatePoolFire.Diameter])

pdLPF = pd.DataFrame(listLPF,columns=['Key','Hole','Weather','Module','Deck','Frequency','Diameter'])
pdLPFGrouped = pdLPF.groupby(['Module','Deck'])
# pdLPFGrouped['Frequency'].sum()
pdLPFGrouped['Frequency'].sum().plot.bar(logy=True)
plt.title('Late pool fire - Frequency per Module-Deck [#/year]')
plt.tight_layout()
plt.show()

is_s04 = pdLPF["Module"] == "S04"
is_ZF = pdLPF["Frequency"] == 0.
pdLPF[is_s04 & is_ZF][['Key','Hole','Frequency','Diameter']]


LPF = {}
for e in lEvent:            
    pv,hole,weather = e.Key.split("\\")
    safety = hole[3:]    
    hole = hole[:2]            
    if e.LatePoolFire != None:        
        LPF[(e.Module,e.Deck,pv,hole,safety,weather)] = [e.LatePoolFire.Frequency,e.LatePoolFire.Diameter]       
        
#Analysis for Radiation distances from PHAST and CMPT formula
pLPF = pd.DataFrame(LPF)
pLPF = pLPF.T
# pLPF.columns=['LPFF','LPFD','D12','D05']
pLPF.columns=['LPFF','LPFD']
pLPF.index.names = ['Module','Deck','IS','Hole','Safety','Weather']


pLPFW = pLPF['LPFD']*pLPF['LPFF']
# pLPFDW = pLPFW.sum(level='IS')/pLPF.sum(level='IS').LPFF #IS-wise weighted jet length
#Option 2
# pLPFDW = pLPFW.sum(level='Module')/pLPF.sum(level='Module').LPFF #IS-wise weighted jet length
# plttitle = 'Frequency-Weighted Late Pool Fire Diameter per Module'
#Option 3
pLPFDW = pLPFW.sum(level='Module')/pLPF.sum().LPFF #IS-wise weighted jet length
plttitle = 'Frequency-Weighted Late Pool Fire Diameter averaged over Process'

pLPFDW.plot.bar()
plt.xticks(rotation=45,ha='right')
plt.title(plttitle)
plt.tight_layout()
plt.show()

#Compare Early and Late pool fire frequencies
pEfrq = pEPF.sum(level=['Module','Deck']).EPFF
pLfrq = pLPF.sum(level=['Module','Deck']).LPFF
pPFF = pd.merge(pEfrq,pLfrq,left_index=True, right_index=True)
pPFF.plot.bar(logy=True)
plt.title('Sum of Early/Late Pool Fire Frequencies [#/year]')
plt.xticks(rotation=45,ha='right')
plt.tight_layout()
plt.show()


pdEPF.groupby('Module').Frequency.sum() + pJF.groupby('Module').JFF.sum()
pJF.sum(level='Module').JFF
pJF.groupby('Module').JFF.sum()
pdEPF.groupby('Module').Frequency.sum()
