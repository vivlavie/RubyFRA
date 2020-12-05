import numpy as np
import math
#PyAnalysis1.py
#for Event analyses
#Process Ara
Area = 'ProcessArea'
# element_dump_filename = 'Bv06_dump_PyExdCrv6'
element_dump_filename = 'Bv06_dump'
element_dump_filename = 'Bv06_dump_PyExdCrv6'


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
with open(element_dump_filename,'rb') as element_dump:
    lEvent = dill.load(element_dump)

F = 0.
for e in lEvent:    
    F += e.Frequency
print("{:20s} Total Leak Frequency {:8.2e}".format(Area,F))

Modules = []
for e in lEvent:    
    if e.Module not in Modules:
        Modules.append(e.Module)



#Leak frequency per release rate range
Module = "P02"
Fl = Fl0 = Fl1 = Fl2 = Fl3 = Fl4 = 0.
for e in lEvent:  
    # if (Module in e.Module) and ('-G' in e.Key):
    if (Module in e.Module) and ('046-' in e.Key):
        # if (e.Discharge != None) and (e.Discharge.ReleaseRate > 1.0):
        f = e.Frequency
        if (e.Discharge.ReleaseRate <= 1.0):
            Fl0 += f
        elif (e.Discharge.ReleaseRate <= 10.0):
            Fl1 += f
        elif (e.Discharge.ReleaseRate <= 50.0):
            Fl2 += f
        elif (e.Discharge.ReleaseRate <= 150.0):
            Fl3 += f
        else:
            Fl4 += f
        Fl += f
print("{:15s} Release Frequency: {:.2e} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e}".format(Area, Fl0,Fl1,Fl2,Fl3,Fl4, Fl))

for e in lEvent:  
    # if ('ME_EOBO' in e.Hole) and ('020-03-01-Li' in e.Key) and ('2.9F' in e.Weather):
    if ('LM_EOBO' in e.Hole) and ('013-01' in e.Key) and ('7.7DF' in e.Weather):
        print(e.Key, e.Hole, e.Weather)
        print(e.Dispersion.DfMax, DMax)

for e in lEvent:  
    if (Module in e.Module) and not ('046-' in e.Key):
        print(e.Key)

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

#Distance of radiation extent to 6.3 kW/m2
from scipy.interpolate import interp1d


#Jet fire frequency
Fj = FjA = FjB = FjC = FjB_Long = FjB_Short = JL = DL = DLA = 0.
for e in lEvent:  
    # if ('S05' in e.Module) and (e.Deck == 'C'):
    if ('S05' in e.Module):
        if (e.JetFire != None):

            tvd = np.array(e.TVD)
            T = tvd[:,0] #t
            m = tvd[:,2] #rr
            J = np.zeros(m.shape)
            for i in range(len(m)):
                J[i] = e.jfscale*jffit(m[i])
            Tfit = interp1d(J,T)
            if (J[-1] > 18):
                Dfit = T[-1]
            else:
                if (J[0] > 18):
                    Dfit = float(Tfit(18))
                else:
                    Dfit = 0.
            if (J[-1] > 4.5):
                DfitA = T[-1]
            else:
                if (J[0] > 4.5):
                    DfitA = float(Tfit(4.5))
                else:
                    DfitA = 0.
            f = e.JetFire.Frequency
            Fj += f
            if e.Deck == "A":
                FjA += f
                DLA += f*DfitA
            elif e.Deck == "B":
                FjB += f
                JL += f*e.JetFire.Length
                DL += f*Dfit
                if e.JetFire.JetLengths[1] > 18:
                    FjB_Long += f
                else:
                    FjB_Short += f
            elif e.Deck == "C": # No Deck C at S05
                FjC += f

print("{:15s} Jet Fire Frequency: {:.2e} {:.2e} {:.2e} (={:.2e}+{:.2e}) {:.2e} Weighted JL(t=0) {:.2e} Weighted Duration Deck A {:.2e}".format(Area, Fj, FjA, FjB, FjB_Long, FjB_Short, FjC, JL/FjB, DLA/FjA))

#Weighted exposure duration from both Deck A and B
(DLA+DL)/(FjB_Long+FjA)

#effect distance as ratio to thet jet length as per 'kW/m2'
# w_rad = {5:1.75, 12.5:1.45, 37.5:1.2}



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

#Fires frequency
Area = "KOD"
for m in Modules:
    Area = m
    Fj = Fep = Flp = 0.
    Dep = Dlp = 0.
    for e in lEvent:
        if Area in e.Module:
            # if "-L" in e.Key:
            #     F += e.Frequency

            if e.JetFire != None:
                Fj += e.JetFire.Frequency

            if e.EarlyPoolFire != None:
                Fep += e.EarlyPoolFire.Frequency
                Dep += PoolFireDuration(e,'Early')*e.EarlyPoolFire.Frequency
            if e.LatePoolFire != None:
                Flp += e.LatePoolFire.Frequency
                Dlp += PoolFireDuration(e,'Late')*e.LatePoolFire.Frequency
    print("{:10s} JF_Direc {:10.2e} PF_Early {:10.2e} PF_Late {:10.2e}".format(m,Fj,Fep,Flp,6*Fj+Fep+Flp))

print("{:15s} Jet Fire Frequency: {:.2e} sec".format(Area, Fj))
print("{:15s} Early Pool Fire Frequency: {:.2e} {:8.1f}sec".format(Area, Fep,Dep/Fep))
print("{:15s} Late Pool Fire Frequency: {:.2e} {:8.1f}sec".format(Area, Flp,Dlp/Flp))

#The total jet (omni-directional) and pool fire frequency:
6*Fj + Fep + Flp


# print("{:15s} Liquid Leak Frequency: {:.2e}".format(Area, F))
# print("{:15s} Pool fire Toal Frequency: Early {:.2e} Late {:.2e}".format(Area, Fe,Fl))
# print("{:15s} Pool fire Toal Frequency: {:.2e}".format(Area, sum(EPFreq.values()) + sum(LPFreq.values())))

# F = Fe = Fl = 0.
# for e in lEvent:
#     if e.EarlyPoolFire is not None:    
#         Fe += e.EarlyPoolFire.Frequency
#     if e.LatePoolFire is not None:    
#         Fl += e.LatePoolFire.Frequency

#Explosion
Fe = 0.
for e in lEvent:
    if Area in e.Module:
        if e.Explosion != None:
            Fe += e.Explosion.Frequency
print("{:15s} Explosion Frequency: {:.2e}".format(Area, Fe))





Fe = Fne = 0.
for e in lEvent:
    if e.Explosion != None:
        print("{:30s}{:10.2e}{:10.2e}{:10.2e}".format(e.Key, e.Frequency, e.PDelIgn*e.PExp_Ign, e.Explosion.Frequency))
        Fe += e.Explosion.Frequency
    else:
        print("{:30s}{:10.2e}{:10.2e} No explosion".format(e.Key, e.Frequency, e.PDelIgn*e.PExp_Ign))
        Fne += e.Frequency
print ("{:10.2e}{:10.2e}{:10.2e}".format(Fe, Fne, Fe+Fne))

F = 0.
Fj = 0.
JFC = 0
D = 0.
Area = 'S01'
Fj2 = 0.
for e in lEvent:  
    if Area in e.Module:
        F += e.Frequency
        Fj2 += e.JetFire.Frequency
        if (e.Discharge != None) and (e.JetFire.Length > 8):
            Fj += e.JetFire.Frequency
            D += e.TVD[-1][0]
            JFC += 1
print("Leak frequency: {:10.2e}".format(F))
print("{:15s} Jet Fire Frequency affecting Workshop: {:.2e}, Average duration: {:.2e}".format(Area, Fj, D/JFC))





#effect distance as ratio to thet jet length as per 'kW/m2'
# w_rad = {5:1.75, 12.5:1.45, 37.5:1.2}
j_rad = [5, 12.5, 27.5]
r_dist = [1.75, 1.45, 1.2]

#Distance of radiation extent to 6.3 kW/m2
from scipy.interpolate import interp1d
f_r_dist = interp1d(j_rad,r_dist)
d_6_3 = float(f_r_dist(6.3))

for mm in ['03', '04', '05']:
    Fj = 0.
    for e in lEvent:  
        if mm in e.Module:
            if d_6_3*e.JetFire.Length > 50:
                Fj += e.JetFire.Frequency
    print("{:15s} Jet Fire Frequency whose initial length is longer than 25m: {:.2e}".format(mm, Fj))




#Large leak of heavy gas frequency analysis
F = 0.
HeavyGasMax = 0.
HeavyGasModules = ['S04','S05']
for e in lEvent:    
    HeavyGasRR = e.Discharge.ReleaseRate * e.Discharge.LiquidMassFraction
    HeavyGasMax = max(HeavyGasMax,HeavyGasRR)
    if ("-L" in e.Key) and (e.Module in HeavyGasModules) and ( HeavyGasRR > 40):
    # if ("-L" in e.Key) and (e.Discharge.ReleaseRate * e.Discharge.LiquidMassFraction > 40):
        print(e.Key, e.Discharge.ReleaseRate)
        F += e.Frequency


#Fires frequency

for m in Modules:
    Area = m
    Fj = Fep = Flp = 0.
    Dep = Dlp = 0.
    mass_release = tvd[:,1] #t
    for e in lEvent:
        if Area in e.Module:
            # if "-L" in e.Key:
            #     F += e.Frequency
            if (e.EarlyPoolFire != None) and ( (mass_release[0]-mass_release[-1]) > 100 ) and :
                Fep += e.EarlyPoolFire.Frequency
                # Dep += PoolFireDuration(e,'Early')*e.EarlyPoolFire.Frequency
            if (e.LatePoolFire != None) and ( (mass_release[0]-mass_release[-1]) > 100 ):
                Flp += e.LatePoolFire.Frequency
                # Dlp += PoolFireDuration(e,'Late')*e.LatePoolFire.Frequency
    print("{:10s} PF_Early {:10.2e} PF_Late {:10.2e}".format(m,Fep,Flp,6*Fj+Fep+Flp))

 for e in lEvent:
    if ("-L" in e.Key):
        if (e.EarlyPoolFire != None):
            if (e.EarlyPoolFire.SEP > 100):
                print(e.Key, e.Module, e.EarlyPoolFire.Diameter)

print(HeavyGasMax)



# fn = 'C\\ProcessArea\\C_'+slugify(Key)+'.png'

fn = 'MAH_jet_v01.png'
fig.savefig(fn) 
plt.close()


#Plot a fire

pool_diameter = [0.5, 2.5, 5, 10, 25, 30, 50, 80,100] #the lowest pool diamter and its radiation distance is added by SHK refering to ther ation at 2.5m
distance_12  = [1.26, 6.3, 9.4, 11, 13.2, 15.6, 25.2, 40.2, 50]
distance_5 = [2.4, 12, 19, 28, 30, 32, 43, 64, 77]
#Distance from the center of pool i.e. radius
f12 = interp1d(pool_diameter,distance_12)
f5 = interp1d(pool_diameter,distance_5)

for e in lEvent:    
# for e in [lEvent[5]]:    
    # if ("-L" in e.Key) and ("LM_EO" in e.Hole) and (e.Weather == "7.7D"):    
        Key = e.Key
        X = e.X
        Y = e.Y        
        
        if e.Deck == 'A':
            img = plt.imread("RubyFPSODeckA.jpg")
        elif (e.Deck == 'B') or (e.Deck == 'C'):
            img = plt.imread("RubyFPSODeckB.jpg")
        elif e.Deck == 'Hull':
            # print('Hull deck?')
            img = plt.imread("RubyFPSOHullDeck.jpg")
        
        if (e.EarlyPoolFire != None) or (e.LatePoolFire != None):
            fig,ax = plt.subplots(3,1,figsize=(8,7.5))
            for axid in [0,1,2]:
                ax[axid].imshow(img, extent=[-11, 252, -32.43, 50.19])
        else:
            fig,ax = plt.subplots(2,1,figsize=(8,5))
            for axid in [0,1]:
                ax[axid].imshow(img, extent=[-11, 252, -32.43, 50.19])

        #Plot Dispersion - LFL
        if (e.Dispersion != None) and (e.Dispersion.DMin != ' '):
            DMin = e.Dispersion.DMin
            DMax = e.Dispersion.DMax
            W = e.Dispersion.W    
            AxisLong = (DMax - DMin)
            AxisShort = W
            DX = (DMax + DMin)*0.5            
            cloud = Ellipse((X-XDir*DX,Y),AxisLong,AxisShort,0)
            # ax.add_artist(cloud)
            ax[0].add_patch(cloud)
            cloud.set_alpha(0.9)
            cloud.set_facecolor('Magenta')
            cloud.set(label='LFL')
        
        #Plot Dispersion - Flash Fire
        if (e.Dispersion != None) and (e.Dispersion.FFMaxDistance != ' '):        
            FFMax = e.Dispersion.FFMaxDistance
            FFW = e.Dispersion.FFWidth
            AxisLong = FFMax
            AxisShort = FFW
            
            cloud = Ellipse((X-XDir*AxisLong*0.5,Y),AxisLong,AxisShort,0)
            # ax.add_artist(cloud)
            ax[0].add_patch(cloud)
            cloud.set_alpha(0.3)
            cloud.set_facecolor('Yellow')
            cloud.set(label='Flash (50\% LFL)')
            ax[0].legend()
        ax[0].set_title(Key)
        


        #Plot Jet Fire
        if e.JetFire != None:
            jfl = e.JetFire.Length
            jfw = 0.12*jfl

            #5, 12, 37.5 kW/m2 contour     
            alpha = 0.8*0.8*0.8     
            w_rad = {5:1.75, 12.5:1.45, 37.5:1.2}
            for rad in [5, 12.5, 37.5]:
                # jfl_rad = np.sqrt(e.JetFire.SEP/rad)+jfl
                jfl_rad = jfl * w_rad[rad]
                jfw_rad = 0.12*jfl_rad

                #Radiation as an ellipse ... the max width is not at the flame tip but at the mid of the flame length
                cloud = Ellipse((X-XDir*jfl_rad*0.5,Y),jfl_rad,jfw_rad,0)                
                ax[1].add_patch(cloud)
                cloud.set_alpha(alpha)
                alpha /= 0.8
                cloud.set_facecolor('Green')
                cloud.set(label="Jet-{:4.1f}kW/m2".format(rad))                
            
            #jet fire shape
            JF = np.array([[X,Y],[X-XDir*jfl,Y+0.5*jfw],[X-XDir*jfl,Y-0.5*jfw]])
            jfpatch = plt.Polygon(JF,color='red')
            ax[1].add_patch(jfpatch)
            jfpatch.set(label="Jet({:3.0f}kW/m2)".format(e.JetFire.SEP))
            ax[1].legend()

        #Plot Early Pool Fire        
        if e.EarlyPoolFire != None:
            EPFD = e.EarlyPoolFire.Diameter
            if e.LatePoolFire != None:            
                LPFD = e.LatePoolFire.Diameter
                PFD = max(EPFD,LPFD)
                SEP = max(e.EarlyPoolFire.SEP,e.LatePoolFire.SEP)
            else:
                PFD = EPFD
                SEP = e.EarlyPoolFire.SEP
            #5, 12, 37.5 kW/m2 contour     
            alpha = 0.8*0.8*0.8     
            for rad in [5, 12.5, 37.5]:
                
                #Total length i.e. diameter is to be calculated
                # D_rad = 2*(np.sqrt(SEP/rad)+PFD/2)
                if PFD > 100:
                    if rad == 37.5:
                        D_rad = 1.*PFD 
                    elif rad == 12.5:
                        D_rad = 1.*PFD
                    elif rad == 5:
                        D_rad = 77/50*PFD                    
                else:
                    if rad == 37.5:
                        D_rad = 1.*PFD 
                    elif rad == 12.5:
                        D_rad = 2*f12(PFD)
                    elif rad == 5:
                        D_rad = 2*f5(PFD)
                # print (rad,D_rad)

                #Radiation as an ellipse ... the max width is not at the flame tip but at the mid of the flame length                
                cloud = Ellipse((X,Y),D_rad,D_rad,0)
                ax[2].add_patch(cloud)
                cloud.set_alpha(alpha)
                alpha /= 0.8
                cloud.set_facecolor('blue')
                cloud.set(label="Pool-{:4.1f}kW/m2".format(rad))                

        if e.LatePoolFire != None:            
            LPF = Ellipse((X,Y),LPFD,LPFD,0)            
            ax[2].add_patch(LPF)
            LPF.set_alpha(1)
            LPF.set_facecolor('cyan')
            LPF.set(label="'Late Pool-{:4.1f}kW/m2".format(e.LatePoolFire.SEP))
        if e.EarlyPoolFire != None:            
            EPF = Ellipse((X,Y),EPFD,EPFD,0)            
            ax[2].add_patch(EPF)
            EPF.set_alpha(0.5)
            EPF.set_facecolor('blue')
            EPF.set(label="'Early Pool-{:4.1f}kW/m2".format(e.EarlyPoolFire.SEP))
        if (e.EarlyPoolFire != None) or (e.LatePoolFire != None):
            ax[2].legend()
        
        
        # ax.clabel(cs2,fmt='%.1e',colors='k',fontsize=14)
        # ax.set_aspect('equal')
        # ax.xaxis.set_major_locator(plt.FixedLocator([120, 141, 168, 193]))
        # ax.xaxis.set_major_formatter(plt.FixedFormatter(['2/3','3/4','4/5','S05']))
        # ax.yaxis.set_major_locator(plt.FixedLocator([-27, -3.1, 3.1, 27]))
        # ax.yaxis.set_major_formatter(plt.FixedFormatter(['ER_S','Tray_S','Tray_P','ER_P']))

        
        # fn = 'C\\ProcessArea\\C_'+slugify(Key)+'.png'
        fn = 'C\\'+folder+'\\C_'+slugify(Key)+'_2.png'

        fig.savefig(fn) 
        plt.show()
        plt.close()


for m in Modules:
    Area = m
    Fj = Fep = Flp = 0.
    Dep = Dlp = 0.
    mass_release = tvd[:,1] #t
    for e in lEvent:
        if Area in e.Module:


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




listJF = []
#To validate the radiation distance model by CMPT: compare JetFire.D12 & D37 
for e in lEvent:
    if (e.JetFire != None):
        jfl = e.JetFire.Length        
        d37 = e.JetFire.D37
        d12 = e.JetFire.D12
        key,hole,weather=e.Key.split("\\")
        listJF.append([key, hole, weather,e.Module,e.JetFire.Frequency,jfl,d12, 1.45*jfl,d37,1.2*jfl])
        if (d37 != 'n/a') and (d12 != 'n/a'):
            r37 = d37/jfl/1.2
            r12 = d12/jfl/1.45
            fmt_jf = "{:30s}\t{:8.1f}\t{:8.1f}\t{:8.1f}\t{:8.1f}%\t{:8.1f}\t{:8.1f}\t{:8.1f}%".format(e.Key,jfl,e.JetFire.D12, 1.45*jfl,r12*100,e.JetFire.D37,1.2*jfl,r37*100)
            print(fmt_jf)
        else:
            print("{:30s}\t{:8.1f} - D37 is not calculated by PHAST".format(e.Key,jfl))


pdJF = pd.DataFrame(listJF,columns=['Key','Hole','Weather','Module','Frequency','Length','D12_PHAST','D12_CMPT','D37_PHAST','D37_CMPT'])
pdJF['LengthFWeighted'] = pdJF['Length']*pdJF['Frequency']
aveLengthFPSO = pdJF['LengthFWeighted'].sum()/pdJF['Frequency'].sum()

import matplotlib.pyplot as plt


pdJFM = pdJF.groupby('Module')
# aveLength = pdJFM['LengthFWeighted'].sum()/pdJFM['Frequency'].sum() #manual way of weighted average
get_wavg = lambda g: np.average(g['Length'],weights=g['Frequency'])
pdJFM.apply(get_wavg).plot.bar()
plt.show()

#Module-wise frequency-weighted jet fire length & duration


listEPF = []
#To validate the radiation distance model by CMPT: compare JetFire.D12 & D37 
for e in lEvent:
    if (e.EarlyPoolFire != None):        
        key,hole,weather=e.Key.split("\\")
        listEPF.append([key, hole, weather,e.Module,e.Deck, e.EarlyPoolFire.Frequency,e.EarlyPoolFire.Diameter])

pdEPF = pd.DataFrame(listEPF,columns=['Key','Hole','Weather','Module','Deck','Frequency','Diameter'])
pdEPFGrouped = pdEPF.groupby(['Module','Deck'])
# pdEPFGrouped['Frequency'].sum()
pdEPFGrouped['Frequency'].sum().plot.bar(logy=True)
plt.show()


listLPF = []
#To validate the radiation distance model by CMPT: compare JetFire.D12 & D37 
for e in lEvent:
    if (e.LatePoolFire != None):        
        key,hole,weather=e.Key.split("\\")
        listLPF.append([key, hole, weather,e.Module,e.Deck, e.LatePoolFire.Frequency,e.LatePoolFire.Diameter])

pdLPF = pd.DataFrame(listLPF,columns=['Key','Hole','Weather','Module','Deck','Frequency','Diameter'])
pdLPFGrouped = pdLPF.groupby(['Module','Deck'])
pdLPFGrouped['Frequency'].sum()
pdLPFGrouped['Frequency'].sum().plot.bar(logy=True)
plt.show()

#To check if none of liquid release results in any pool fire
#The results is that is the case for S04, S03, & P04
RelLiqFreq = dict.fromkeys(Modules,0.)
NPFF = dict.fromkeys(Modules,0.) #No pool fire frequency
print("Liquid discharge but no pool fire")
for e in lEvent:
    if ("-L" in e.Key) and (e.Discharge.LiquidMassFraction > 0.):
        RelLiqFreq[e.Module] += e.Frequency
    # if ("-L" in e.Key):
    # if ("020-03-03" in e.Key):
    if ("-L" in e.Key) and (e.Discharge.LiquidMassFraction > 0.) and (e.EarlyPoolFire == None):
        print("{:40s} {:10s} {:8.1f}".format(e.Key, e.Module, e.Discharge.LiquidMassFraction))     
        NPFF[e.Module] += e.Frequency
        # break
for m in Modules:
    print(m,     RelLiqFreq[m] - NPFF[m])    




JF_Cube=['025-02-01-G\LM_EOBO\7.7DFO_Jet', '024-02-G-DC\LM_EOBO\7.7DFO_Jet','024-03-01-G\LM_EOBO\14.5DFO_Jet','020-02-01-Li\ME_EOBN\7.7DFO_Jet','025-01-01-G\LM_EOBO\7.7DFO_Jet',
'020-03-01-Li\LM_EOBO\7.7DFO_Jet','020-05-01-G\MA_EXBX\14.5DFO_Jet','020-02-01-Li\ME_EOBN\7.7DFO_Jet','027-01-G-DA\LM_EXBO\7.7DFO_Jet','027-01-G-DB\ME_EOBO\2.9FFO_Jet',
'027-01-G-DB\LM_EOBO\7.7DFO_Jet','024-02-G-DA\LM_EOBO\7.7DFO_Jet','020-02-01-Li\ME_EOBN\7.7DFO_Jet','027-01-G-DA\ME_EOBO\7.7DFO_Jet','023-03-01-G-DB\LM_EOBO\2.9FFO_Jet',
'013-01-C2-G\LM_EOBO\7.7DFO_Jet','020-05-01-G\LA_EXBO\2.9FFO_Jet','027-01-G-DB\LM_EOBO\7.7DFO_Jet','023-02-01-G-DB\LM_EOBO\14.5DFO_Jet',
'045-01-G\MA_EOBX\7.7DFO_Jet','045-01-G\ME_EOBO\2.9FFO_Jet','024-02-G-DA\LM_EOBX\14.5DFO_Jet','027-03-G\ME_EXBO\2.9FFO_Jet','025-02-01-G\LM_EOBO\7.7DFO_Jet',
'013-01-N2-G\LM_EOBO\7.7DFO_Jet','027-02-G\LM_EOBO\7.7DFO_Jet','023-01-01-G-DB\LM_EOBO\2.9FFO_Jet','045-01-G\MA_EOBX\7.7DFO_Jet']
for j in JF_Cube:
    scn,hole_weather = j.split("\\")
    for e in lEvent:
        if (scn in e.Key) and (hole_weather[:7] in e.Key):
            break
    print("{:20s}{:10s}{:8.1f}".format(scn, hole_weather[:7], e.TVD[2,0]))



#To check effect of vapor fraction of 5% on the ignition probability

IgnProb = {}
IgnProb['023-01']=[0.00050007132890309,0.00120241794400881,0.0121099047419176,0.0133941985161508]
IgnProb['023-02']=[0.000500009730687852,0.00112728346839808,0.0114299096696923,0.0124894860472479]
IgnProb['045-01']=[0.000500092602037187,0.00214457708972392,0.0159824850536888,0.015707932976315]
IgnProb['024-01']=[0.000500049238017383,0.00135755735655058,0.0127160815903571,0.0137890702655233]
IgnProb['024-02']=[0.000500092602037187,0.00177918758871906,0.014582259955286,0.0152682624925632]
IgnProb['024-03']=[0.000500092602037187,0.00187273284139629,0.0148864999002958,0.0153843297225681]
IgnProb['025-01']=[0.000500066037325877,0.00155843931004843,0.0138293983168883,0.0148186113336855]
IgnProb['025-02']=[0.000500076839107057,0.001583080559352,0.0138624647673917,0.0148213944403002]
IgnProb['027-01']=[0.000500092602037187,0.00141234355179696,0.0130284920753651,0.014098525976969]
IgnProb['027-02']=[0.000500092602037187,0.00167828753177232,0.0139540201469402,0.0148124734230255]
IgnProb['023-03']=[0.000500074056674887,0.00134865578738282,0.0127446106070196,0.013756233693594]
IgnProb['023-04']=[0.000500014189107731,0.00114561328550342,0.0121656333278765,0.0136982452845986]
IgnProb['024-01']=[0.000500056109924782,0.00140883072048731,0.013214286058242,0.0148872691410471]
IgnProb['045-02']=[0.000500077778731805,0.00141477302560298,0.012839767974552,0.0139304583821647]
IgnProb['020-01']=[0.000500081597970954,0.00198507390262245,0.0159033776172752,0.0159730079240214]
IgnProb['020-02']=[0.000500044233038885,0.00127252532590165,0.0126086841965461,0.0139360824741458]
IgnProb['013-01']=[0.000500081597970954,0.00125827417928573,0.0124609025752969,0.0136146062289823]
IgnProb['020-03']=[0.000500037598733664,0.00117389329582008,0.0122180292562774,0.0135140045120202]
IgnProb['020-04']=[0.000500009985715746,0.00113070577630583,0.0120388416453966,0.0131069382140561]
IgnProb['027-03']=[0.000500087509063428,0.00137082995866607,0.0129208012390835,0.0142136300044506]
IgnProb['General']=[0.000502,0.001677,0.017112,0.017916]


IgnitionModel = "OLF"
p = []

for e in lEvent:    
    # study,pv,pv_hole,hole,weather = e.Key.split("\\")
    pv,hole,weather = e.Key.split("\\")        
    key_hole = pv +"-"+ hole
    rr0 = e.Discharge.ReleaseRate         
    
    # Adjustment for liquid mass fraction
    if "-L" in key_hole:
        rr = rr0 * (1-min(e.Discharge.LiquidMassFraction,0.95))        #Total rate considerin the addtion vapour 5% only for the full liquid release
        rr2 = rr0 * (1 - e.Discharge.LiquidMassFraction + 0.05)         #Total rate x Vapour fraction with additional 5%
        if (rr2 > 1 and rr < 1) or (rr2 > 10 and rr < 10)  or (rr2 > 30 and rr < 30):
            print(key_hole, rr0, rr, rr2)

        if rr < 1:
            e.PImdIgn = 0.0005
            e.PDelIgn = p[0] - e.PImdIgn
            # e.PDelIgn = p[0]
        elif rr < 10:
            e.PImdIgn = 0.001
            e.PDelIgn = p[1] - e.PImdIgn        
            # e.PDelIgn = p[1]
        elif rr < 30.:
            e.PImdIgn = 0.01        
            e.PDelIgn = p[2] - e.PImdIgn        
            # e.PDelIgn = p[2]
        elif rr >= 30.:        
            e.PImdIgn = 0.01
            e.PDelIgn = p[3] - e.PImdIgn        

PImdIgnMax = 0.
PDelIgnMax = 0.
for e in lEvent:
    PImdIgnMax = max(PImdIgnMax, e.PImdIgn)
    PDelIgnMax = max(PDelIgnMax, e.PDelIgn)


for (e1,e2) in zip(lEvent, lEvent2):
    if e1.Key == e2.Key:
        if e1.PImdIgn == e2.PImdIgn:
            continue
        else:
            print(e1.Key, e2.PImdIgn, e2.PImdIgn)
        if e1.PDelIgn == e2.PDelIgn:
            continue
        else:
            print(e1.Key, e2.PDelIgn, e2.PDelIgn)

for e in lEvent:
    if ('020-03-01-L' in e.Key):
        if (e.EarlyPoolFire != None):
            print(e.Key, e.EarlyPoolFire.Duration)

f_max = 0.
for k in Frequency:
    if f_max < max(Frequency[k]):
        f_max = max(Frequency[k])
        is_max = k
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



pdF_L = pdF.index.str.contains("131-01")
pdF_L = pdF[pdF_L]
pdF_LM = pdF_L.groupby(Modules)
pdF_LM.SM.max()
pdF_LM.SM.sum()
pdF_LM.SM.describe()


# pdF = pd.DataFrame([Frequency,Modules,Deck])
# pdF.set_index(['SM','ME','ME','LM','Module','Deck'])


print_a_event(lEvent,'020-03-01-Li','SM_EOBO','2.9F')
e_MP = e_sm

print_a_event(lEvent,'020-04-01-L','SM_EOBN','2.9F')
# 342 020-04-01-L\SM_EOBN\2.9F
e_LP = lEvent[342]

print_a_event(lEvent,'020-02-01-Li','SM_EOBN','2.9F')
# 192 020-02-01-Li\SM_EOBN\2.9F
e_HP = lEvent[192]


for e in lEvent:     
    if e.Weather == "2.9F" and "EOBX" in e.Hole:
        pv,scenario,weather = e.Key.split('\\')
        hole,safety = scenario.split('_')
        dd1 = e.Discharge.Duration
        dd2 = e.TVD[-1,0]
        mm2 = e.TVD[-1,1]
        rr2 = e.TVD[-2,2]        
        mm0 = e.TVD[0,1]
        rr0 = e.TVD[0,2]
        dd0 = mm0/rr0
        # print("{:30s} {:8.1f} {:8.1f} {:8.0f}%".format(e.Key,dd1,dd0,dd0/dd1*100))
        # mm0/rr0
        print("{:30s} {:8.1f} {:8.1f} {:8.0f}% {:8.1f} {:8.1f}".format(e.Key,dd1,dd2,dd2/dd1*100,mm2,rr2))


JF = {}
JF2 = {}
for e in lEvent:            
    pv,hole,weather = e.Key.split("\\")
    safety = hole[3:]    
    hole = hole[:2]
    if e.JetFire != None:
        jfl = e.JetFire.Length
        JF[e.Key] = [pv,hole,safety,weather,e.Module,e.Deck,e.JetFire.Frequency,jfl,e.JetFire.D37,1.2*jfl, e.JetFire.D12,1.45*jfl,e.JetFire.D04,1.75*jfl]        
        JF2[(e.Module,e.Deck,e.Key)] = [e.JetFire.Frequency,jfl,1.2*jfl, 1.45*jfl,1.75*jfl]        
        # fmt_jf = "{:10.2e}\t{:8.1f}\t{:8.1f}\t{:8.1f}\t{:8.1f}\t{:8.1f}\t".format(e.JetFire.Frequency,e.JetFire.SEP,jfl,e.JetFire.D37,e.JetFire.D12,e.JetFire.D04)
#Analysis for Radiation distances from PHAST and CMPT formula
pJF = pd.DataFrame(JF)
pJF = pJF.T
pJF.columns=['IS','Hole','Safety','Weather','Module','Deck','JFF','JFL','D37p','D37c','D12p','D12c','D04p','D05c']
pJF = pJF.replace('n/a',np.nan)

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
    if e.EarlyPoolFire != None:
        epfd = e.EarlyPoolFire.Diameter
        EPF[(e.Module,e.Deck,e.Key)] = [e.EarlyPoolFire.Frequency,epfd,float(f12(epfd)),float(f5(epfd))]        
        # EPF[(e.Module,e.Deck,e.Key)] = [pv,hole[:2],hole[3:],weather,e.EarlyPoolFire.Frequency,epfd,f12(epfd),f5(epfd)]        
        # EPF[e.Key] = [pv,hole[:2],hole[3:],weather,e.Module,e.Deck,e.EarlyPoolFire.Frequency,epfd,float(f12(epfd)),float(f5(epfd))]        
        # fmt_jf = "{:10.2e}\t{:8.1f}\t{:8.1f}\t{:8.1f}\t{:8.1f}\t{:8.1f}\t".format(e.JetFire.Frequency,e.JetFire.SEP,jfl,e.JetFire.D37,e.JetFire.D12,e.JetFire.D04)
#Analysis for Radiation distances from PHAST and CMPT formula
pEPF = pd.DataFrame(EPF)
pEPF = pEPF.T
pEPF.columns=['EPFF','EPFD','D12','D05']
pEPF.index.names = ['Module','Deck','Key']

pEPF.groupby('Module').agg({'EPFF':['sum','mean'],'EPFD':'mean'})
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

#S02


DC = {}
for e in lEvent:            
    pv,hole,weather = e.Key.split("\\")
    safety = hole[3:]    
    hole = hole[:2]
    if e.JetFire != None:
        rr = e.Discharge.ReleaseRate
        dd = e.Discharge.Duration
        DC[(e.Module,e.Deck,pv,hole,safety,weather)] = [e.Frequency, rr, e.Discharge.LiquidMassFraction,dd]        
        
#Analysis for Radiation distances from PHAST and CMPT formula
pDC = pd.DataFrame(DC)
pDC = pDC.T
pDC.columns=['Frequency','ReleaseRate','LiqFrac','Duration']
pDC.index.names = ['Module','Deck','IS','Hole','Safety','Weather']

# c1 = pDC['Module'] == 'S03'
# c2 = pDC['Module'] == 'P03'

pDC['VapRelRate'] = pDC['ReleaseRate'] * (1-pDC['LiqFrac'])
c3 = pDC['VapRelRate'] > 40

pDC.loc[['S03','P03']][c3].sum()

#Frequency of havign flammable gas from S03 and P03 to GTG exhaust
pDC.loc[['S03','P03']][c3].sum()/6*0.5/2
# Frequency        0.000263
