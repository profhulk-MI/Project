
# coding: utf-8

# In[ ]:


Skip to content
Search or jump to…

Pull requests
Issues
Marketplace
Explore
 
@profhulk-MI 
Learn Git and GitHub without any code!
Using the Hello World guide, you’ll start a branch, write comments, and open a pull request.


profhulk-MI
()
Project
1
00
Code
Issues
Pull requests
Actions
Projects
Wiki
Security
Insights
Settings
Project/Project
@profhulk-MI
profhulk-MI Create Project
Latest commit a6c88b5 on May 12, 2020
 History
 1 contributor
1025 lines (912 sloc)  38.4 KB
  
#modelling of the motors

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import statistics as st
from pandas import DataFrame
#from keras.models import Sequential
#from keras.layers import Dense
#from keras.utils.np_utils import to_categorical
#from keras import layers
#from keras.layers import LSTM
from sklearn import preprocessing
import math
import difflib
#TORQUE Te=M(II-II)
 # Rotor's speed
    #Assuming the rotor at max reaches 95% of stators speed, 0.95*3000=1900
    #Assumed no of poles p=2
    #syn_speed=3000


class Pump1(object):
    
    #df=pd.read_excel('C:/green/New folder/Central park energy and treated flow rate.xlsx',parse_cols="A:D",sheetname="Sheet2")
    #df=df[(df.T!=0).all()]
    def __init__(self,Q,pd1,pd2,D1,D2,statichead):
        self.Q=Q
        self.pd1=pd1  #Inlet pipe size
        self.pd2=pd2  #Outlet pipe size
        self.D1=D1    #inlet Dia impeller
        self.D2=D2    #Outlet Dia impeller
        self.statichead=statichead #Nominal Distance between Suction and Delivery tank
        
    def veltriangle(Q,D2,D1):
       Qm3s=Q/3600
       Qls=Qm3s*1000
       
       #Assuming that the slip would be linear and lies bwn 8-0 with 0-5%
       #Assuming the syn speed to be 3000
       slipcent=(((0.95-0.99)/8)*Qls)+0.99
       Rotorspeed=slipcent*3000
       vs=(4*Qm3s)/(math.pi*pd1*pd1)
       vd=(4*Qm3s)/(math.pi*pd2*pd2)
       u1=(math.pi*D1*Rotorspeed)/60
       u2=(math.pi*D2*Rotorspeed)/60
       #Assuing beta to be 22.5
       #assuming that the spped of the flow would be twice the current
       vw2=u2-((2*vd)/(math.tan(math.radians(22.5))))
       vr1=math.sqrt(((6*vs)**2)+(u1**2))
       vr2=math.sqrt(((6*vd)**2)+((u2-vw2)**2))
       Tneed=(1000*Qm3s*vw2*(D2/2))
       T=Tneed*(1+(Q/8))
       #print('Torque',T)
       print('\n****Hydarulic Details Pump1****\n'
             'Torque',T,
             '\nFLow',Q,
             '\nRotor Speed',Rotorspeed,
             '\nSlip Percentage',(1-slipcent)*100)
       return T,Rotorspeed,vw2,u2,vr1,vr2,vs
        
    def motor_losses(P,Q):
        Qls=(Q*1000)/3600
        slipcent=(((0.95-0.99)/8)*Qls)+0.99
        Speed=3000*slipcent
        # Iron losses and Stator losses assumed to be 2%
        Iron_Losses=P*0.02
        #print('Iron_Losses',Iron_Losses)
        # Rotor Losses assumed to be s*P
        Motor_losses=(P-Iron_Losses)*(1-slipcent)
        #print('Motor_losses',Motor_losses)
        # Windage Losses assumed to be 2%
        Wind_Stray=(Motor_losses)*0.02
        #print('Wind_Stray Losses',Wind_Stray)
        # Total Motor Losses
        Mlosses=(Iron_Losses+Motor_losses+Wind_Stray)*1.2
        #print('Total Motor Losses',Mlosses)
        # Final Power to the rotor 
        Rotor_power=P-Mlosses
        #print('Rotor_power',Rotor_power)
        BSP=Rotor_power*0.8
        print('\n*****Motor Losses Summary Pump1*****\n'
              'Inlet Power\t\t',P,
              '\nTotal Motor Losses\t\t',Mlosses,
              '\nRotor Losses\t\t',Rotor_power,
              '\nBSP',BSP)
        #Brake Shaft power assuming the bearing losses to be 2%
       
        return Mlosses,BSP
    
    def pump_losses(P,Q,T,D1,D2):
        # This accounts for the losses in th volute and impeller
        Qm3s=Q/3600
        Qls=Qm3s*1000
        
        #Assuming that the slip would be linear and lies bwn 8-0 with 0-5%
        #Assuming the syn speed to be 3000
        slipcent=(((0.95-0.99)/8)*Qls)+0.99
        Rotorspeed=slipcent*3000
        Motor_losses,BSP=Pump1.motor_losses(P,Q)
        #Torque lost is assumed to be some percent of of T needed
        T,Rotorspeed,vw2,u2,vr1,vr2,vs1=Pump1.veltriangle(Q,D2,D1)
        T_lost=(T*(Q/6))

        ############ IMPELLER LOSSES #############
                          #*****Impeller Inlet Shock losses
        maxflow=math.pi*D1*b*(6*vs1)
        h1=(0.05*((Qm3s-maxflow)**2))/(2*9.81)
        print('Impeller Inlet Shock losses',h1)
                         #****** Impeller Friction Losses
        hyd_rad=(0.05*((math.pi*D2)/4)*math.sin(math.radians(40)))/(0.05+(((math.pi*D2)/4)*math.sin(math.radians(40))))
        h2=(0.05*(D2-D1)*math.pow(np.sum([vr1,vr2]),2))/(2*math.sin(math.radians(40))*hyd_rad*4*9.81)
        print('Impeller Friction Losses',h2)
                        #*****Impeller Diffusion losses
        h3=(0.25*(math.pow((vr1-vr2),2)))/(2*9.81)
        print('Impeller Diffusion losses',h3)
                        #****** Volute friction losses
        volute_throat_velo=Qm3s/((math.pi*0.116*0.116)/4)
           #volute flow coefficient(C),  volute circumferential length(le),volute_mean_diameter
        volute_circumferential_length=(math.pi*D1)/8
        volute_mean_diameter=D1/8
        C=1+(0.02*(volute_circumferential_length/volute_mean_diameter))
        h4=(C*math.pow(volute_throat_velo,3))/(2*9.81)
        print('Volute friction losses',h4)
        # *****Disk Friction losses
        w=(2*math.pi*Rotorspeed)/60
        h5=(0.05*1000*(w**3)*math.pow((D2/2),5))/(math.pow(10,9)*Qm3s)
        print('Disk Friction losses',h5)
        #return T_lost
        h_lost=np.sum(h1+h2+h3+h4+h5)
        return T_lost,Motor_losses,BSP,h_lost
       
  
    def heat(P,Q,T,Overall_effi,Rotorspeed,m,t,na,Tmax_a,time_elasped):
        #**** COMPUTATION OF HEAT GENERATED IN MOTOR AND PUMP
        if na>=0:
            time_elasped=time_elasped
        else:time_elasped=0
        # Heat generated in the motor
        T_lost,Mlosses,BSP,h_lost=Pump1.pump_losses(P,Q,T,D1,D2)
        P_lost=T_lost*((2*math.pi*Rotorspeed)/60)
        # Assuming that the 90% and 80% of Motorlosses & T_lost resp generated heat
        Heat=(Mlosses*0.8)+(P_lost*0.6)
        H_area=(2*math.pi*0.116*0.386)+(2*math.pi*0.116*0.116)
        #**** Heat Flux
        Heat_gen=Heat/H_area
        # Heat generated in the pump
        dt=BSP*(((1-Overall_effi))/(4.2*(Q/3600)*1000))
        Total_heat_gen=Heat_gen+dt
        #print('Total_heat_gen',Total_heat_gen)
        
        #***** AMBIENT TEMPERATURE VARIATION
        # Avergae temp distribution of chennai contains max average and min temp
        Temp_dist={'Jan':[30,24,19],'Feb':[31,26,21],'Mar':[33,28,24],'Apr':[34,30,26],
                   'May':[38,33,28],'Jun':[34,30,27],'Jul':[33,29,26],'Aug':[33,29,25],
                   'Sep':[32,29,25],'Oct':[31,28,25],'Nov':[30,26,23],'Dec':[29,25,22]}

        a,b,c=Temp_dist[m]
        if 0<=t<=4:
            Avg_temp=c
        elif 10<=t<=14:
            Avg_temp=a
        elif 5<=t<10:
            gr=(a-c)/5
            Avg_temp=((t-4)*gr)+c
        elif 15<=t<=23:
            gr=(a-b)/9
            Avg_temp=(a-(t-14)*gr)
        #Average Wind Speed
        #Vel_dist={'Jan':5.64,'Feb':4.92,'Mar':8.64,'Apr':10.53,'May':10.87,'Jun':13.53,'Jul':12.06,
         #         'Aug':10.74,'Sep':8.8,'Oct':8,'Nov':7.2,'Dec':7.19}
        
        
        #****** HEAT DISSIPATION - CONVECTION AND RADIATION
        
        #** FORCED CONVECTION
        #** Evaluation of heat transfer coefficient
        # Assuming that the flow is continous and turbulent i.e Re>5x10000, Hence finding prandtl and Rey no
        w=(2*math.pi*Rotorspeed)/60
        Rey=(((w*(0.110/2))*0.6)*0.380)/0.00001578
        if Rey<=500000:
            Nu=0.664*math.pow(Rey,0.5)*math.pow(0.707,0.33)
        else:Nu=(0.037*math.pow(Rey,0.8)-871)*math.pow(0.707,0.33)
        h=(0.02623/0.380)*Nu
        #print('h',h)
        Max_temp=(Total_heat_gen/h)+Avg_temp
        #print('Max_temp initial',Max_temp)
        # Raditional losses assuming that the natural is equal to radiation 
        # Assuming that the emmissivity constant is 0.75 colour blue
        boltz=5.67e-08
        Heat_loss1=2*0.75*H_area*boltz*((math.pow(Max_temp+273,4))-(math.pow(Avg_temp+273,4)))
        Heat_lossp=2*(0.75*H_area*boltz*((math.pow(Tmax_a+273,4))-(math.pow(Avg_temp+273,4))))
        Heat_loss2=Heat_loss1+(Heat_lossp*time_elasped)
    
        #print('Heat loss due to radiation and natural convection',Heat_loss)
        # Assuimg that forced convection removes 20% of the heat generated 
        Total_heat_gen=Total_heat_gen-(Total_heat_gen*0.25)-Heat_loss2
        if Total_heat_gen>0:
            Total_heat_gen=Total_heat_gen
        else:Total_heat_gen=0
        Max_temp=(Total_heat_gen/h)+Avg_temp
         # In order to prevent temp calculations to go below 0c
        if Max_temp>=Avg_temp:
            Max_temp=Max_temp
        else: Max_temp=Avg_temp
        #print('Max_temp',Max_temp)
        #print('Total_heat_gen after natural convection and radiation',Total_heat_gen)
        #print('Avg_temp',Avg_temp)
        print('\n*****Total Heat Losses Summary Pump1*****\n'
              'Total Heat Generated',Total_heat_gen,
              '\nHeat Transfer coeficients',h,
              '\nHeat loss due to radiation and natural convection',Heat_loss1,
              '\nHeat loss due to radiation and natural convection',Heat_loss2,'\t after time',time_elasped,
              '\nMax Temp attained',Max_temp,
              '\nTotal_heat_gen after natural convection and radiation',Total_heat_gen,
              '\nAverage Temp',Avg_temp)
        return Max_temp,a,c
    
    def head(Q):
        #Fitted an exponential curve  for 2 variables
        y=(39.38*math.exp(-0.03078*Q))+(-0.2239*math.exp(0.9178*Q))
        return y
    
    
    
    def fric(Q,l,pd2):
        A=(math.pi*(0.25*pd2**2))
        v=(Q/3600)/(A)
        #print('Velocity',v)
        abs_vis=0.001
        density_water=1000
        Rey=(v*pd2*density_water)/(abs_vis)  
        #print('Reynolds No',Rey)
        # Hydraulic Diameter in meters  inches
        dh=pd2/4
        # Roughness coefficient
        k=0.00015
        if Rey>=4500:
            #Turbulent Flow
            #GRADIENT DESCENT for approximating Friction Coefficient
            x=1
            while(1):
                y=((2.51/(Rey*math.sqrt(x)))+((k/dh)/3.72))-(math.exp(-(0.5/math.sqrt(x))))
                x=x+0.01*y
                y=round(y,2)
                #print('x',x,'y',y)
                if y==0:
                    break
        else: x=64/Rey
             #Laminar Flow
             #HEAD LOSS due to Friction: DARCY WEISBACH EQUATION
        
        #Darcy Weishbach Approximation
        maj_loss= 4*x*(l/dh)*(v**2/19.62)
        #print('Major Losses',maj_loss)
        #Chezy Formula
        m=pd2/4
        C=60
        i=(((1/20)**2)*(1/m))*l
        # Loss of head due to sudden enlargement
        l1=0
        # Loss of head due to sudden contraction
        l2=0
        # Loss of head at the entrance to a pipe
        l3=0
        # Loss of head due at the exit of pipe
        l4=0
        # Loss of head due to bend in the pipe
        l5=3
        # Loss of head due to various pipe fittings
        l6=3
        min_loss=l1*(((v)**2)/(19.62))+l2*0.5*(v**2/19.62)+l3*0.5*(v**2/19.62)+l4*(v**2/19.62)+l5*(v**2/19.62)+l6*(v**2/19.62)
        head_loss_fric=maj_loss+min_loss
        print('\n****Friction losses Summary of Pump 1***\n'
              'Major Losses',maj_loss,
              '\nMinor Losses',min_loss,
              '\nTotal Friction Losses',head_loss_fric)
        return head_loss_fric
    
    
    def effi(Q,pd1,pd2,D1,D2,statichead,m,t,na,Tmax_a,time_elasped):
        #Converting flow in m3/hr to m3/s and l/s respectively
        
        #Assuing beta to be 22.5
        #assuming that the spped of the flow would be twice the current
        T,Rotorspeed,vw2,u2,vr1,vr2,vs=Pump1.veltriangle(Q,D2,D1)
        P=((2*math.pi*Rotorspeed)/60)*T
        #print('vw2',vw2,'u2',u2,'P',P)
        head_loss_fric=Pump1.fric(Q,statichead,pd2)
        weightpersec=(1000*9.81*Q)/3600
        waterpower=weightpersec*(statichead+head_loss_fric)
        WDpersec=(weightpersec*vw2*u2)/9.81
        Mano_eff=waterpower/WDpersec
        #print('Manometric Efficiency',Mano_eff*100)
        
        
        T_lost,Motor_losses,BSP,h_lost=Pump1.pump_losses(P,Q,T,D1,D2)
        #print('waterpower',waterpower,'BSP',BSP)
        Mech_eff=WDpersec/BSP
        #print('Mechanical Efficiency',Mech_eff*100)
        Overall_effi=Mech_eff*Mano_eff
        #print('Overall Efficiency',Mech_eff*Mano_eff*100)
        T1,a1,c1=Pump1.heat(P,Q,T,Overall_effi,Rotorspeed,m,t,na,Tmax_a,time_elasped)
        
        # Disk Friction coefficient
        #kd=10.6*math.pow(N,1.9)*math.pow(Q_ratio,-0.32)*math.pow(Rey,-0.2)
        #power_loss_diskfric=kd*1000*(N**3)*(D2**5)
        return Mano_eff,Mech_eff,Overall_effi,T1,P,a1,c1
    

class Pump2(object):
    
    #df=pd.read_excel('C:/green/New folder/Central park energy and treated flow rate.xlsx',parse_cols="A:D",sheetname="Sheet2")
    #df=df[(df.T!=0).all()]
    def __init__(self,Q,pd1,pd2,D1,D2,P,statichead):
        self.Q=Q
        self.pd1=pd1  #Inlet pipe size
        self.pd2=pd2  #Outlet pipe size
        self.D1=D1    #inlet Dia impeller
        self.D2=D2    #Outlet Dia impeller
        self.P=P      #Power Consumption
        self.statichead=statichead #Nominal Distance between Suction and Delivery tank
        
         
    def veltriangle(Q,D2,D1):
       Qm3s=Q/3600
       Qls=Qm3s*1000
       
       #Assuming that the slip would be linear and lies bwn 8-0 with 0-5%
       #Assuming the syn speed to be 3000
       slipcent=(((0.95-0.99)/8)*Qls)+0.99
       Rotorspeed=slipcent*3000
       vs=(4*Qm3s)/(math.pi*pd1*pd1)
       vd=(4*Qm3s)/(math.pi*pd2*pd2)
       u1=(math.pi*D1*Rotorspeed)/60
       u2=(math.pi*D2*Rotorspeed)/60
       vw2=(u2-((2*vd)/(math.tan(math.radians(22.5)))))*0.98
       vr1=math.sqrt(((6*vs)**2)+(u1**2))
       vr2=math.sqrt(((6*vd)**2)+((u2-vw2)**2))
       #Assuing beta to be 22.5
       #assuming that the spped of the flow would be twice the current
       
       Tneed=(1000*Qm3s*vw2*(D2/2))
       T=Tneed*(1+(Q/6))
       print('\n****Hydarulic Details Pump2****\n'
             'Torque',T,
             '\nFLow',Q,
             '\nRotor Speed',Rotorspeed,
             '\nSlip Percentage',(1-slipcent)*100)
       return T,Rotorspeed,vw2,u2,vr1,vr2,vs
        
    def motor_losses(P,Q):
        Qls=(Q*1000)/3600
        slipcent=(((0.95-0.99)/8)*Qls)+0.99
        Speed=3000*slipcent
        # Iron losses and Stator losses assumed to be 2%
        Iron_Losses=P*0.02
        #print('Iron_Losses',Iron_Losses)
        # Rotor Losses assumed to be s*P
        Motor_losses=(P-Iron_Losses)*(1-slipcent)
        #print('Motor_losses',Motor_losses)
        # Windage Losses assumed to be 2%
        Wind_Stray=(Motor_losses)*0.02
        #print('Wind_Stray Losses',Wind_Stray)
        # Total Motor Losses
        Mlosses=(Iron_Losses+Motor_losses+Wind_Stray)*1.2
        #print('Total Motor Losses',Mlosses)
        # Final Power to the rotor 
        Rotor_power=P-Mlosses
        #print('Rotor_power',Rotor_power)
        
        #Brake Shaft power assuming the bearing losses to be 2%
        BSP=Rotor_power*0.92
        print('\n*****Motor Losses Summary Pump2*****\n'
              'Inlet Power\t\t',P,
              '\nTotal Motor Losses\t\t',Mlosses,
              '\nRotor Losses\t\t',Rotor_power,
              '\nBSP',BSP)
        return Mlosses,BSP
    
    def pump_losses(P,Q,T,D2,D1):
        # This accounts for the losses in th volute and impeller
        Qm3s=Q/3600
        Qls=Qm3s*1000
        
        #Assuming that the slip would be linear and lies bwn 8-0 with 0-5%
        #Assuming the syn speed to be 3000
        slipcent=(((0.95-0.99)/8)*Qls)+0.99
        Rotorspeed=slipcent*3000
        Motor_losses,BSP=Pump2.motor_losses(P,Q)
        #Torque lost is assumed to be some percent of of T needed
        T,Rotorspeed,vw2,u2,vr1,vr2,vs=Pump2.veltriangle(Q,D2,D1)
        T_lost=(T*(Q/6))

        ############ IMPELLER LOSSES #############
                          #*****Impeller Inlet Shock losses
        maxflow=math.pi*D1*b*(6*vs)
        h1=(0.5*((Qm3s-maxflow)**2))/(2*9.81)
        print('Impeller Inlet Shock losses',h1)
                         #****** Impeller Friction Losses
        hyd_rad=(0.05*((math.pi*D2)/4)*math.sin(math.radians(40)))/(0.05+(((math.pi*D2)/4)*math.sin(math.radians(40))))
        h2=(0.05*(D2-D1)*math.pow(np.sum([vr1,vr2]),2))/(2*math.sin(math.radians(40))*hyd_rad*4*9.81)
        print('Impeller Friction Losses',h2)
                        #*****Impeller Diffusion losses
        h3=(0.25*(math.pow((vr1-vr2),2)))/(2*9.81)
        print('Impeller Diffusion losses',h3)
                        #****** Volute friction losses
        volute_throat_velo=Qm3s/((math.pi*0.116*0.116)/4)
           #volute flow coefficient(C),  volute circumferential length(le),volute_mean_diameter
        volute_circumferential_length=(math.pi*D1)/8
        volute_mean_diameter=D1/8
        C=1+(0.02*(volute_circumferential_length/volute_mean_diameter))
        h4=(C*math.pow(volute_throat_velo,3))/(2*9.81)
        print('Volute friction losses',h4)
        # *****Disk Friction losses
        w=(2*math.pi*Rotorspeed)/60
        h5=(0.05*1000*(w**3)*math.pow((D2/2),5))/(math.pow(10,9)*Qm3s)
        print('Disk Friction losses',h5)
        #return T_lost
        h_lost=np.sum(h1+h2+h3+h4+h5)
        return T_lost,Motor_losses,BSP,h_lost
  
    def heat(P,Q,T,Overall_effi,Rotorspeed,m,t,nb,Tmax_b,time_elasped):
        #**** COMPUTATION OF HEAT GENERATED IN MOTOR AND PUMP
        # Heat generated in the motor
        if nb>=0:
            time_elasped=time_elasped
        else:time_elasped=0
            
        T_lost,Mlosses,BSP,h_lost=Pump2.pump_losses(P,Q,T,D1,D2)
        P_lost=T_lost*((2*math.pi*Rotorspeed)/60)
        # Assuming that the 90% and 80% of Motorlosses & T_lost resp generated heat
        Heat=(Mlosses*0.8)+(P_lost*0.6)
        H_area=(2*math.pi*0.116*0.386)+(2*math.pi*0.116*0.116)
        #**** Heat Flux
        Heat_gen=Heat/H_area
        # Heat generated in the pump
        dt=BSP*(((1-Overall_effi))/(4.2*(Q/3600)*1000))
        Total_heat_gen=Heat_gen+dt
        #print('Total_heat_gen',Total_heat_gen)
        
        #***** AMBIENT TEMPERATURE VARIATION
        # Avergae temp distribution of chennai contains max average and min temp
        Temp_dist={'Jan':[30,24,19],'Feb':[31,26,21],'Mar':[33,28,24],'Apr':[34,30,26],
                   'May':[38,33,28],'Jun':[34,30,27],'Jul':[33,29,26],'Aug':[33,29,25],
                   'Sep':[32,29,25],'Oct':[31,28,25],'Nov':[30,26,23],'Dec':[29,25,22]}
        
        a,b,c=Temp_dist[m]
        if 0<=t<=4:
            Avg_temp=c
            
        elif 10<=t<=14:
            Avg_temp=a
        elif 5<=t<10:
            gr=(a-c)/5
            Avg_temp=((t-4)*gr)+c
        elif 15<=t<=23:
            gr=(a-b)/9
            Avg_temp=(a-(t-14)*gr)
        #Average Wind Speed
        #Vel_dist={'Jan':5.64,'Feb':4.92,'Mar':8.64,'Apr':10.53,'May':10.87,'Jun':13.53,'Jul':12.06,
        #          'Aug':10.74,'Sep':8.8,'Oct':8,'Nov':7.2,'Dec':7.19}
        
        #print('Avg_temp',Avg_temp)
        #****** HEAT DISSIPATION - CONVECTION AND RADIATION
        
        #** FORCED CONVECTION
        #** Evaluation of heat transfer coefficient
        # Assuming that the flow is continous and turbulent i.e Re>5x10000, Hence finding prandtl and Rey no
        w=(2*math.pi*Rotorspeed)/60
        Rey=(((w*(0.110/2))*0.6)*0.380)/0.00001578
        if Rey<=500000:
            Nu=0.664*math.pow(Rey,0.5)*math.pow(0.707,0.33)
        else:Nu=(0.037*math.pow(Rey,0.8)-871)*math.pow(0.707,0.33)
        h=(0.02623/0.380)*Nu
        #print('h',h)
        Max_temp=(Total_heat_gen/h)+Avg_temp
        #print('Max_temp initial',Max_temp)
        # Raditional losses assuming that the natural is equal to radiation 
        # Assuming that the emmissivity constant is 0.75 colour blue
        boltz=5.67e-08
        Heat_loss1=2*(0.75*H_area*boltz*((math.pow(Max_temp+273,4))-(math.pow(Avg_temp+273,4))))
        Heat_lossp=2*(0.75*H_area*boltz*((math.pow(Tmax_b+273,4))-(math.pow(Avg_temp+273,4))))
        Heat_loss2=Heat_loss1+(Heat_lossp*time_elasped)
        #print('Heat loss due to radiation and natural convection',Heat_loss)
        # Assuimg that forced convection removes 20% of the heat generated 
        Total_heat_gen=Total_heat_gen-(Total_heat_gen*0.3)-Heat_loss2
        if Total_heat_gen>0:
            Total_heat_gen=Total_heat_gen
        else:Total_heat_gen=0
        Max_temp=(Total_heat_gen/h)+Avg_temp
        
        # In order to prevent temp calculations to go below 0c
        if Max_temp>=Avg_temp:
            Max_temp=Max_temp
        else: Max_temp=Avg_temp
        #print('Max_temp',Max_temp)
        #print('Total_heat_gen after natural convection and radiation',Total_heat_gen)
        print('\n*****Total Heat Losses Summary Pump2*****\n'
              'Total Heat Generated',Total_heat_gen,
              '\nHeat Transfer coeficients',h,
              '\nHeat loss due to radiation and natural convection',Heat_loss1,
              '\nHeat loss due to radiation and natural convection',Heat_loss2,'\t after time',time_elasped,
              '\nMax Temp attained',Max_temp,
              '\nTotal_heat_gen after natural convection and radiation',Total_heat_gen,
              '\nAverage Temp',Avg_temp)
        
        return Max_temp,a,c
    
    def head(Q):
        #Fitted an exponential curve  for 2 variables
        y=(39.38*math.exp(-0.03078*Q))+(-0.2239*math.exp(0.9178*Q))
        return y
    
    
    
    def fric(Q,l,pd2):
        A=(math.pi*(0.25*pd2**2))
        v=(Q/3600)/(A)
        #print('Velocity',v)
        abs_vis=0.001
        density_water=1000
        Rey=(v*pd2*density_water)/(abs_vis)  
        #print('Reynolds No',Rey)
        # Hydraulic Diameter in meters  inches
        dh=pd2/4
        # Roughness coefficient
        k=0.00015
        if Rey>=4500:
            #Turbulent Flow
            #GRADIENT DESCENT for approximating Friction Coefficient
            x=1
            while(1):
                y=((2.51/(Rey*math.sqrt(x)))+((k/dh)/3.72))-(math.exp(-(0.5/math.sqrt(x))))
                x=x+0.01*y
                y=round(y,2)
                #print('x',x,'y',y)
                if y==0:
                    break
        else: x=64/Rey
             #Laminar Flow
             #HEAD LOSS due to Friction: DARCY WEISBACH EQUATION
        
        #Darcy Weishbach Approximation
        maj_loss= 4*x*(l/dh)*(v**2/19.62)
        #print('Major Losses',maj_loss)
        #Chezy Formula
        m=pd2/4
        C=60
        i=(((1/20)**2)*(1/m))*l
        # Loss of head due to sudden enlargement
        l1=0
        # Loss of head due to sudden contraction
        l2=0
        # Loss of head at the entrance to a pipe
        l3=0
        # Loss of head due at the exit of pipe
        l4=0
        # Loss of head due to bend in the pipe
        l5=8
        # Loss of head due to various pipe fittings
        l6=5    
        min_loss=l1*(((v)**2)/(19.62))+l2*0.5*(v**2/19.62)+l3*0.5*(v**2/19.62)+l4*(v**2/19.62)+l5*(v**2/19.62)+l6*(v**2/19.62)
        head_loss_fric=maj_loss+min_loss
        print('\n****Friction losses Summary of Pump 2***\n'
              'Major Losses',maj_loss,
              '\nMinor Losses',min_loss,
              '\nTotal Friction Losses',head_loss_fric)
        #print('Head loss Friction',head_loss_fric)
        
        return head_loss_fric
    
    
    def effi(Q,pd1,pd2,D1,D2,statichead,m,t,nb,Tmax_b,time_elasped):
        #Converting flow in m3/hr to m3/s and l/s respectively
        
        #Assuing beta to be 22.5
        #assuming that the spped of the flow would be twice the current
        T,Rotorspeed,vw2,u2,vr1,vr2,vs1=Pump2.veltriangle(Q,D2,D1)
        P=((2*math.pi*Rotorspeed)/60)*T
        #print('vw2',vw2,'u2',u2,'P',P)
        head_loss_fric=Pump2.fric(Q,statichead,pd2)
        weightpersec=(1000*9.81*Q)/3600
        waterpower=weightpersec*(statichead+head_loss_fric)
        WDpersec=(weightpersec*vw2*u2)/9.81
        Mano_eff=waterpower/WDpersec
        #print('Manometric Efficiency',Mano_eff*100)
        
        
        T_lost,Motor_losses,BSP,h_lost=Pump2.pump_losses(P,Q,T,D1,D2)
        # In order to show case the dec the eff of two by some percent
        #dec=0.2
        #BSP=BSP*(1-dec)
        #print('waterpower',waterpower,'BSP',BSP)
        Mech_eff=WDpersec/BSP
        #print('Mechanical Efficiency',Mech_eff*100)
        Overall_effi=Mech_eff*Mano_eff
        #print('Overall Efficiency',Mech_eff*Mano_eff*100)
        T1,a2,c2=Pump2.heat(P,Q,T,Overall_effi,Rotorspeed,m,t,nb,Tmax_b,time_elasped)
        
        # Disk Friction coefficient
        #kd=10.6*math.pow(N,1.9)*math.pow(Q_ratio,-0.32)*math.pow(Rey,-0.2)
        #power_loss_diskfric=kd*1000*(N**3)*(D2**5)
        
        return Mano_eff,Mech_eff,Overall_effi,T1,P,a2,c2
    
    
    '''
    def Msr(Ns,Nr,r,l,u_u0,lg):
        a=[]
        b=[]
        d=[]
        for i in range(0,361):
            if i<=180:
                i=i
            else:i=180-i
            j=-math.radians(i)
            k=math.radians(180-i)
            c1=(j+k)
            a.append(c1)
            p=abs(i-60)
            m=-math.radians(p)
            n=math.radians(180-p)
            c2=math.pi*r*(m+n)
            b.append(c2)
            q=abs(i-120)
            w=-math.radians(q)
            e=math.radians(180-q)
            c3=math.pi*r*(w+e)
            d.append(c3)
            total=((Ns*Nr*u_u0*l)/(2*lg))*((np.mean(a)+np.mean(b)+np.mean(d)))
            Msr=(total*3)
        return Msr
    
    '''
    
def grades(Q,T):
    Ii=1
    count=0
    p=1
    while(1):
        Te=Torque(Ii,Q)
        if p==1:                     
            If=Ii-0.1*(T-Te)
        else: If=Ii+0.1*(T-Te)
        Ii=If
        y=np.around((T-Te),2)
        #print('T',T,' ','Te',Te,' ','y',y,' ','count',count)
        count+=1
        if y==0:
            break
        elif y!=0 and count>=200:
            p=0
            count=0
            Ii=1
            
        #print(Ii)
    return Ii
    
def Torque(I,Q):    
    Ns=12   
    Nr=12    
    r=0.115
    l=0.216
    u_u0=0.0063
    lg=0.001
    a=np.linspace(0,360,100)
    b=np.linspace(0,3000,100)
    i_s_a=[]
    i_s_b=[]
    i_s_c=[]
    i_r_a=[]
    i_r_b=[]
    i_r_c=[]
    i_alpha_r=[]
    i_beta_r=[]
    # 2 machine with 50 hz frequency supply
    
    # Converting Q in m3/hr to l/s
    Q=(Q*1000)/3600
    speed_syn=3000
    slip_per=(((0.96-0.99)/5)*(Q))+0.99  
    speed_rot=slip_per*speed_syn
    w_s=2*math.pi*50
    w_r=2*math.pi*50*slip_per
            
    for i in a:
        #slip Frequency=P*N/120
        # 50 Hz supply
        #here the load acting on the machine is acounted and assumed that rotor speed is 99% and 95% at no and full load.
        # Here the copper and windage losses are accounted assumed to be 10% of total losses i.e. 10%
        # Also the effect of load in slip is considered as a function of added load
        s_a=math.cos(math.radians(w_s*i))
        r_a=math.cos(math.radians(w_r)*i)
        s_b=math.cos(math.radians(w_s*i+120))
        r_b=math.cos(math.radians(((w_r)*i)+120))
        s_c=math.cos(math.radians(w_s*i+240))
        r_c=math.cos(math.radians(((w_r)*i)+240))
        i_s_a.append(abs(s_a))
        i_s_b.append(abs(s_b))
        i_s_c.append(abs(s_c))
        i_r_a.append(abs(r_a))
        i_r_b.append(abs(r_b))
        i_r_c.append(abs(r_c))
        
    Current_stator=np.array([[np.mean(i_s_a)],[np.mean(i_s_b)],[np.mean(i_s_c)]])
    Current_rotor=np.array([[np.mean(i_r_a)],[np.mean(i_r_b)],[np.mean(i_r_c)]])
    Current_stator_2phase=(2/math.sqrt(3))*I*Current_stator
    Current_rotor_2phase=(2/math.sqrt(3))*I*Current_rotor
    Current_sr_2phase=np.concatenate((Current_stator_2phase,Current_rotor_2phase))
    
    
    #MUTUAL FLUX LINKAGES CALCULATIONS
    q=np.linspace(0,360,30)
    for i in q:  
        tran_matrix=np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,0,0,0,0],[0,0,0,math.cos(math.radians(i)),-math.sin(math.radians(i)),0],[0,0,0,math.sin(math.radians(i)),math.cos(math.radians(i)),0],[0,0,0,0,0,0]])
        Current_sr_2phase_pseudo=np.matmul(tran_matrix,Current_sr_2phase)
        #print('i',i,Current_sr_2phase_pseudo[3],' ',Current_sr_2phase_pseudo[4])
        i_alpha_r.append(abs(Current_sr_2phase_pseudo[3]))
        i_beta_r.append(abs(Current_sr_2phase_pseudo[4]))
    #Mutal_ind=Msr(Ns,Nr,r,l,u_u0,lg)
    
    #Considering 4 pole Machine
    P1=(Current_sr_2phase[1]*np.mean(i_alpha_r))
    P2=(Current_sr_2phase[0]*np.mean(i_beta_r))
    Mutal_ind=2*((Ns*Nr*u_u0*l)/(2*lg))*r*0.003402473*1.5
    if I>=10:
        Torque=(Mutal_ind*(P1+P2))
    else: Torque=(Mutal_ind*(P1+P2)) +0.5
    return Torque
a=13.73
Q=4.29555555
Torque(a/3,Q)

grades(4.295555556, 2.690275227)

   
    
def ratio(Q,pd1,pd2,D1,D2,b,statichead,m,t,na,nb,Tmax_a,Tmax_b,time_elasped1,time_elasped2):
    #Converting the efficieny in the range of 5
    #Converting the temp raise in the range of 5
    #I=grades(Q,T)
    mo1,m1,o1,Temp1,P1,a1,c1=Pump1.effi(Q,pd1,pd2,D1,D2,statichead,m,t,na,Tmax_a,time_elasped1)
    mo2,m2,o2,Temp2,P2,a2,c2=Pump2.effi(Q,pd1,pd2,D1,D2,statichead,m,t,nb,Tmax_b,time_elasped2)
    Overeff1.append(o1)
    Overeff2.append(o2)
    TMax1.append(Temp1)
    TMax2.append(Temp2)
    
    #mo2,m2,o2=Pump2.effi(Q,pd1,pd2,D1,D2,statichead)
    print('\n********PUMP1******** \t\t\t********PUMP2********\nManometric Eff=',round(mo1,2),
          '\t\t\tManometric Eff=',round(mo2,2),'\nMechanical Eff=',round(m1,2),'\t\t\tMechanical Eff=',round(m2,2),
          '\nOverall Eff=',round(o1,2),'\t\t\tOverall Eff=',round(o2,2))
    Over_eff1p=round(5*np.mean([Overeff1]),2)
    Over_eff2p=round(5*np.mean([Overeff2]),2)
    temp_1p=round(5-((((np.mean([TMax1])-np.min([TMax1]))/np.max([TMax1]))*5)),2)
    temp_2p=round(5-((((np.mean([TMax2])-np.min([TMax1]))/np.max([TMax1]))*5)) ,2)
    Points1=round(st.mean([Over_eff1p,temp_1p]),2)
    Points2=round(st.mean([Over_eff2p,temp_2p]),2)  
    MT1=round((24*Points1)/np.sum([Points1+Points2]),0)
    MT2=round((24*Points2)/np.sum([Points1+Points2]),0)
    return MT1,MT2,o1,o2,Temp1,Temp2,P1,P2


def decide(pd1,pd2,b,D1,D2,statichead,m):
        average_inlet=[4.13,3.45,3.76,4.37,3.96,2.32,3.00,3.87,4.43,4.50,4.01,3.72,3.50,4.42,3.97,3.96,
                   3.71,3.52,3.63,3.89,4.30,3.19,3.45,4.15]
       
        Power1=[]
        Power2=[] 
        PowerDiff=[]
        Eff=[]
        Effr=[]
        action_plan=[]
        rein_action_plan=[]
        time_elasped1=0
        time_elasped2=0
        slope=np.diff(average_inlet)
        slope=np.insert(slope,0,0,0)
        local_max=np.cumsum(slope)
        slope[0]=1
        T_a_award=0
        na=-1
        nb=-1
        Tmax_a=0
        Tmax_b=0
        T_b_award=0
        local_min=np.cumsum(slope)
        local_max=pd.Series(local_max)
        local_min=pd.Series(local_min)
        count_max=local_max[local_max.T>=0].index
        count_min=local_min[local_min.T<0].index
        slope=np.delete(slope,0)   
        count1=0
        count2=0
        v=0


        for i in range(0,len(average_inlet)):  
            Q=np.random.normal(average_inlet[i],0)
            t=i
            time_elasped1=t-T_a_award-1
            time_elasped2=t-T_b_award-1
            R1,R2,o1,o2,T1,T2,P1,P2,=ratio(Q,pd1,pd2,D1,D2,b,statichead,m,t,na,nb,Tmax_a,Tmax_b,time_elasped1,time_elasped2)
            Tmax_a=T1
            Tmax_b=T2
            Power1.append(round(P1,2))
            Power2.append(round(P2,2))
            PowerDiff.append(abs(P1-P2))
            choices=[1,2]
            ap=0
            bp=0
            turn=[P1,P2]
            p=0
           
            # To simulate human works
            countH=t+1
            if countH % 4!=0:
                a=0
            else:a=1
            HWorks.append(turn[a])

            
            #identification of local maxima and minima
            #Prior  Average

            ap+=math.exp(Q-np.average(average_inlet))
            bp+=math.exp(-(Q-np.average(average_inlet))) 

            print('\nAfter prior ap',ap,'','After prior bp',bp)    
            
            #Posterior
            #Rule 1 for Power Measurement, favours the higher efficiency punp one and accounts for slope info
            #Checking for any local max or min
            
         
            local_max_ava= t in count_max
            local_min_ava= t in count_min
            
            if count_max.size!=0:
                k=np.abs(t-count_max).argmin()
                k1=count_max[k]
                if t==k1:
                    count_max=np.delete(count_max,k)
                else: count_max=count_max
                if k1==0 or (t-k1)==0:
                    u=0
                else:u=k1-t
            else:u=0
            
            if count_min.size!=0:
                    q=np.abs(t-count_min).argmin()
                    q1=count_min[q]
                    if t==q1:
                        count_min=np.delete(count_min,q)
                    else:count_min=count_min
                    if q1==0 or (t-q1)==0:
                        w=0
                    else:w=q1-t
            else: w=0
            slopec=slope[t:t+3]
            #slopec=slopec.reshape(-1,1)
            #slopec=preprocessing.normalize(slopec)
            slopec=np.ndarray.tolist(slopec)
            slopec_size=np.size(slopec)
            ao1=((Power1[t])/np.max(Power1))
            bo1=((Power2[t])/np.max(Power2))
            #bo1=(1-ao1)
            if slopec_size==3:
                slopequotient=((slopec[0]+0.75*slopec[1]+0.5*slopec[2])/(2.25))
            elif slopec_size==2:
                slopequotient=((slopec[0]+0.75*slopec[1])/(1.75))
            else:slopequotient=0
                            
            if local_max_ava==True:
                ap+=np.random.normal(4,0.1)
            else: ap+=math.exp(ao1/3)+((1-math.exp(-u))/(1+math.exp(-u)))+math.exp(-slopequotient)
            #print(ap)
            if local_min_ava==True:
                bp+=np.random.normal(4,0.1)
                #Tip: Here the SLope Quotient sign in postive in the ap ans nagative in the bp, the reason being the changing slope
                # the postive prefers and vice versa
            else:bp+=math.exp(bo1/3)+((1-math.exp(-w))/(1+math.exp(-w)))+math.exp(slopequotient)
            #print(bp)
            print('\nAfter power ap',ap,'',' bp',bp)  
                                    
            # Rule2 Heat 
            #Have NOrmalized the range to 0-40 and taking Ln as the rating function

            ap+=math.log((abs(T1-80)/80)*30)+(math.exp(1-((R1-count1)/R1))*3) 
            bp+=math.log((abs(T2-80)/80)*30)+(math.exp(1-((R2-count1)/R2))*3)
            print('\nAfter Heat ap',ap,'',' bp',bp)    

            if ap>bp:
                p=1
                T_a_award=t
                count1+=1
                na=na+round((1-math.exp(-2*t))/(1+math.exp(-2*t)),0)+round(math.exp(-t),0)
                MWorks.append(P1)
            elif ap<bp:
                p=2
                T_b_award=t
                count2+=1 
                nb=nb+round((1-math.exp(-2*t))/(1+math.exp(-2*t)),0)++round(math.exp(-t),0)
                MWorks.append(P2)
            print('\n*****Working counts*****\n'
                  'Alloted Work Pump 1:',R1,'\t\tExecuted',count1,
                  '\nAlloted Work Pump 2:',R2,'\t\tExecuted',count2)
            #Eff.append(((((np.sum([HWorks]))-(np.sum([MWorks]))))/np.sum([HWorks]))*100)
        
            action_plan.append(p)
            
            #choices
            f=np.random.choice(choices)
            
            if Q>np.mean(average_inlet):
                Q1.append(Q)
            else:Q0.append(Q)
            
            if v==0 and f==1:
                choice_a1.append(f)
            elif v!=0:
                choice_a1.append(v)
            if v==0 and f==2:
                choice_a2.append(f)
            elif v!=0:
                choice_a2.append(v)
             
         
            #Decision 
            choices.remove(f)
            if np.size(choice_a1)!=0 and np.size(choice_a2)!=0:
                if Q>np.mean(average_inlet):
                        a=(np.size(pr_a1)/np.size(choice_a1))*(np.size(a1_Q1)/np.size(Q1))
                        b=(np.size(pr_a2)/np.size(choice_a2))*(np.size(a2_Q1)/np.size(Q1))
                        if f==1 and a>b:
                            v=f
                        else:v=np.random.choice(choices)
                        #b=((1-(np.size(pr_a1)))/np.size(choice_a1))*(np.size(a1_Q1)/np.size(Q1))
                elif Q<np.mean(average_inlet):
                        a=(np.size(pr_a2)/np.size(choice_a2))*(np.size(a2_Q0)/np.size(Q0))
                        b=(np.size(pr_a1)/np.size(choice_a1))*(np.size(a1_Q0)/np.size(Q0))
                        #b=((1-(np.size(pr_a1)))/np.size(choice_a1))*(np.size(a1_Q0)/np.size(Q0))   
                        #b=1-a
                        if f==2 and a>b:
                            v=f
                        else:v=np.random.choice(choices)
            elif np.size(choice_a2)==0 or np.size(choice_a1)==0:
                v=v=np.random.choice(choices)
            
                    
              #Reward Function    
            if v==1 and ap>bp:
                pr_a1.append(v)
            if v==2 and ap<bp:
                pr_a2.append(v)

            
            if v==1 and Q>np.mean(average_inlet):
                a1_Q1.append(v)
            elif v==1 and Q<np.mean(average_inlet):
                a1_Q0.append(v)
            if v==2 and Q<np.mean(average_inlet):
                a2_Q0.append(v)
            elif v==2 and Q>np.mean(average_inlet):
                a2_Q1.append(v)
            
            if v==1:
                MrWorks.append(P1)
            else:MrWorks.append(P1)
            print('a',a,'b',b,
                  '\nQ',Q,
                  '\nv',v,
                  '\nf',f)
            rein_action_plan1.append(v)
            rein_action_plan.append(v)
        Effi=(((HWorks[t]-MWorks[t])/HWorks[t])*100)
        Eff.append(Effi)
        Effe=(((HWorks[t]-MrWorks[t])/HWorks[t])*100)
        Effr.append(Effe)
        print('\nEfficiency',Effi,
              '\nEffficiency reinforced',Effe)
        return action_plan,HWorks,MWorks,Eff,Effr
    
a=input('1. Default  2. Customized')
if a=='1':
        statichead=30
        pd1=0.05
        b=0.05
        m='Dec'
        D1=0.102
        D2=0.153
        pd2=0.04
#else:statichead,pd1,b,D1,D2,pd2,m,t=input('Enter Values for 1.Static Head 2.Pipe Size Inlet 3. Impeller Width 4. Inlet Dia of Impeller 5. Outlet Dia of Impeller 6. Pipe size of Outlet 7. Current Month 8.Current time(hr)')
pr_a1=[]
pr_a2=[]
nr_a1=[]
nr_a2=[]
a1_Q0=[]
a2_Q1=[]
Q1=[]
Q0=[]
choice_a1=[]
choice_a2=[]
a1_Q1=[]
a2_Q0=[]
Overeff1=[]
Overeff2=[]
TMax1=[]
TMax2=[]
HWorks=[]
MWorks=[]
MrWorks=[]
rein_action_plan1=[]

#Q_distperday=[4.13,3.45,3.76,4.37,3.96,2.32,3.00,3.87,4.43,4.50,4.01,3.72,3.50,4.42,3.97,3.96,
#                   3.71,3.52,3.63,3.89,4.30,3.19,3.45,4.15]
#np.random.shuffle(Q_distperday)
#for i in Q_distperday:
#    Q=np.random.normal(loc=Q_distperday[0],scale=1)
p,H,M,E,Er=decide(pd1,pd2,b,D1,D2,statichead,m)
q=difflib.SequenceMatcher(None,p,rein_action_plan1)
print('Concurrence Ratio',q.ratio())
print('Efficency Ratio',(np.sum(Er)/np.sum(E)))
#l=Pump1(Q,pd1,pd2,D1,D2,statichead)
© 2021 GitHub, Inc.
Terms
Privacy
Security
Status
Help
Contact GitHub
Pricing
API
Training
Blog
About

