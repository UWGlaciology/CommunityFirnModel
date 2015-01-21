'''
Created on Aug 23, 2013

@author: Max
'''
import matplotlib.pyplot as plt
import os


def makeplots(plotting,Z_P,phi,gas_meas,meas_depth,meas_conc,ResultsPlace,
              diffu_full_Sev,diffu_full_fre,diffu_full_sch,diffu_full_Christo, meas_uncert=None):
        # Plotting
    
    phi=(phi-1)*1000
    meas_conc=(meas_conc-1)*1000
    meas_uncert=meas_uncert*1000
    
    if meas_uncert is None:
        pass
    
    if plotting == 'off':
        print "Plotting is turned off."
    
    elif plotting == 'on':
        
        
        fig1=plt.figure(1)
        plt.clf()
        #plt.plot(Z_P,phi[:,-8],'b')
        #plt.plot(Z_P,phi[:,-4],'r')
        #plt.plot(Z_P,phi[:,-6],'g')
        #plt.plot(Z_P,phi[:,-7],'c')
        plt.plot(Z_P,phi[:,-5],'k')

        #plt.plot(Z_P,gas_meas,'r')
        if meas_uncert is None:
            plt.plot(meas_depth,meas_conc,'k.',markersize=12)
        else:
            plt.errorbar(meas_depth,meas_conc,yerr=meas_uncert,xerr=None,fmt='.')

            
        plt.xlabel('Depth (m)')
        plt.ylabel('$\delta^{15}N_{2}$ (per mil)')
        #plt.ylabel('${CO}_{2}$ (ppm)')
        #plt.title('Concentration of $CO_{2}$ in firn at NEEM')
        #plt.legend(('Max\'s model','Measurements'),loc=3)
        plt.grid()
        
        ax = plt.gca()
        ax.ticklabel_format(useOffset=False)
        
        fname=os.path.join(ResultsPlace,'PIRE_fig3.eps')
        plt.savefig(fname,dpi=100)
        
        fig2=plt.figure(2)
        plt.clf()
        plt.plot(Z_P,diffu_full_Sev,'r')
        plt.plot(Z_P,diffu_full_fre,'b')    
        plt.plot(Z_P,diffu_full_sch,'g')    
        plt.plot(Z_P,diffu_full_Christo,'m')
        plt.xlabel('Depth (m)')
        plt.ylabel('Diffusivity ($m^{2} yr^{-1})$')
        #plt.title('Diffusivity with depth for different parameterizations')
        plt.legend(('Severinghaus et al., 2001','Freitag et al., 2002','Schwander, 1988','Buizert,2011 (tuned from data)'))  
        plt.grid()

        fname=os.path.join(ResultsPlace,'ESS524_fig2.eps')
        plt.savefig(fname,dpi=100)    
#         plt.ion()
        plt.show()
          
        
        
        #
        #fig3=plt.figure(3)
        #plt.clf()
        #plt.plot(time_yr,gas_org,'b')
        #plt.xlabel('Time (yrs)')
        #plt.ylabel('CO$_{2}$ (ppm)')
        ##plt.title('Northern Hemisphere Atmospheric CO$_{2}$ Concentration')
        #plt.grid()
        #fig3.set_size_inches(10,5)
        #fname=os.path.join(ResultsPlace,'ESS524_fig3.eps')
        #
        #plt.savefig(fname,dpi=100)    
        #plt.show()  
        #files = []
        #fig2 = plt.figure()
        #ax = fig2.add_subplot(111)
    #    for i in range(aa):  # nt frames
    #        plt.plot(Z_P,phi_toplot[:,i])
    #        plt.axis((Z_P[0],Z_P[-1],270,400))
    #        plt.xlabel('Depth (m)')
    #        plt.ylabel('$CO_{2}$ (ppm)')
    #        plt.title('Concentration of $CO_{2}$ in firn at NEEM')
    #
    #        fname = 'saved_figs/'+str('%03d' %i)+'.png'
    #        plt.savefig(fname,dpi=100)
    #        plt.clf()
    
    