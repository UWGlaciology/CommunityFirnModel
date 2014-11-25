

T10m=244.0

### Initial BCO,LIZ, and DIP ###
#Vc = (6.95e-4)*T10m-0.043 #Martinerie et al., 1994, Eq. 2: critical pore volume at close off
#bcoMartRho = ((Vc+1/(917.0*(1e-3)))**-1)*1000 # Martinerie density at close off
bcoMartRho = 1/( 1/(917.0) + T10m*6.95E-7 - 4.3e-5) # Martinerie density at close off; see Buizert thesis (2011), Blunier & Schwander (2000), Goujon (2003)
#bcoAgeMart = min(age[rho>=bcoMartRho])/c['sPerYear'] # close-off age from Martinerie
#bcoDepMart = min(z[rho>=(bcoMartRho)])
#bcoAgeMartAll.append(bcoAgeMart) #age at the 815 density horizon
#bcoDepMartAll.append(bcoDepMart) #this is the 815 close off depth

# bubble close-off age and depth assuming rho_crit = 815kg/m^3
#bcoAge815 = min(age[rho>=(c['rho2'])])/c['sPerYear'] #close-off age where rho=815 kg m^-3
#bcoDep815 = min(z[rho>=(c['rho2'])]) #depth of 815 horizon
#bcoAge815All.append(bcoAge815) #age at the 815 density horizon
#bcoDep815All.append(bcoDep815) #this is the 815 close off depth

phiC = 1-bcoMartRho/917.0; #porosity at close off from mart

s=1-bcoMartRho/917.0;

phiclosed=0.37*s*(s/phiC)**(-7.6)
