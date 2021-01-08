from stochastic_cosmo.constants import *
import numpy as np

def r_hor(H0,c,G,Mchirp,Omega,fref,R):
    '''calculates the horizon distance assuming a euclidian static universe.
    Inputs must all have astropy units.'''
    Omegaf = Omega*np.power(fref,-2./3)
    GpiM53 = np.power(G*np.pi*Mchirp,-5./3)
    Hoc = H0**2 *(c**3)
    
    rhor = (9/(8 * R)) * Hoc * GpiM53 * Omegaf
    return rhor

#this assumes equal masses and a power law for each component mass
def r_hor_marg_Mchirp(H0,c,G,Omega,fref,R0):
    '''calculates the horizon distance assuming a euclidian static universe.
    Inputs must all have astropy units.'''
    
    Omegaf = Omega*np.power(fref,-2./3)
    Gpi53 = np.power(np.pi,-5./3)
    Hoc = H0**2 * (c**3)
    
    prefac = R0/(2**(1./3))
    power = (2 * gamma + 1)
    norm = power /(np.power(m_max*u.Msun,power) - np.power(m_min*u.Msun,power))    
    RofM = norm * prefac \
    * (np.power(m_max*u.Msun,power)*np.power(m_max*u.Msun*G,5./3)-np.power(m_min*u.Msun,power)*np.power(m_min*u.Msun*G,5./3)) \
    / (power+5./3)
    #(you have to do some finegling to get the 5/3 thing to cancel units with the mass in G properly but this does the trick^)
    
    rhor = (9./(8 *RofM)) * Hoc * Gpi53 * Omegaf
    return rhor
