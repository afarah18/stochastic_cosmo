from scipy import integrate as inte
from stochastic_cosmo.constants import *


#mass distributions
def one_chirp_mass():
    '''instead of a distribution in mass, we just pretend all systems have equal mass components 
    and chirp mass 30 Msun'''
    Mchirp = 30
    massfac = np.power(Mchirp,5./3.)
    return massfac

def equal_masses():
    '''lets again assume the equal masses model that does not evolve in redshift 
    (eventually we will want to evolve it in redshift so that zf has an astrophysical rather that cosmological meaning)
    '''
    prefac = 1./(2**(1./3))
    power = (2 * gamma + 1)
    norm = power /(np.power(m_max,power) - np.power(m_min,power)) 
    massfac = norm * prefac \
    * (np.power(m_max,power)*np.power(m_max,5./3)-np.power(m_min,power)*np.power(m_min,5./3)) \
    / (power+5./3)
    return massfac

def plaw_selection():
    ''' this is a power law mass distribution with a selection function 
    for the mass ratio that is proportional to (m2/m1)** beta_q'''
    def int_over_m2(m1_arr,power_of_m2):
        return (np.power(m1_arr,power_of_m2+1) - np.power(m_min,power_of_m2+1))/ (power_of_m2 +1)
    def int_over_m1(m1_arr,m2_int_arr,power_of_m1):
        integrand = m2_int_arr * np.power(m1_arr,power_of_m1) 
        return inte.trapz(integrand)
    possible_masses = np.linspace(m_min,m_max)
    p_m1_m2_integral_A= int_over_m2(possible_masses,gamma+beta_q)
    p_m1_m2_integral_B = int_over_m1(possible_masses,p_m1_m2_integral_A,gamma-beta_q)
    p_norm = 1./p_m1_m2_integral_B
        
    def p_m1_m2(m2,m1):
        f_q = np.power(m2/m1,beta_q)
        p_m1_m2 = f_q * np.power(m1 * m2, gamma)
        return p_m1_m2
    def mass_integrand(m2,m1):
        return p_m1_m2(m2,m1) * (m1 * m2) / np.power(m1 + m2, 1./3)
    def double_integral(integrand):
        int_over_m2 = np.zeros(len(possible_masses)) 
        for i in range(len(int_over_m2)):
            # does it make sense for the upper limit and m1 argument in p to be the same?
            int_over_m2[i] = inte.quad(integrand, m_min, possible_masses[i],args=possible_masses[i])[0]
        int_over_m1 = inte.trapz(int_over_m2)
        return int_over_m1
    
    #p_norm = 1./ double_integral(p_m1_m2) #this works too! its just slower.
    integral = double_integral(mass_integrand) 
    
    return p_norm * integral

def model_A():
    def p_m1(m1):
        return np.power(m1, gamma)
    def p_m1_m2(m2,m1):
        C = 1./(m1-m_min)
        if m2>m1:
            return 0
        elif m2<=m1:
#             print C * p_m1(m1)
            return C * p_m1(m1)
    def mass_integrand(m2,m1):
        return p_m1_m2(m2,m1) * (m1 * m2) / np.power(m1 + m2, 1./3)
    def double_integral(integrand):
        int_over_m2 = np.zeros(len(possible_masses)) 
        for i in range(len(int_over_m2)-1): #dont include actual m_min because then the normilization C is infinite
            int_over_m2[i] = inte.quad(integrand, m_min, m_max,args=possible_masses[i+1])[0]
        int_over_m1 = inte.trapz(int_over_m2)
        return int_over_m1
    
    possible_masses = np.linspace(m_min,m_max)
    integral = double_integral(mass_integrand) 
    return integral    

#redshift evolution
def int_over_z(z,Omega_m,Omega_Lam,Omega_r, zp=0,alpha=0,beta=0):
    ''' putting in values for zp, alpha, and beta gives the redshift evolution in Callister et al 2020'''
    Omega_k = 1 - (Omega_m + Omega_Lam + Omega_r)
    H_z = np.sqrt(Omega_Lam + Omega_r * (1 + z)**4. + Omega_m * (1 + z)**3. + Omega_k * (1 + z)**2. )
    
    z_evo = np.power(1+z,alpha)/(1+ np.power((1+z)/(1+zp),beta+alpha))
    C = 1 + np.power(1+zp,-1*(alpha+beta))
    
    integrand = z_evo * C / (H_z*np.power(1+z,4./3.))
    return integrand

#final fractional energy density
def OmegaGW(zf,Omega_m,fref,massdist=equal_masses,z_evo_params=[0,0,0]):
    fpiG = np.power(fref*np.pi*G,2./3)/3.
    rhoc2 = ( c**2 ) * ( H0**2 ) * 3 /( 8 * np.pi * G )
    massfac = massdist() * np.power(u.Msun,5./3)
    integral = inte.quad(int_over_z,0,zf,
                         args=(Omega_m, 1-Omega_m,0,z_evo_params[0],
                               z_evo_params[1],z_evo_params[2])
                        )[0]
    integral *= 1./H0
    
    return (integral * Rate * massfac * fpiG / rhoc2).to(u.m/u.m)
