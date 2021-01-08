import numpy as np
from astropy import units as u
from astropy import constants as const

#constants and units
H0 = 67.9 * u.km/u.s/u.Mpc
c = const.c
G = const.G.to(u.km**2 * u.Mpc /u.Msun/u.s**2)
chirpM = 26*u.Msun
rho_crit = (3/8/np.pi) * H0**2 /G # has dimensions of density (mass/volume)

#from O2 Stochastic paper
fref = 25. /u.s
OmegaGW_limit = 4.8e-8
#From O2 rates paper
Rate = 53.2 /(u.Gpc**3 * u.year)
t_obs = 1.5 *u.year


#mass distribution parameters - here from Picky Partners paper
m_min = 6.7 #minimum allowed component mass
m_max = 41.9 #maximum allowed component mass
gamma = -1.4 #power law spectral index for individual component masses
beta_q = 4 #spectral index for selection function: f(q) \prop q^beta_q
