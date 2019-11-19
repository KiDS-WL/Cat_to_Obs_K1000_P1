#!/usr/bin/env python
# From Hendrik
import math, string, sys, os

import scipy
import scipy.integrate

def norm(k_vec): # the norm of a 3d vector
    return math.sqrt(k_vec[0]**2+k_vec[1]**2+k_vec[2]**2)

def W_k(k_vec): # the Fourier transform of the survey volume
    a=k_vec[0]*l[0]/2
    b=k_vec[1]*l[1]/2
    c=k_vec[2]**2*l[2]**2/2
    return exp(-c)*math.sin(a)/a*math.sin(b)/b

def f_k(k,R): # the Fourier transform of a spherical top-hat with radius R
    y=R*k
    return 3/y**3*(math.sin(y)-y*math.cos(y))

class Cosmology:
    """This class computes various cosmological quantities like comoving,
    angular diameter, luminosity distance, lookback time etc.. Distance
    definitions are from Hogg 1999, astro-ph/9905116.
    """
    
    def __init__(self, omega_m=0.27, omega_l=0.73, h=0.7, Gamma=0.2, n_s=1.0, sigma_8=0.81):
        self.omega_m = omega_m
        self.omega_l = omega_l
        self.omega_k = 1. - self.omega_m - self.omega_l
        self.h = h
        self.c = 2.99792458E8 # speed of light in m/s
        self.pc = 3.085678E16 # parsec in metres
        self.G = 6.673E-11 # Gravitational constant
        self.M_sun = 1.98892E30 # solar mass in kg
        self.H_0 = self.h * 100. * 1.E3 / 1.E6 / self.pc # Hubble constant in SI units
        self.dh = 3000./self.h   # Hubble distance (Hogg eq. 4) in Mpc.
        self.th = 9.78e9/self.h  # Hubble time in years
        self.th_sec = 3.09e17/self.h # Hubble time in seconds
        self.Gamma=Gamma # should be calculated by gamma=omega_m*h*exp(-omega_b*(1 + sqrt(2*h)/omega_m))
        self.n_s=n_s
        self.sigma_8=sigma_8
        self.norm_int=1/(2*math.pi)**3 * 4*math.pi * scipy.integrate.quad(lambda k: k**2*self.P_L(k)*f_k(k,8.0)**2, 0, scipy.Inf)[0]
        self.A=self.sigma_8**2/self.norm_int
        self.ro_0=2.77786E11 # critical density in M_sun/Mpc**3
        self.dlnsigma_dlnM=(math.log(self.sigma_M(10.**15))-math.log(self.sigma_M(10.**5)))/(math.log(15)-math.log(5))
        return
    

    def Ez(self, z):
        """E(z) function of Hogg's equation 14"""
        e = math.sqrt(self.omega_m*(1+z)**3 + self.omega_k*(1+z)**2 \
                      + self.omega_l)
        return e
        

    def ooEz(self, z):
        """Returns 1/E(z), E(z) being Hogg's eq. 14."""
        return 1./self.Ez(z)


    def ooEzopz(self, z):
        """Returns 1/(E(z)*(1+z)), E(z) being Hogg's eq. 14."""
        return 1./(self.Ez(z)*(1+z))

    
    def dcom_los(self, z1, z2):
        """Returns the line of sight comoving distance between objects at
        redshifts z1 and z2, z2>z1. Value is in Mpc/h"""
        if z1>=z2:
            print "z2 must be greater than z1"
            return -1
        dclos = self.dh * scipy.integrate.quad(self.ooEz, z1, z2)[0]
        return dclos

    def dcom_tra(self, z1, z2):
        """Returns the transverse comoving distance (proper motion distance)
        between objects at redshift z1 and z2."""
        dcl = self.dcom_los(z1, z2)
        if self.omega_k == 0.0:
            dct = dcl
        elif self.omega_k > 0:
            dct = self.dh / math.sqrt(self.omega_k) \
                  * math.sinh(math.sqrt(self.omega_k)*dcl/self.dh)
        else:
            dct = self.dh / math.sqrt(math.fabs(self.omega_k)) \
                  * math.sin(math.sqrt(math.fabs(self.omega_k))*dcl/self.dh)
        return dct


    def dang(self, z1, z2):
        """Returns the angular diameter distance between objects at
        redshift z1 and z2."""
        dct = self.dcom_tra(z1, z2)
        return dct/(1+z2)


    def dlum(self, z1, z2):
        """Returns the luminosity distance between objects at
        redshift z1 and z2.

        WARNING!                                          WARNING!  
                   This function is untested for z1>0!
        WARNING!                                          WARNING! 
        """
        dct = self.dcom_tra(z1, z2)
        return (1+z2)/(1+z1) * dct


    def covol(self, z):
        """Returns the comoving volume element d V_c in a solid angle
        d Omaga at redshift z."""
        da = self.dang(0, z)
        return self.dh * (1+z)**2 * da**2 / self.Ez(z)

    def tlook(self, z):
        """This function returns the lookback time in units of the
        Hubble time. The Hubble time can be accessed as the attributes
        th (in years) or th_sec (in seconds)."""
        tl = scipy.integrate.quad(self.ooEzopz, 0, z)[0]
        return tl

    def DM(self, z1, z2):
        """Returns the distance modulus between objects at
        redshift z1 and z2.
        """
        x=self.dlum(z1,z2)
        return 5*math.log(x/1.e-5)/math.log(10)

    def rho_crit(self, z1):
        """Returns the critical density at z1 in SI units.
        """
        return 3*(self.Ez(z1)*self.H_0)**2/(8*math.pi*self.G)

    def Sigma_crit(self, z1, z2):
        """Returns the critical surface mass density for lenses at z1 and sources at z2 in SI units.
        """
        return self.c**2/(4*math.pi*self.G)*self.dang(0.,z2)/(self.dang(0.,z1)*self.dang(z1,z2))/(1.E6*self.pc)*self.h

    ########## Power spectrum and mass function #############

    def T_k(self, k): # the Transfer function
        q=k/self.Gamma
        T=math.log(1+2.34*q)/(2.34*q)*(1+3.89*q+(16.1*q)**2+(5.46*q)**3+(6.71*q)**4)**(-0.25)
        return T

    def H_sqd(self, a1): # the Hubble parameter
        H=(100.*self.h)**2*(self.omega_m/(a1**3)+self.omega_l)
        return H

    def D_plus(self, a2): # the growth factor
        def func(x):
            return 1/(self.omega_m/x+self.omega_l*x**2)**1.5
        integral=scipy.integrate.quad(func,0,a2)
        integral_0=scipy.integrate.quad(func,0,1)
        D_a=math.sqrt(self.H_sqd(a2))/100.*integral[0]
        D_0=math.sqrt(self.H_sqd(1))/100.*integral_0[0]
        return D_a/D_0
    
    def D_plus2(self, a2): # the growth factor
        om = self.omega_m/(a2+self.omega_m*(1.-a2)+self.omega_l*a2*(a2*a2-1.))
        ol = self.omega_l*a2*a2*a2/(a2+self.omega_m*(1.-a2)+self.omega_l*a2*(a2*a2-1.))
        g1  = 5./2.*self.omega_m/(self.omega_m**(4./7.)-self.omega_l+(1+self.omega_m/2.0)*(1.0+self.omega_l/70.0))
        g  = 5./2.*om/(om**(4./7.)-ol+(1+om/2.0)*(1.0+ol/70.0))
        return a2*g/g1
    
    def P_L(self, k): # the linear CDM power spectrum
        P=self.T_k(k)**2*k**self.n_s
        return P

    def P_L_norm(self, k): # the normalised, linear CDM power spectrum
        P=self.A*self.T_k(k)**2*k**self.n_s
        return P

    def P_L_norm_z(self, k, z): # the normalised, linear CDM power spectrum
        P=self.A*self.T_k(k)**2*k**self.n_s*self.D_plus(1/(1+z))
        return P

    def d_ln_P_L_norm(self, k): # derivative of the normalised, linear CDM power spectrum
        P=(math.log(self.P_L_norm(k+k/1000.))-math.log(self.P_L_norm(k-k/1000.)))/(math.log(k+k/1000.)-math.log(k-k/1000.))
        return P

    def d_ln_P_L_norm_z(self, k,z): # derivative of the normalised, linear CDM power spectrum
        P=(math.log(self.P_L_norm_z(k+k/1000.,z))-math.log(self.P_L_norm_z(k-k/1000.,z)))/(math.log(k+k/1000.)-math.log(k-k/1000.))
        return P

    def Delta_sq_L_norm(self, k): # the normalised, linear, dimensionless CDM power spectrum
        P=self.A*self.T_k(k)**2*k**self.n_s*k**3/(2*math.pi**2)
        return P

    def Delta_sq_L_norm_z(self, k,z): # the normalised, linear, dimensionless CDM power spectrum
        P=self.A*self.T_k(k)**2*k**self.n_s*k**3/(2*math.pi**2)*self.D_plus(1/(1+z))
        return P

    def sigma_M(self, M):
        def func(k,R):
            return k**2*self.P_L_norm(k)*f_k(k,R)
        R=(M/self.ro_0*3/4/math.pi)**(1/3.)
        integrand=scipy.integrate.quad(func, 0, scipy.Inf, args=(R), limit=50000)[0]
        return R #1/(2*math.pi**2)*integrand

    def Jenkins(self, M):
        return 0.315*self.ro_0/M**2*self.dlnsigma_dlnM*math.exp(-math.sqrt((0.61-math.log(self.sigma_M(M)))**2)**3.8)

    def f96(self, x, n_eff): # Peacock and Dodds 1996 fitting formula
        A_c=0.482*(1.+n_eff/3.)**(-0.947)
        B_c=0.226*(1.+n_eff/3.)**(-1.778)
        alpha_c=3.310*(1.+n_eff/3.)**(-0.244)
        beta_c=0.862*(1.+n_eff/3.)**(-0.287)
        V_c=11.55*(1.+n_eff/3.)**(-0.423)
        g=5./2.*self.omega_m*(self.omega_m**(4./7.)-self.omega_l+(1+self.omega_m/2)*(1+self.omega_l/70))**(-1)
        return x*((1+B_c*beta_c*x+(A_c*x)**(alpha_c*beta_c))/(1+((A_c*x)**alpha_c*g**3/(V_c*x**0.5))**beta_c))**(1/beta_c)

    def Delta_sq_NL_PD96_norm(self, k_L): # the normalised, non-linear, dimensionless CDM power spectrum from Peacock and Dodds 1996
        n_eff=self.d_ln_P_L_norm(k_L/2.)
        return self.f96(self.Delta_sq_L_norm(k_L), n_eff)

    def Delta_sq_NL_PD96_norm_z(self, k_L,z): # the normalised, non-linear, dimensionless CDM power spectrum from Peacock and Dodds 1996
        n_eff=self.d_ln_P_L_norm_z(k_L/2.,z)
        return self.f96(self.Delta_sq_L_norm_z(k_L,z), n_eff)

    def P_NL_PD96_norm(self, k): # the normalised, non-linear CDM power spectrum from Peacock and Dodds 1996
        return self.Delta_sq_NL_PD96_norm(k)*((k/self.k_L_over_k_NL_PD96(self.Delta_sq_NL_PD96_norm(k)))**3/(2*math.pi**2))**(-1)

    def P_NL_PD96_norm_z(self, k, z): # the normalised, non-linear CDM power spectrum from Peacock and Dodds 1996
        return self.Delta_sq_NL_PD96_norm_z(k,z)*((k/self.k_L_over_k_NL_PD96(self.Delta_sq_NL_PD96_norm_z(k,z)))**3/(2*math.pi**2))**(-1)

    def k_L_over_k_NL_PD96(self, Delta):
        return (1+Delta)**(-1./3.)

