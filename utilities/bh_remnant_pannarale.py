import numpy as np
import matplotlib.pyplot as plt
from . import m_ejecta_disk_lambda as medl
#import sympy
from scipy.optimize import newton


def m_torus_func(mbh,mns,x,lamb):
	m_out = medl.M_out_f18_lambda(mbh,mns,x,lamb)
	m_ej = medl.Mej_Kawaguchi16_lambda(mbh,mns,x,0.,lamb)
	m_disk = m_out - m_ej
	return max(m_disk,0.)
	
def e_func(r,x):
	s1 = (r**2 - 2.*r + x*np.sqrt(r))/(r*np.sqrt(r**2 - 3.*r + 2.*x*np.sqrt(r)))	
	return s1
	
def nu_f(mbh,mns):
	return mbh*mns/(mbh+mns)**2
	
def Mfunc(mbh,mns):
	return mbh+mns

def Z1(x):
    return(1.+np.power(1.-x**2,1./3.)*(np.power(1.+x,1./3.)+np.power(1.-x,1./3.)))
def Z2(x):
    return(np.sqrt(3.*(x**2)+Z1(x)**2))	
def R_isco_func(m,x):
	return((3.+Z2(x)-np.sign(x)*np.sqrt((3.-Z1(x))*(3.+Z1(x)+2.*Z2(x)))))  

def nu_func(nu):
	if(nu<=0.16):
		return 0.
	elif(nu>=2./9. and nu<=0.25):
		return 1.
	elif(nu>0.16 and nu<2./9.):
		return 0.5*(1. - np.cos( np.pi*(nu-0.16)/(2./9. - 0.16) ))

def mnsb_func(mns,lamb):
	comp = medl.compactness(lamb)
	return medl.m_to_mb(mns,comp)
	
def lz_func(r,x):
	s1 = (r**2 - 2.*x*np.sqrt(r) + x**2) / ( np.sqrt(r)*np.sqrt(r**2 - 3.*r + 2.*x*np.sqrt(r)) )
	return s1

#def spin_fin(mbh,mns,xi,lamb):
#	xf = sympy.Symbol('xf',positive=True,real=True)
#	sol = sympy.nsolve(xf - ( (xi*mbh**2 + lz_func(R_isco_func(mbh,xf),xf)*mbh*((1.-nu_func(nu_f(mbh,mns)))*mns + nu_func(nu_f(mbh,mns))*mnsb_func(mns,lamb) - m_torus_func(mbh,mns,xi,lamb))) / (( Mfunc(mbh,mns)*(1. - (1.-e_func(R_isco_func(mbh,xi),xi))*nu_f(mbh,mns)) - e_func(R_isco_func(mbh,xf),xf)*m_torus_func(mbh,mns,xi,lamb) )**2) ),xf,0.)
#	return (sympy.functions.re(sol))


def spin_fin(mbh,mns,xi,lamb):
	func = lambda xf: xf - ( (xi*mbh**2 + lz_func(R_isco_func(mbh,xf),xf)*mbh*((1.-nu_func(nu_f(mbh,mns)))*mns + nu_func(nu_f(mbh,mns))*mnsb_func(mns,lamb) - m_torus_func(mbh,mns,xi,lamb))) / (( Mfunc(mbh,mns)*(1. - (1.-e_func(R_isco_func(mbh,xi),xi))*nu_f(mbh,mns)) - e_func(R_isco_func(mbh,xf),xf)*m_torus_func(mbh,mns,xi,lamb) )**2) )
	return (newton(func,xi))


		
