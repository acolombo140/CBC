import numpy as np  
  
def M_dyn_r18_lambda(mns1,mns2,lamb1,lamb2):
    """
    Mass M_dyn [Msun] of the dynamical ejecta from NSNS merger. Fit by [Radice et al. 18].
    Function of lambda parameter of NS (dimensionless quadrupolar tidal deformability).
    Input: 
         mns1 = Primary NS mass [Msun]
         mns2 = Secondary NS mass [Msun]
         lamb1 = Primary NS dimensionless quadrupolar tidal deformability[-]
         lamb2 = secondary NS dimensionless quadrupolar tidal deformability[-]		
    """
    alpha = -0.657
    beta = 4.254
    gamma = -32.61
    delta = 5.202
    n = -0.773

    comp1 = compactness(lamb1)
    mnsb1 = m_to_mb(mns1,comp1)
    comp2 = compactness(lamb2)
    mnsb2 = m_to_mb(mns2,comp2)  
    
    f12 = ( alpha * (mns2/mns1)**(1./3.) * (1.-2.*comp1)/comp1 + beta * (mns2/mns1)**n + gamma * (1. - mns1/mnsb1) ) * mnsb1
    f21 = ( alpha * (mns1/mns2)**(1./3.) * (1.-2.*comp2)/comp2 + beta * (mns1/mns2)**n + gamma * (1. - mns2/mnsb2) ) * mnsb2
    
    mfit = 1.e-3 * ( f12 + f21 + delta )
    return np.maximum(mfit,0)
    
def M_dyn_d17_lambda(mns1,mns2,lamb1,lamb2):
    """
    Mass M_dyn [Msun] of the dynamical ejecta from NSNS merger. Fit by [Dietrich & Ujevic 17].
    Function of lambda parameter of NS (dimensionless quadrupolar tidal deformability).
    Input: 
         mns1 = Primary NS mass [Msun]
         mns2 = Secondary NS mass [Msun]
         lamb1 = Primary NS dimensionless quadrupolar tidal deformability[-]
         lamb2 = secondary NS dimensionless quadrupolar tidal deformability[-]		
    """
    alpha = -1.35695
    beta = 6.11252
    gamma = -49.43355
    delta = 16.1144
    n = -2.5484

    comp1 = compactness(lamb1)
    mnsb1 = m_to_mb(mns1,comp1)
    comp2 = compactness(lamb2)
    mnsb2 = m_to_mb(mns2,comp2)  
    
    f12 = ( alpha * (mns2/mns1)**(1./3.) * (1.-2.*comp1)/comp1 + beta * (mns2/mns1)**n + gamma * (1. - mns1/mnsb1) ) * mnsb1
    f21 = ( alpha * (mns1/mns2)**(1./3.) * (1.-2.*comp2)/comp2 + beta * (mns1/mns2)**n + gamma * (1. - mns2/mnsb2) ) * mnsb2
    
    mfit = 1.e-3 * ( f12 + f21 + delta )
    return np.maximum(mfit,0)    
    

def M_disk_r18_lambda(mns1,mns2,lamb1,lamb2):
    """
    Mass M_disk [Msun] of the disk ejecta from NSNS merger. Fit by [Radice et al. 18].
    Function of lambda parameter of NS (dimensionless quadrupolar tidal deformability).
    Input: 
         mns1 = Primary NS mass [Msun]
         mns2 = Secondary NS mass [Msun]
         lamb1 = Primary NS dimensionless quadrupolar tidal deformability[-]
         lamb2 = secondary NS dimensionless quadrupolar tidal deformability[-]		
    """
    alpha = 0.084
    beta = 0.127
    gamma = 567.1
    delta = 405.14
    
    lamb_tilde = 16./13. * ( (mns1 + 12.*mns2) * mns1**4 * lamb1 + (mns2 + 12.*mns1) * mns2**4 * lamb2 ) / (mns1 + mns2)**5 
    
    mfit = alpha + beta * np.tanh( (lamb_tilde - gamma) / delta )
    return np.maximum(mfit,1.e-3)

def v_dyn_r18_lambda(mns1,mns2,lamb1,lamb2):
    """
    Velocity [c] of the dynamical ejecta from NSNS merger. Fit by [Radice et al. 18].
    Function of lambda parameter of NS (dimensionless quadrupolar tidal deformability).
    Input: 
         mns1 = Primary NS mass [Msun]
         mns2 = Secondary NS mass [Msun]
         lamb1 = Primary NS dimensionless quadrupolar tidal deformability[-]
         lamb2 = secondary NS dimensionless quadrupolar tidal deformability[-]		
    """
    alpha = -0.287
    beta = 0.494
    gamma = -3.

    comp1 = compactness(lamb1)
    comp2 = compactness(lamb2)
    
    g12 = alpha * (mns1 / mns2) * (1. + gamma * comp1)
    g21 = alpha * (mns2 / mns1) * (1. + gamma * comp2)
    
    vfit = g12 + g21 + beta
    return np.maximum(vfit,0.)

def v_dyn_d17_lambda(mns1,mns2,lamb1,lamb2):
    """
    Velocity [c] of the dynamical ejecta from NSNS merger. Fit by [Dietrich & Ujevic 17].
    Function of lambda parameter of NS (dimensionless quadrupolar tidal deformability).
    Input: 
         mns1 = Primary NS mass [Msun]
         mns2 = Secondary NS mass [Msun]
         lamb1 = Primary NS dimensionless quadrupolar tidal deformability[-]
         lamb2 = secondary NS dimensionless quadrupolar tidal deformability[-]		
    """
    comp1 = compactness(lamb1)
    comp2 = compactness(lamb2) 
   
    alpha = -0.219479
    beta = 0.444836
    gamma = -2.67385    
    g12 = alpha * (mns1 / mns2) * (1. + gamma * comp1)
    g21 = alpha * (mns2 / mns1) * (1. + gamma * comp2)    
    v_rho = g12 + g21 + beta
    
    alpha = -0.315585
    beta = 0.63808
    gamma = -1.00757    
    g12 = alpha * (mns1 / mns2) * (1. + gamma * comp1)
    g21 = alpha * (mns2 / mns1) * (1. + gamma * comp2)    
    v_z = g12 + g21 + beta   
    
    vfit = np.sqrt(v_rho**2 + v_z**2)
    
    return np.maximum(vfit,0.)

def m_to_mb(mns,comp):
    return mns * ( 1. + (0.6 * comp / (1. - 0.5 * comp)) )


def compactness(lamb):
	return 0.36 - 0.0355*np.log(lamb) + 0.000705*np.log(lamb)**2




def l_td (m1,m2,l1,l2):
	q = m2/m1
	return 16./13. * ( (12.*q + 1.)*l1 + (12. + q)*q**4 * l2 ) / ( 1. + q )**5
def l_Om(lt,l0,a,mi,mj,b):
	return (lt/l0)**a * (mj/mi)**b
def x_Om(mi,mj,lt,l0,a,b):
	r = 2. * ( 1./(1. + mj/mi) + 1./l_Om(lt,l0,a,mi,mj,b) - 1. )
	if(np.isscalar(mi)):
		if(r<0.):
			return 0.
		elif(r>1.):
			return 1.
	else:
		r[r<0.] = 0.
		r[r>1.] = 1.	
	return r
    
def m_disc_Om(m1,m2,l1,l2,l0=245.,a=0.097,b=0.241):
	lt = l_td(m1,m2,l1,l2)
	x1 = x_Om(m1,m2,lt,l0,a,b)
	x2 = x_Om(m2,m1,lt,l0,a,b)
	md1 = 0.25 * (2. + x1) * (x1 - 1.)**2 * m1
	md2 = 0.25 * (2. + x2) * (x2 - 1.)**2 * m2
	#return np.maximum(1.e-3,md1+md2) #lower limit
	return np.maximum(md1+md2,0.)


def mdyn_fit_Kruger(mns1,mns2,lamb1,lamb2):
    a = -9.3335
    b = 114.17
    c = -337.56
    n = 1.5465	    
    comp1 = compactness(lamb1)
    comp2 = compactness(lamb2) 
    f1 = 1.e-3 * mns1 * (a/comp1 + b*(mns2/mns1)**n + c*comp1)  
    f2 = 1.e-3 * mns2 * (a/comp2 + b*(mns1/mns2)**n + c*comp2)  
    return np.maximum(f1+f2,0.)

def mdisk_fit_Kruger(m2,lamb2):
	a = -8.1324
	c = 1.4820
	d = 1.7784
	comp2 = compactness(lamb2) 
	return m2 * np.maximum(5.e-4,(a*comp2+c)**d)
    
    
