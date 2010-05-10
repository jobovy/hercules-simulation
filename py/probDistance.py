import scipy as sc
def kullbackLeibler(p,q,dx=1):
    """
    NAME:
       kullbackLeibler
    PURPOSE:
       Calculate the Kullback-Leibler divergence D(p||q)
    INPUT:
       p - probability density at points i
       q - another probability density at points i
       dx - distance between points i (can be an array)
    OUTPUT:
       D(p||q)
    HISTORY:
       2010-05-09 - Written - Bovy (NYU)
    """
    return sc.sum(p*dx*sc.log(p/q))

def symmKullbackLeibler(p,q,dx=1):
    """
    NAME:
       symmKullbackLeibler
    PURPOSE:
       Calculate the symmetrixed Kullback-Leibler divergence D(p||q)+D(q||p)
    INPUT:
       p - probability density at points i
       q - another probability density at points i
       dx - distance between points i (can be an array)
    OUTPUT:
       D(p||q) +D(q||p)
    HISTORY:
       2010-05-09 - Written - Bovy (NYU)
    """
    return kullbackLeibler(p,q,dx)+kullbackLeibler(q,p,dx)

def jensenShannon(p,q,dx):
    """
    NAME:
       jensenShannon
    PURPOSE:
       Calculate the Jensen-Shannon divergence D(p||q)
    INPUT:
       p - probability density at points i
       q - another probability density at points i
       dx - distance between points i (can be an array)
    OUTPUT:
       D_JS(p||q)
    HISTORY:
       2010-05-09 - Written - Bovy (NYU)
    """
    m= (p+q)/2.
    return 0.5*kullbackLeibler(p,m,dx)+0.5*kullbackLeibler(q,m,dx)
