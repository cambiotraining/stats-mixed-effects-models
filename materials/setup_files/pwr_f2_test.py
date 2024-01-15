# Originally developed by Emer Jones (MRC Cognition and Brain Sciences Unit, Cambridge University)
# adapted by M van Rongen (Nov 2022)

def pwr_f2_test(u=None,v=None,f2=None,sig_level=None,power=None):
    """
    Arguments:
        u: numerator degrees of freedom (#number of groups/parameters - 1)
	v: denominator degrees of freedom (#number of observations - #number of groups/parameters 
        f2: Cohen's f2 effect size value
	sig_level: significance level we're working with
	power: desired power of test 
    
    Only four out of the five parameters must be specified. The function returns the missing fifth value.
    """
    if power is None:
        power = 1 - ncfdtr(u , v , f2*(u+v+1) , f.isf(q=sig_level , dfn=u , dfd=v))
    
    elif u is None:
        def findu(u):
           return 1 - ncfdtr(u , v , f2*(u+v+1) , f.isf(q=sig_level , dfn=u , dfd=v)) - power
        u= brenth(findu , 1+1e-10, 100)
    
    elif v is None:
        def findv(v):
           return 1 - ncfdtr(u , v , f2*(u+v+1) , f.isf(q=sig_level , dfn=u , dfd=v)) - power
        v= brenth(findv , 1+1e-10, 1e9)
    
    elif f2 is None:
        def findf2(f2):
           return 1 - ncfdtr(u , v , f2*(u+v+1) , f.isf(q=sig_level , dfn=u , dfd=v)) - power
        f2= brenth(findf2 , 1e-7, 1e7)
    
    elif sig_level is None:
        def findsig(sig_level):
           return 1 - ncfdtr(u , v , f2*(u+v+1) , f.isf(q=sig_level , dfn=u , dfd=v)) - power
        sig_level= brenth(findsig , 1e-10, 1-1e10)
    
        
    #return {"u is":u,"v is":v,"f2 is":f2,"sig_level is":sig_level,"power is":power,"num_obs is":int(ceil(u))+int(ceil(v))+1}
    return print("Power analysis results: \n u is: {},\n v is: {},\n f2 is: {},\n sig_level is: {},\n power is: {},\n num_obs is: {}".format(u,v,f2,sig_level,power,int(ceil(u))+int(ceil(v))+1))
