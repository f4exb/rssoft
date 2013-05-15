import math


# binomial coefficient (combination of k in n)
def binomial(n,kk):
    return math.factorial(n)/(math.factorial(kk)*math.factorial(n-kk))


# weighted degree of bivariate polynomial
def wdegree(P,a,b):
    wdmax=0
    m_list = P.monomials()
    for mono in m_list:
        mono_wdeg=mono.degree(X)*a+mono.degree(Y)*b
        if mono_wdeg > wdmax:
            wdmax = mono_wdeg
    return wdmax

    
# find leading monomial with respect to weighted inverse lexical order
def wlm(P,a,b):
    Ml=1*X^0*Y^0
    wdmax=0
    m_list = P.monomials()
    for mono in m_list:
        mono_wdeg=mono.degree(X)*a+mono.degree(Y)*b
        if mono_wdeg > wdmax: 
            wdmax = mono_wdeg
            Ml = mono
        elif mono_wdeg == wdmax:
            if mono.degree(X) < Ml.degree(X):
                Ml = mono
    return Ml
    

# (u,v)-ith Hasse derivative of bivariate polynomial 
def hasse(P,u,v):
    H=0*X^0*Y^0
    cm_list = list(P)
    
    for cm in cm_list:
        aij=cm[0]
        i=cm[1].degrees()[0]
        j=cm[1].degrees()[1]
        if i>=u and j>=v:
            H += binomial(i,u)*binomial(j,v)*aij*X^(i-u)*Y^(j-v)
    return H
    
    
# transform 0 order monomial to element. Used when evaluating a polynomial    
def poly2element(P):
    return P.constant_coefficient()
    
    
# print element as a power of alpha    
def print_apow(e):
    if e == 0:
        return "0"
    elif e == 1:
        return "a^0"
    else:    
        return "a^%s" % e.log_repr()
        

# print polynomial with coefficiens as powers of alpha
def print_apow_poly(P):
    poly_str = ""
    cm_list = list(P)
    for icm,cm in enumerate(cm_list):
        c = cm[0]
        M = cm[1]
        if icm != 0:
            poly_str += " + "
        poly_str += "%s*%s" % (print_apow(c),M)
    return poly_str


# set G, lodG and calc
def init_G():
    global G, calc, lodG
    inclod=1
    lod=0
    
    for i in range(dY+1):
        if i == 0:
            G=[1*X^0*Y^0]
        else:
            G.append(Y^i)
        calc.append(True)
        lodG.append(lod)
        inclod += k-1
        lod += inclod
        

# outer iteration    
def process_point(ip):
    for u in range(m[ip]):
        for v in range(m[ip]-u):
            process_hasse(x[ip],y[ip],u,v)
    
    
# inner iteration    
def process_hasse(xx,yy,u,v):
    ig_lodmin=0
    lodmin=0
    ig_wdmin=0
    wdmin=0
    first_hnn=True
    h_xy_list = []
    G_next = []
    lodG_next = []
    ind=""
    zero_hasse = True # all Hasse derivatives to be considered are zero
    global it, G, lodG, calc

    print "it=%d x=%s y=%s u=%d v=%d" % (it, print_apow(xx), print_apow(yy), u, v)
    
    # Hasse derivatives calculation
    for ig,g in enumerate(G):
        if calc[ig]:
            h=hasse(g,u,v)
            h_xy = poly2element(h(X=xx,Y=yy))
            h_xy_list.append(h_xy)
            wd=wdegree(g,1,k-1)
            
            if h_xy == 0:
                ind="-"
            else:
                zero_hasse = False
                ind="+"
                if first_hnn:
                    lodmin=lodG[ig]
                    wdmin=wd
                    ig_lodmin=ig
                    ig_wdmin=ig
                    first_hnn = False
                if lodG[ig] < lodmin:
                    lodmin=lodG[ig]
                    ig_lodmin=ig
                if wd<wdmin:
                    wdmin=wd
                    ig_wdmin=ig
        else:
            ind="x"
            h_xy_list.append(0)
                
        print "%s g_%d,%d=%s" % (ind, it, ig, print_apow_poly(g))
        if calc[ig]:
            print "  D_%d,%d=%s" % (it, ig, print_apow(h_xy))
            #print "  wd=%d" % wd
            #print "  Lm=%s" % wlm(g,1,k-1)
            print "  lod=%d" % lodG[ig]
        else:            
            print "  lod=%d" % lodG[ig]
    
    if zero_hasse:
        print "All Hasse derivatives are 0 so G_%d=G_%d" % (it+1,it)
    else:
        print "Minimal LOD polynomial g_%d,%d (ig_wdmin=%d)" % (it, ig_lodmin, ig_wdmin)   

    # compute next values in G
    for ig,g in enumerate(G):
        if calc[ig]:
            if h_xy_list[ig] == 0:
                G_next.append(g)
                lodG_next.append(lodG[ig])
            else:
                if ig == ig_lodmin:
                    G_next.append(h_xy_list[ig]*(X-xx)*g)
                    mX = (wlm(g,1,k-1)).degree(X) # find leading monomial's X power
                    mY = (wlm(g,1,k-1)).degree(Y) # find leading monomial's Y power
                    lodG_next.append(lodG[ig_lodmin]+(int(mX)/int(k-1))+1+mY) # new leading order by sliding one position of X powers to the right
                else:
                    G_next.append(h_xy_list[ig]*G[ig_lodmin]-h_xy_list[ig_lodmin]*g)
                    lodG_next.append(max(lodG[ig],lodG[ig_lodmin])) # new leading order is the max of the two
            if lodG_next[ig] > Cm: # Complexity reduction, skip polynomial processing if its lod is too big (bigger than multiplicity cost)
                calc[ig] = False
        else:
            lodG_next.append(lodG[ig])
            G_next.append(g)
                
    G = G_next
    lodG = lodG_next
    it += 1
    print ""
    

# outcoming G vector
def final_G():

    ig_lodmin=0
    lodmin=lodG[0]
    ig_wdmin=0
    wdmin=wdegree(G[0],1,k-1)
    
    print "it=%d final result" % it
    
    for ig,g in enumerate(G):    
        if lodG[ig] < lodmin:
            lodmin=lodG[ig]
            ig_lodmin=ig
        wd=wdegree(G[ig],1,k-1)             
        if wd < wdmin:
            wdmin=wd
            ig_wdmin=ig 
        print "o g_%d,%d=%s" % (it, ig, print_apow_poly(g))
        #print "  wd=%d" % wd
        #print "  Lm=%s" % wlm(g,1,k-1)
        print "  lod=%d" % lodG[ig]
        
    print "Minimal LOD polynomial g_%d,%d (ig_wdmin=%d)" % (it, ig_lodmin, ig_wdmin)   
    print "Q=%s" % print_apow_poly(G[ig_lodmin])     
            

# Run the GS iterative algorithm            
def run_iterations():
    for ip in range(len(x)):
        process_point(ip)

        
#=======================================================================================

# GF(8), primitive element, GF(8)[X,Y] and polynomial variables definition         

Fq.<a> = GF(2^3)
R = PolynomialRing(Fq,2,'X,Y')
X,Y = R.gens()


# example initialisation

lodG=[] # leading orders of polynomials in G
calc=[] # skip calculation if false 
        # (complexity reduction - Li Chen's improvement on original GS algorithm)
it=0 # Construction iteration

nb_example=1

if nb_example == 0:
    # Li Chen's example
    dY=5 # maximum power of Y in G0
    x=[  1,  a,a^3,a^2,a^6,a^4,a^5] # Interpolation point x coordinate
    y=[a^5,a^3,a^4,  0,a^6,a^2,a^2] # Interpolation point y coordinate
    m=[  2,  2,  2,  2,  2,  2,  2] # Multiplicity vector corresponding to
                                    # interpolation points (hard decision list decoding hence constant)
    Cm=21 # Multiplicity cost
    k=2 # k as in RS(n,k)

else:
    # Gross and al's example
    dY=2 # maximum power of Y in G0
    x=[  1,  a,a^2,a^3,a^3,a^4,a^5,a^6] # Interpolation point x coordinate
    y=[  0,a^3,a^3,  0,  1,  1,a^2,  0] # Interpolation point y coordinate
    m=[  2,  1,  2,  1,  1,  2,  2,  1] # Multiplicity vector corresponding to 
                                        # interpolation points and associated multiplicity
    Cm=16 # Multiplicity matrix cost
    k=5 # k as in RS(n,k)

    
# execution

def run():
    init_G()
    run_iterations()
    final_G()
