# Node in the Roth-Ruckenstein's algorithm    
class RR_Node:
    def __init__(self,parent,coeff,Q):
        self.id=t
        self.Q=Q
        self.ry_set=set()
        if parent is not None:
            self.parent=parent
            self.coeff=coeff
            self.deg=parent.deg+1
        else:
            self.parent=None
            self.coeff=None
            self.deg=-1
        

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
        

# Divides bivariate polynomial by the greatest power of X that divides it
def star(P):
    return P/P.gcd(X^(P.degree(X)))
    
    
# Finds the roots of P(0,Y) of a bivariate polynomial    
def rootsy(P):
    Py=Sy(P(X=0)) # cast to univariate polynomial in Y to find roots
    return Py.roots()

    
# Run the Roth-Ruckenstein's algorithm    
def rr_run():
    u = RR_Node(None,None,Q)
    rr_dfs(u)
    
    
# Roth-Ruckenstein's recursive node processing in deep first search strategy
def rr_dfs(u):
    global t

    #print u.id, u.deg, u.Q(Y=0) == 0, u.coeff
    
    Ry=rootsy(u.Q)
    
    for ry in Ry:
        if ry not in u.ry_set: 
            u.ry_set.add(ry)
            Qv=star((u.Q)(Y=X*Y+ry[0]))
            # Optimization: anticipate behaviour at child node
            if Qv(Y=0) == 0: 
                if u.deg < k-1:
                    return u.coeff*X^u.deg+ry[0]*X^(u.deg+1) # trace back this route from node v
                else:
                    return u.coeff*X^u.deg # trace back this route from node u
            elif u.deg == k-1 and Qv(Y=0) != 0:
                return None # cancel this route 
            # construct child node
            else: 
                t += 1
                v = RR_Node(u,ry[0],Qv)
                fpart_v = rr_dfs(v) # recursive call
                # unroll child node
                if u.deg == -1: # Root node collects results
                    if fpart_v is not None:
                        F.append(fpart_v)
                else:
                    if fpart_v == None:
                        return None
                    else:
                        return u.coeff*X^u.deg + fpart_v
                                

# Display the resulting list of polynomials 
def rr_final():
    print "Done!"
    for f in F:
        print print_apow_poly(f)
        
        
#=======================================================================================       

# GF(8), primitive element, GF(8)[X,Y], GF(8)[Y] and polynomial variables definition   
      
Fq.<a> = GF(2^3)
R = PolynomialRing(Fq,2,'X,Y')
X,Y = R.gens()
Sy = PolynomialRing(Fq,'Y')


# initialisation

t=0 # nodes but root node count
F=[] # Result set of f(X) polynomials

nb_example=2;

if nb_example == 0:
    # Li Chen's example:
    k=2 # k as in RS(n,k)
    Q=a^0+a^4*X^2+a^2*X^4+a^5*Y^2+a^4*X^2*Y^2 
elif nb_example == 1:
    # Gross and al's example
    k=5 # k as in RS(n,k)
    Q=a^2*X^5*Y+a^2*X^9+a^4*Y^2+a^5*X^4*Y+a^4*X^8+a^4*X^3*Y+X^7+a^3*X^6+a^3*X*Y+a^4*X^5+a^4*Y+a^3*X^4+a^5*X^3+a^2*X^2+a^2*X+a 
else:
    # Gross and al's example - calculated by previous interpolation
    k=5 # k as in RS(n,k)
    P=a^2*X^5*Y+a^2*X^9+a^4*Y^2+a^5*X^4*Y+a^4*X^8+a^4*X^3*Y+X^7+a^3*X^6+a^3*X*Y+a^4*X^5+a^4*Y+a^3*X^4+a^5*X^3+a^2*X^2+a^2*X+a 
    Q=P/a

    
# execution

def run():
    print "Q=%s" % print_apow_poly(Q)
    rr_run()
    rr_final()
