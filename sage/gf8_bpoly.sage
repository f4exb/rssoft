# print element as a power of alpha    
def print_apow(e):
    if e == 0:
        return "0"
    elif e == 1:
        return "a^0"
    else:    
        return "a^%s" % e.log_repr()
        

# print element as binary     
def print_abin(e):
    if e == 0:
        return "0"
    else:    
        return "%s" % e.int_repr()
        

# print polynomial with coefficients as powers of alpha
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
        

# print polynomial with coefficients in binary representation 
def print_abin_poly(P):
    poly_str = ""
    cm_list = list(P)
    for icm,cm in enumerate(cm_list):
        c = cm[0]
        M = cm[1]
        if icm != 0:
            poly_str += " + "
        poly_str += "%s*%s" % (print_abin(c),M)
    return poly_str

        
# Divides bivariate polynomial by the greatest power of X that divides it
def star(P):
    return P/P.gcd(X^(P.degree(X)))

        
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
        

# GF(8), primitive element, GF(8)[X,Y], GF(8)[Y] and polynomial variables definition   
      
Fq.<a> = GF(2^3)
R = PolynomialRing(Fq,2,'X,Y')
X,Y = R.gens()
Sy = PolynomialRing(Fq,'Y')
Sx = PolynomialRing(Fq,'X')

# execution

def run():
    P = X*Y + X^2 + a*X^3 + 1
    Q = a*X*Y + Y^2
    print "P = %s" % print_abin_poly(P)
    print "Q = %s" % print_abin_poly(Q)
    print "P+Q = %s" % print_abin_poly(P+Q)
    print "P+a = %s" % print_abin_poly(P+a)
    print "Q+a = %s" % print_abin_poly(Q+a)
    print ""
    U = 1 + X
    print "U = %s" % print_abin_poly(U)
    print "P*U = %s" % print_abin_poly(P*U)
    print "Q/Y = %s" % (Q/Y)
    print "Q/a = %s" % print_abin_poly(Q/a)
    print ""
    print "P(a,a^2) = %s" % (P(a,a^2))
    print "P(a,1) = %s" % (P(a,1))
    print "P(1,a^2) = %s" % print_abin(P(1,a^2))
    print ""
    print "P(X,0) = %s" % Sx(P(Y=0))
    print "P(0,Y) = %s" % Sy(P(X=0))
    print "Q(X,0) = %s" % Sx(Q(Y=0))
    print "Q(0,Y) = %s" % Sy(Q(X=0))
    print ""
    P1 = (X^2)*P
    print "P1(X,Y) = %s" % print_abin_poly(P1)
    print "P1*(X,Y) = %s" % star(P1)
    print "P*(X,Y)  = %s" % star(P)
    print ""
    print "P(X,Y)^[0,0] = %s" % hasse(P,0,0)
    print "P(X,Y)^[1,0] = %s" % hasse(P,1,0)
    print "P(X,Y)^[0,1] = %s" % hasse(P,0,1)
    print "P(X,Y)^[1,1] = %s" % hasse(P,1,1)
    print "P(X,Y)^[2,0] = %s" % hasse(P,2,0)
    print "P(X,Y)^[0,2] = %s" % hasse(P,0,2)
    
