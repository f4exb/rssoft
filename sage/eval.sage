# print element as a power of alpha    
def print_apow(e):
    if e == 0:
        return "0"
    elif e == 1:
        return "1"
    else:    
        return "a^%s" % e.log_repr()


def run():
    for i in range(7):
        print "f1(a^%d)=%s" % (i,print_apow(f1(X=a^i)))
        
    print ""
        
    for i in range(7):
        print "f2(a^%d)=%s" % (i,print_apow(f2(X=a^i)))

    print ""
        
    print Q
    print Q.roots()

Fq.<a> = GF(2^3)
R = PolynomialRing(Fq,'X')
X = R.gen()

f1=X^4+a^3*X^3+a^4*X^2+a^4*X+a
f2=a^3*X^4+a*X^3+a^4*X^2+a^3*X+a^3
Q=a^2 + a^6*X + a^3*X^2

run()
