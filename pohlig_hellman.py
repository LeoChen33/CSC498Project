from sage.all import *
import random

def PohligHellman(P, Q):
    zList, conjList, rootList = [], [], []
    n = P.order()
    for temp_factor in n.factor():
        P0 = (ZZ(n/temp_factor[0]))*P
        conjList.append(0)
        rootList.append(temp_factor[0]^temp_factor[1])
        for i in range(temp_factor[1]):
            Qpart = Q
            for j in range(1,i+1):
                Qpart = Qpart - (zList[j-1]*(temp_factor[0]^(j-1))*P)
            Qi = (ZZ(n/(temp_factor[0]^(i+1))))*Qpart
        zList.insert(i,discrete_log(Qi,P0,operation='+'))
        conjList[-1] = conjList[-1] + zList[i]*(temp_factor[0]^i)
    return crt(conjList,rootList)

if __name__ == "__main__":
    E = EllipticCurve(GF(7919), [234,75])
    P = E.gens()[0]
    n = P.order()

    Q = random.randint(n // 2, n) * P
    PohligHellman(P, Q)