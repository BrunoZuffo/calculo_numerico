import numpy as np


def CalculoCondutancia():
    pi=np.pi
    Ak=2.5*10**(-7)
    mi=0.001
    Dk=np.sqrt(4*Ak/pi)
    kk=(pi*Dk**4)/(128*mi)
    Lk=float(input("Digite o comprimento do cano:"))
    Ck=kk/Lk #Lk é o comprimento do cano

    print(f"A condutancia é: {Ck}")
    return Ck
CalculoCondutancia()
