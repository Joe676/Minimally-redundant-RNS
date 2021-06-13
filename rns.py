from timeit import default_timer as timer
from random import randint

def inverse_mod(a, n): #multiplicative inverse of a (mod n)
    u = 1
    w = a
    x = 0
    z = n 
    while w != 0:
        if w < z:
            u, x = x, u
            w, z = z, w
        q = w//z
        u = u - q*x
        w = w - q*z
   
    if z==1 and x < 0:
        x+=n
    return x

class RNS: # Regular RNS factory
    def __init__(self, base=[255, 256, 257]):
        self.base = base.copy()
        self.M = 1 #dynamic range limit
        self.M_ex = [] #dynamic range excluding each base modulo
        self.inverses = [] #inverses of the above
        for m in self.base:
            self.M *= m
        for m in self.base:
            self.M_ex.append(self.M//m)
        for i in range(len(self.base)):
            self.inverses.append(inverse_mod(self.M_ex[i], self.base[i]))
    def int_to_rns(self, n):
        return RNS_Number(system=self, num = n)
    def residues_to_num(self, vals):
        return RNS_Number(system=self, vals=vals)

class RNS_Number:
    def __init__(self, system, vals=0, num=0):
        self.system = system
        if vals==0:
            self.from_int(num)
        else:
            self.values = vals.copy()
    
    def to_int(self): #CRT conversion to integer number
        out = 0
        for i, v in enumerate(self.values):
            m = self.system.base[i]
            out = (out + self.system.M_ex[i] * ((v * self.system.inverses[i])%m)) % self.system.M
        return out

    def from_int(self, num):
        self.values = [num%n for n in self.system.base]
    
    def checkCompatibility(self, other):
        #Raise errors if the numbers have incompatible bases
        if self.system == other.system:
            return
        if len(self.values) != len(other.values):
            raise ValueError("RNS numbers not compatible for the operation")
        for i in range(len(self.system.base)):
            if self.system.base[i] != other.system.base[i]:
                raise ValueError("RNS numbers not compatible for the operation")

    # Overwriting operators
    def __add__(self, other):
        self.checkCompatibility(other)
        out = self.values.copy()
        for i,v in enumerate(other.values):
            out[i] += v
            out[i] %= self.system.base[i]
        return RNS_Number(system = self.system, vals=out)
    
    def __sub__(self, other):
        self.checkCompatibility(other)
        out = self.values.copy()
        for i,v in enumerate(other.values):
            out[i] -= v
            out[i] %= self.system.base[i]
        return RNS_Number(system = self.system, vals=out)
    
    def __mul__(self, other):
        self.checkCompatibility(other)
        out = self.values.copy()
        for i,v in enumerate(other.values):
            out[i] *= v
            out[i] %= self.system.base[i]
        return RNS_Number(system = self.system, vals=out)
    
    def __div__(self, other):
        self.checkCompatibility(other)
        return RNS_Number(system = self.system, num = self.to_int()//other.to_int())

    def __str__(self):
        return str(self.values)
    
class Redundant_RNS:
    def __init__(self, base=[255, 257, 259]):
        self.base = base.copy()
        self.base.insert(0, 2) #redundant modulo 2
        self.M = 1 #dynamic range limit
        self.M_ex = [] #dynamic range excluding each base modulo
        self.M_ex_R = [] #dynamic range excluding each base modulo (for improved CRT calculations)
        self.inverses = [] #inverses of M_ex
        self.inversesK = [] #inverses of base moduli (mod last modulo)
        for i in range(1, len(self.base)):
            self.M *= self.base[i]
        for i in range(1, len(self.base)):
            self.M_ex.append(self.M//self.base[i])
        for i in range(1, len(self.base)):
            self.M_ex_R.append(inverse_mod(self.M_ex[i-1]//self.base[-1], self.base[i]))
        for i in range(1, len(self.base)):
            self.inverses.append(inverse_mod(self.M_ex[i-1], self.base[i]))
        for i in range(1, len(self.base)-1):
            self.inversesK.append(inverse_mod(self.base[i], self.base[-1]))
        
        self.LUT = [dict() for _ in self.base] #Look Up Table of R and X values
        for i in range(1, len(self.LUT)):
            m = self.base[i]
            for x in range(m):
                R = 0
                if i < len(self.inversesK):
                    R = (-((self.M_ex_R[i-1]*x)%self.base[i])*self.inversesK[i-1])%self.base[-1]
                else:
                    R = (self.M_ex_R[i-1]*x)%m
                X = (self.inverses[i-1]*x)%self.base[i]%2
                self.LUT[i][x] = (R, X)

    def int_to_rns(self, n):
        return Redundant_RNS_Number(system=self, num = n)

    def residues_to_num(self, vals):
        return Redundant_RNS_Number(system=self, vals=vals)

class Redundant_RNS_Number(RNS_Number):
    def __init__(self, system, vals=0, num=0):
        self.system = system
        if vals==0:
            self.from_int(num)
        else:
            self.values = vals.copy()
        
    def to_int(self): #improved CRT conversion
        out = 0
        RX = [self.system.LUT[i][self.values[i]] for i in range(1, len(self.system.base))] #getting the appropriate data from the lookup table
        Rs, Xs = list(zip(*RX)) #unpacking the data
        numOfOverflows = sum(Rs)//self.system.base[-1] #ro
        delta = (self.values[0] + sum(Xs) + numOfOverflows % 2) % 2 

        for i, v in enumerate(self.values[1:]):
            m = self.system.base[i+1]
            out += self.system.M_ex[i] * ((v * self.system.inverses[i])%m)
        return out-(numOfOverflows+delta)*self.system.M

    def from_int(self, num):
        self.values = [num%n for n in self.system.base]
    
    def checkCompatibility(self, other):
        #Raise errors if the numbers have incompatible bases
        if self.system == other.system:
            return
        if len(self.values) != len(other.values):
            raise ValueError("RNS numbers not compatible for the operation")
        for i in range(len(self.system.base)):
            if self.system.base[i] != other.system.base[i]:
                raise ValueError("RNS numbers not compatible for the operation")

    # Overwriting operators
    def __add__(self, other):
        self.checkCompatibility(other)
        out = self.values.copy()
        for i,v in enumerate(other.values):
            out[i] += v
            out[i] %= self.system.base[i]
        return Redundant_RNS_Number(system = self.system, vals=out)
    
    def __sub__(self, other):
        self.checkCompatibility(other)
        out = self.values.copy()
        for i,v in enumerate(other.values):
            out[i] -= v
            out[i] %= self.system.base[i]
        return Redundant_RNS_Number(system = self.system, vals=out)
    
    def __mul__(self, other):
        self.checkCompatibility(other)
        out = self.values.copy()
        for i,v in enumerate(other.values):
            out[i] *= v
            out[i] %= self.system.base[i]
        return Redundant_RNS_Number(system = self.system, vals=out)
    
    def __div__(self, other):
        self.checkCompatibility(other)
        return Redundant_RNS_Number(system = self.system, num = self.to_int()//other.to_int())

    def __str__(self):
        return str(self.values)
            


if __name__ == '__main__':
    base = [5, 7, 9, 11]
    rrns = Redundant_RNS(base = base)
    rns = RNS(base = base)
    
    a = int(input('a:'))

    A = rns.int_to_rns(a)

    Ar = rrns.int_to_rns(a)

    print(A, Ar)

    start = timer()
    for _ in range(1000000):
        A.to_int()
    end = timer()
    print(end-start)

    start = timer()
    for _ in range(1000000):
        Ar.to_int()
    end = timer()
    print(end-start)
