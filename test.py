import unittest
from rns import Redundant_RNS_Number, Redundant_RNS, inverse_mod
from random import randint

class TestInverseMod(unittest.TestCase):

    def test_inverse_mod(self):
        self.assertEqual(inverse_mod(2, 5), 3)
        self.assertEqual(inverse_mod(4, 11), 3)
        self.assertEqual(inverse_mod(7, 13), 2)

class TestRnsMethods(unittest.TestCase):
    
    def test_init(self):
        r = Redundant_RNS([5, 7, 9, 11])
        self.assertEqual(r.base, [2, 5, 7, 9, 11])
        self.assertEqual(r.M, 5*7*9*11)
        self.assertEqual(r.M_ex, [7*9*11, 5*9*11, 5*7*11, 5*7*9])
        self.assertEqual(r.inverses, [inverse_mod(r.M_ex[i-1], r.base[i]) for i in range(1, len(r.base))])

    def test_int_to_rns(self):
        b = [5, 7, 9, 11]
        r = Redundant_RNS(b)
        b.insert(0, 2)
        for _ in range(10):
            x = randint(0, r.M)
            self.assertEqual(r.int_to_rns(x).values, [x % m for m in b])
    
    def test_residues_to_num(self):
        r = Redundant_RNS([5, 7, 9, 11])
        self.assertEqual(r.residues_to_num(vals = [1, 4, 6, 8, 10]).values, [1, 4, 6, 8, 10])
    
class TestRnsNumberMethods(unittest.TestCase):
    def test_init(self):
        b = [5, 7, 9, 11]
        r = Redundant_RNS(b)
        a = Redundant_RNS_Number(r, num=1)
        self.assertEqual(a.system, r)
        self.assertEqual(a.values, [1, 1, 1, 1, 1])
    
    def test_to_int(self):
        b = [5, 7, 9, 11]
        r = Redundant_RNS(b)
        a = Redundant_RNS_Number(r, num=1)
        for n in range(1, r.M):
            self.assertEqual(a.to_int(), n)
            a = Redundant_RNS_Number(r, num=n+1)
    
    def test_from_int(self):
        b = [5, 7, 9, 11]
        r = Redundant_RNS(b)
        a = Redundant_RNS_Number(r, num=1)
        for n in range(r.M):
            a.from_int(n)
            self.assertEqual(a.values, [n%m for m in r.base])

    def test_check_compatibility(self):
        r1 = Redundant_RNS([5, 7, 9, 11])
        r2 = Redundant_RNS([5, 7, 9, 17])
        a = r1.int_to_rns(50)
        b = r2.int_to_rns(50)
        with self.assertRaises(ValueError):
            a + b
    
    def test_operators(self):
        b = [5, 7, 9, 11]
        r = Redundant_RNS(b)
        a = r.int_to_rns(50)
        b = r.int_to_rns(25)
        #__eq__
        self.assertFalse(a==b)
        self.assertTrue(a==a)
        #__ne__
        self.assertTrue(a!=b)
        self.assertFalse(a!=a)
        #__add__
        self.assertEqual(a+b, Redundant_RNS_Number(r, num=75))
        #__sub__
        self.assertEqual(a-b, Redundant_RNS_Number(r, num=25))
        #__mul__
        self.assertEqual(a*b, Redundant_RNS_Number(r, num=25*50))
        #__div__
        self.assertEqual(a//b, Redundant_RNS_Number(r, num=2))
        self.assertEqual(a/b, Redundant_RNS_Number(r, num=2))
        #__lt__
        self.assertTrue(b<a)
        self.assertFalse(a<b)
        #__le__
        self.assertTrue(b<=a)
        self.assertTrue(a<=a)
        self.assertFalse(a<=b)
        #__gt__
        self.assertTrue(a>b)
        self.assertFalse(b>a)
        #__ge__
        self.assertTrue(a>=b)
        self.assertTrue(a>=a)
        self.assertFalse(b>=a)


if __name__ == '__main__':
    unittest.main()
