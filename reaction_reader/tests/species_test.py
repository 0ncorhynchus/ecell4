import unittest

import species
from decorator2 import create_species


class SpeciesTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def test_creation(self):
        sp1 = species.Species()
        su1 = species.Subunit("A")
        su1.add_modification("l", binding=1)
        su1.add_modification("r", state="u")
        sp1.add_subunit(su1)
        su2 = species.Subunit("B")
        su2.add_modification("bs", state="p", binding=1)
        sp1.add_subunit(su2)

        self.assertTrue(sp1.match(sp1))
        self.assertEqual(sp1, sp1)
        self.assertEqual(str(sp1), "A(l^1,r=u).B(bs=p^1)")
        self.assertEqual(sp1, create_species("B(bs=p^3).A(l^3,r=u)"))

    def test_properties1(self):
        sp1 = create_species("A(m3,~m1,m2,m5=s5,m4=s4)")

        self.assertEqual(str(sp1), "A(~m1,m2,m3,m4=s4,m5=s5)")
        self.assertTrue(sp1.match(sp1))
        self.assertEqual(sp1.get_binding_stride(), 0)
        self.assertEqual(sp1.num_bindings(), 0)
        self.assertEqual(sp1.num_subunits(), 1)
        self.assertEqual(sp1.num_subunits("A"), 1)

    def test_properties2(self):
        sp1 = create_species("A(l^4,r^1).A(l^4,r^2).B(l^1,r^5).B(l,r^2).C(bs^5)")

        self.assertEqual(sp1.get_binding_stride(), 5)
        self.assertEqual(sp1.num_bindings(), 4)

        sp1.sort()
        self.assertEqual(str(sp1), "A(l^1,r^2).A(l^1,r^3).B(l,r^3).B(l^2,r^4).C(bs^4)")
        self.assertEqual(sp1.num_subunits(), 5)
        self.assertEqual(sp1.num_subunits("A"), 2)
        self.assertEqual(sp1.num_subunits("C"), 1)
        self.assertEqual(sp1.get_binding_stride(), 4)
        self.assertEqual(sp1.num_bindings(), 4)

    def test_matches1(self):
        sp1 = create_species("A(l^4,r^1).A(l^4,r^2).B(l^1,r^5).B(l^2,r).C(bs^5)")

        self.assertEqual(len(create_species("A").match(sp1)), 2)
        self.assertEqual(len(create_species("B").match(sp1)), 2)
        self.assertEqual(len(create_species("C").match(sp1)), 1)
        self.assertEqual(len(create_species("D").match(sp1)), 0)

        self.assertEqual(len(create_species("B(r^_)").match(sp1)), 1)
        self.assertEqual(len(create_species("B(l^_)").match(sp1)), 2)
        self.assertEqual(len(create_species("_(r^_)").match(sp1)), 3)
        self.assertEqual(len(create_species("A(r^1).B(l^1)").match(sp1)), 2)
        self.assertEqual(len(create_species("C(bs)").match(sp1)), 0)

    def test_matches2(self):
        sp1 = create_species("X(c,a=p,b=u)")

        self.assertEqual(len(create_species("_").match(sp1)), 1)
        self.assertEqual(len(create_species("X").match(sp1)), 1)
        self.assertEqual(len(create_species("X(a,b)").match(sp1)), 1)
        self.assertEqual(len(create_species("X(~c)").match(sp1)), 0)
        self.assertEqual(len(create_species("X(c^_)").match(sp1)), 0)
        self.assertEqual(len(create_species("X(b=u)").match(sp1)), 1)

    def test_matches3(self):
        sp1 = create_species("_(ps=u)")

        self.assertEqual(len(sp1.match(create_species("A(ps=u)"))), 1)
        self.assertEqual(len(sp1.match(create_species("A(ps1=u,ps2=u)"))), 0)

        sp2 = create_species("A(ps1=u,ps2=u,ps=(ps1,ps2))")
        self.assertEqual(len(sp1.match(sp2)), 2)
        self.assertEqual(
            len(sp1.match(create_species("A(ps1=u,ps2=p,ps=(ps1,ps2))"))), 1)

        sp3 = create_species("_(_1=u,ps=(_1,))")
        self.assertEqual(len(sp3.match(sp2)), 2)
        self.assertEqual(
            set(context.get("_1") for context in sp3.match(sp2)),
            set(["ps1", "ps2"]))

        sp4 = create_species("_(_1=u,_2=u,ps=(_1,_2))")
        self.assertEqual(len(sp4.match(sp2)), 2)
        self.assertEqual(
            set((context.get("_1"), context.get("_2"))
                for context in sp4.match(sp2)),
            set([("ps1", "ps2"), ("ps2", "ps1")]))


if __name__ == '__main__':
    unittest.main()
