from ecell4.core import *
from ecell4.ode import *

def run():
    volume = 0.5
    sp1, sp2, sp3 = Species("A"), Species("B"), Species("C")
    k = 1.0
    rr1 = create_unbinding_reaction_rule(sp1, sp2, sp3, k)

    m = NetworkModel()
    m.add_species(sp1)
    m.add_species(sp2)
    m.add_species(sp3)
    m.add_reaction_rule(rr1)

    w = ODEWorld(volume)
    w.add_species(sp1)
    w.add_species(sp2)
    w.add_species(sp3)
    w.add_molecules(sp1, 60)

    target = ODESimulator(m, w)

    next_time = 0.0
    dt = 0.01

    print "t = %g\t A = %g\t B = %g\t C = %g" % (
        target.t(), w.num_molecules(sp1), w.num_molecules(sp2),
        w.num_molecules(sp3))
    for i in range(200):
        next_time += dt
        target.step(next_time)
        print "t = %g\t A = %g\t B = %g\t C = %g" % (
            target.t(), w.num_molecules(sp1), w.num_molecules(sp2),
            w.num_molecules(sp3))
        

if __name__ == "__main__":
    run()

