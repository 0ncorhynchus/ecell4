#!/usr/bin/env python

'''
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py rev.3.out 1.25e-2 5000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py rev.2.out 1.25e-3 4000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py rev.1.out 1.25e-4 2000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py rev.0.out 1.25e-5 2000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py rev.-1.out 1.25e-6 2000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py rev.-2.out 1.25e-7 2000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py rev.-3.out 1.25e-8 1000000 &
'''

import sys
from egfrd import *
from bd import *

def run(outfilename, T, N):
    print outfilename

    outfile = open(outfilename, 'w')

    for i in range(N):
        d, t = singlerun(T)
        outfile.write('%.18g\n' % d)
        outfile.flush()
        #print d, t
        assert d == 0 or t == T

    outfile.close()



def singlerun(T):

    w = World(1e-3, 3)
    s = EGFRDSimulator(w)
    #s.set_user_max_shell_size(1e-6)
    #s = BDSimulator(w)

    sigma = 5e-9
    r0 = sigma
    D = 1e-12
    D_tot = D * 2

    tau = sigma * sigma / D_tot

    kf = 100 * sigma * D_tot
    koff = 0.1 / tau

    m = ParticleModel()

    A = m.new_species_type('A', D, sigma/2)
    B = m.new_species_type('B', D, sigma/2)
    C = m.new_species_type('C', D, sigma/2)

    r1 = create_binding_reaction_rule(A, B, C, kf)
    m.network_rules.add_reaction_rule(r1)

    r2 = create_unbinding_reaction_rule(C, A, B, koff)
    m.network_rules.add_reaction_rule(r2)

    s.set_model(m)

    s.place_particle(A, [0,0,0])
    s.place_particle(B, [(float(A['radius']) + float(B['radius']))+1e-23,0,0])

    end_time = T
    s.step()

    while 1:
        next_time = s.get_next_time()
        if next_time > end_time:
            s.stop(end_time)
            break
        s.step()

    
    if len(s.particle_pool[C.id]) != 0:
        return 0, s.t

    distance = s.distance_between_particles(A.id, B.id)

    return distance, s.t
    
if __name__ == '__main__':
    run(sys.argv[1], float(sys.argv[2]), int(sys.argv[3]))
