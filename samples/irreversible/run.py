#!/usr/bin/env python

'''
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py irr.3.out 1.25e-2 20000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py irr.2.out 1.25e-3 20000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py irr.1.out 1.25e-4 7000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py irr.0.out 1.25e-5 5000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py irr.-1.out 1.25e-6 2000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py irr.-2.out 1.25e-7 2000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py irr.-3.out 1.25e-8 1000000 &
'''

import sys
from egfrd import *
#from bd import *

def run( outfilename, T, N ):

    outfile = open( outfilename, 'w' )

    for i in range( N ):
        d, t = singlerun1( T )
        outfile.write( '%.18g\n' % d )
        outfile.flush()
        #print i
        #print d, t
        assert d == 0 or t == T


    outfile.close()



def singlerun1( T ):

    w = World(1e-3, 3)
    s = BDSimulator(w)

    #s.setMaxShellSize( 1e-6 )


    sigma = 5e-9
    r0 = sigma
    D = 2e-12
    D_tot = D
    kf = 10 * sigma * D_tot

    m = ParticleModel()

    A = m.new_species_type( 'A', 0.0, sigma/2 )
    B = m.new_species_type( 'B', D, sigma/2 )
    C = m.new_species_type( 'C', 0.0, sigma/2 )

    r1 = createBindingReactionRule( A, B, C, kf )
    m.network_rules.add_reaction_rule( r1 )

    s.setModel( m )

    particleA = s.placeParticle( A, [0,0,0] )
    particleB = s.placeParticle( B, [(float(A['radius']) + float(B['radius']))+1e-23,0,0] )

    endTime = T
    s.step()

    while 1:
        nextTime = s.getNextTime()
        if nextTime > endTime:
            s.stop( endTime )
            break
        s.step()
        if s.lastReaction:
            print 'reaction'
            return 0.0, s.t

    distance = s.distance_between_particles(particleB, particleA)

    return distance, s.t


def singlerun2( T ):

    w = World(1e-3, 3)
    s = EGFRDSimulator(w)

    #s.setUserMaxShellSize( 1e-7 )
    #s.setUserMaxShellSize( 1e-3 )

    sigma = 5e-9
    r0 = sigma
    D = 1e-12
    D_tot = D * 2

    kf = 100 * sigma * D_tot

    m = ParticleModel()

    A = m.new_species_type( 'A', D, sigma/2 )
    B = m.new_species_type( 'B', D, sigma/2 )
    C = m.new_species_type( 'C', D, sigma/2 )

    r1 = createBindingReactionRule( A, B, C, kf )
    m.network_rules.add_reaction_rule( r1 )

    s.setModel( m )

    particleA = s.placeParticle( A, [0,0,0] )
    particleB = s.placeParticle( B, [float(A['radius']) + float(B['radius'])+1e-23,0,0] )

    endTime = T

    while 1:
        s.step()
        if s.lastReaction:
            #print 'reaction'
            return 0.0, s.t

        nextTime = s.getNextTime()
        if nextTime > endTime:
            s.stop( endTime )
            break

    distance = s.distance_between_particles(particleA, particleB) 

    return distance, s.t

    

if __name__ == '__main__':
    run( sys.argv[1], float( sys.argv[2] ), int( sys.argv[3] ) )
