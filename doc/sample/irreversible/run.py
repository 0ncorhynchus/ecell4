#!/usr/bin/env python

from egfrd import *

def run( outfilename, N ):
    print outfilename

    outfile = open( outfilename, 'w+' )

    T = .1

    for i in range( N ):
        d, t = singlerun( T )
        outfile.write( '%g\n' % d )

        print d, t
        assert d == 0 or t == T

    outfile.close()



def singlerun( T ):

    s = EGFRDSimulator()
    s.setCellSize( 1e-3 )

    s.setMaxShellSize( 1e-6 )

    A = Species( 'A', 0.0, 5e-8 )
    s.addSpecies( A )
    B = Species( 'B', 1e-11, 5e-8 )
    s.addSpecies( B )
    C = Species( 'C', 0.0, 5e-8 )
    s.addSpecies( C )
    
    r1 = BindingReactionType( A, B, C, 1e6 / N_A )
    s.addReactionType( r1 )
    
    s.placeParticle( A, [0,0,0] )
    particle = s.placeParticle( B, [1e-7+1e-18,0,0] )

    endTime = T

    s.step()

    while 1:
        nextTime = s.scheduler.getTopEvent().getTime()
        if nextTime > endTime:
            s.stop( endTime )
            break
        s.step()
        if s.populationChanged():
            print 'reaction'
            t = s.t
            del s
            return 0.0, t

    distance = s.distance( particle.getPos(), [0,0,0] )
    t = s.t
    del s
    return distance, t
    


if __name__ == '__main__':
    run( sys.argv[1], int( sys.argv[2] ) )
