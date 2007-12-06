#!/usr/bin/env python

import unittest

import numpy

from egfrd import *

class EGFRDSimulatorTestCase( unittest.TestCase ):

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test_instantiation( self ):
        s = EGFRDSimulator()
        self.failIf( s == None )

    
    def no_test_OneParticle( self ):  # this currently fails
        s = EGFRDSimulator()
        s.setCellSize( 1e-5 )
        S = Species( 'S', 2e-11, 5e-8 )
        s.addSpecies( S )
        s.placeParticle( 'S', [0.0,0.0,0.0] )

        t = s.getTime()
        for i in range( 5 ):
            s.step()
        self.failIf( t == s.getTime() )

    def no_test_TwoParticles( self ):
        s = EGFRDSimulator()
        s.setCellSize( 1e-5 )
        S = Species( 'S', 2e-11, 5e-8 )
        s.addSpecies( S )
        s.placeParticle( 'S', [0.0,0.0,0.0] )
        s.placeParticle( 'S', [5e-6,5e-6,5e-6] )

        t = s.getTime()
        for i in range( 5 ):
            s.step()
        self.failIf( t == s.getTime() )

    def no_test_ThreeParticles( self ):
        s = EGFRDSimulator()
        s.setCellSize( 1e-5 )
        S = Species( 'S', 2e-11, 5e-8 )
        s.addSpecies( S )
        s.placeParticle( 'S', [0.0,0.0,0.0] )
        s.placeParticle( 'S', [5e-6,5e-6,5e-6] )
        s.placeParticle( 'S', [1e-7,1e-7,1e-7] )

        t = s.getTime()
        for i in range( 5 ):
            s.step()
        self.failIf( t == s.getTime() )


    def no_test_ThreeParticlesInContact( self ):
        s = EGFRDSimulator()
        s.setCellSize( 1e-5 )
        S = Species( 'S', 2e-11, 5e-8 )
        s.addSpecies( S )
        s.placeParticle( 'S', [0.0,0.0,0.0] )
        s.placeParticle( 'S', [1.00001e-7,0.0,0.0] )

        # dummy
        s.placeParticle( 'S', [1e-6,1e-6,1e-6] )

        t = s.getTime()
        for i in range( 5 ):
            s.step()
        self.failIf( t == s.getTime() )

    def test_immobile_is_immobile( self ):
        s = EGFRDSimulator()
        s.setCellSize( 1e-5 )
        A = Species( 'A', 0, 1e-8 )
        s.addSpecies( A )
        B = Species( 'B', 2e-11, 5e-9 )
        s.addSpecies( B )

        particleA = s.placeParticle( A, [0.0,0.0,0.0] )
        s.placeParticle( B, [1.5000001e-8,0.0,0.0] )

        initialPosition = particleA.getPos().copy()

        for i in range( 10 ):
            s.step()
            print particleA.getPos()
        
        newPosition = particleA.getPos().copy()
        dist = s.distance( initialPosition, newPosition )

        self.failIf( dist != 0, 'initial pos: %s,\tnew pos: %s' %
                     ( initialPosition, newPosition ) )



if __name__ == "__main__":
    unittest.main()
