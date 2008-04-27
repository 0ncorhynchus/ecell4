#!/usr/env python

import weakref

import math

import numpy


from utils import *
from surface import *

from gfrdbase import *
import _gfrd

DEFAULT_DT_FACTOR = 1e-5

def calculateBDDt( speciesList, factor = DEFAULT_DT_FACTOR ):

    D_list = []
    radius_list = []
    for species in speciesList:
        if species.pool.size != 0:
            D_list.append( species.D )
            radius_list.append( species.radius )
    D_max = max( D_list ) * 2  # max relative diffusion speed
    radius_min = min( radius_list )
    sigma_min = radius_min * 2

    dt = factor * sigma_min ** 2 / D_max  
    print 'bd dt = ', dt

    return dt

class BDSimulatorCoreBase( object ):
    
    '''
    BDSimulatorCore borrows the following from the main simulator:
    - speciesList
    - reactionTypes list (both 1 and 2)
    
    '''

    def __init__( self, main ):

        self.main = weakref.proxy( main )

        self.particleList = []

        self.t = 0.0
        self.dt = 0.0

        self.stepCounter = 0

        self.reactionEvents = 0

        self.speciesList = self.main.speciesList
        self.getReactionType1 = self.main.getReactionType1
        self.getReactionType2 = self.main.getReactionType2
        self.applyBoundary = self.main.applyBoundary

        self.populationChanged = False

        self.P_acct = {}

    def initialize( self ):
        self.determineDt()

    def clearParticleList( self ):
        self.particleList = []

    def addToParticleList( self, particle ):
        self.particleList.append( particle )

    def removeFromParticleList( self, particle ):
        self.particleList.remove( particle )


    def getNextTime( self ):
        return self.t + self.dt

    def stop( self, t ):
        # dummy
        self.t = t

    def determineDt( self ):

        self.dt = calculateBDDt( self.speciesList.values() )


    def getP_acct( self, rt, D, sigma ):

        try:
            return self.P_acct[ rt ]

        except KeyError:
            I = _gfrd.I_bd( sigma, self.dt, D )
            p = rt.k * self.dt / ( I * 4.0 * numpy.pi )
            if not 0.0 <= p < 1.0:
                raise RuntimeError,\
                    'Invalid acceptance ratio (%s) for reaction %s.' \
                    % ( p, rt )
            self.P_acct[ rt ] = p
            return p


    def step( self ):

        self.stepCounter += 1

        self.propagate()

        self.t += self.dt


    def propagate( self ):
        
        self.particlesToStep = self.particleList[:]

        #print particleList

        random.shuffle( self.particlesToStep )
        while len( self.particlesToStep ) != 0:
            particle = self.particlesToStep.pop() # last one
            self.propagateParticle( particle )


    def propagateParticle( self, particle ):

        species = particle.species
        rt1 = self.attemptSingleReactions( species )

        if rt1:
            try:
                self.fireReaction1( particle, rt1 )
            except NoSpace:
                print 'fireReaction1 rejected.'

            return


        D = species.D
        if D == 0.0:
            return

        displacement = drawR_free( self.dt, D )

        newpos = particle.pos + displacement
        
        closest, dist = self.getClosestParticle( newpos, 
                                                 ignore = [ particle, ] )
        
        if closest and dist <= species.radius:  # collision
            species2 = closest.species

            rt = self.main.reactionTypeMap2.get( ( species, species2 ) )
            k = rt.k

            if k != 0.0:
                radius12 = species.radius + species2.radius
                D12 = species.D + species2.D

                p = self.getP_acct( rt, D12, radius12 )

                rnd = numpy.random.uniform()
                if p > rnd:
                    print 'fire reaction2'
                    try:
                        self.fireReaction2( particle, closest, rt )
                    except NoSpace:
                        print 'fireReaction2 move rejected'
                    return
                else:
                    print 'reaction2 reject'
                    return

            else:
                print 'collision move reject'
                return


        self.applyBoundary( newpos )

        try:
            self.moveParticle( particle, newpos )
        except NoSpace:
            print 'propagation move rejected.'
        return


    def attemptSingleReactions( self, species ):

        reactionTypes = self.getReactionType1( species )
        if not reactionTypes:
            return None  # no reaction

        k_array = [ rt.k * self.dt for rt in reactionTypes ]
        k_array = numpy.add.accumulate( k_array )
        k_max = k_array[-1]

        rnd = numpy.random.uniform()
        if k_max < rnd:
            return None

        i = numpy.searchsorted( k_array, rnd )

        return reactionTypes[i]


    def fireReaction1( self, particle, rt ):
        
        oldpos = particle.pos.copy()

        if len( rt.products ) == 0:
            
            self.removeParticle( particle )
            
        elif len( rt.products ) == 1:
            
            productSpecies = rt.products[0]

            if not self.checkOverlap( oldpos, productSpecies.radius,
                                      ignore = [ particle, ] ):
                print 'no space for product particle.'
                raise NoSpace()
                
            self.removeParticle( particle )
            self.createParticle( productSpecies, oldpos )
            
        elif len( rt.products ) == 2:
            
            productSpecies1 = rt.products[0]
            productSpecies2 = rt.products[1]
            
            D1 = productSpecies1.D
            D2 = productSpecies2.D
            D12 = D1 + D2
            
            radius1 = productSpecies1.radius
            radius2 = productSpecies2.radius
            radius12 = radius1 + radius2

            for i in range( 100 ):

                rnd = numpy.random.uniform()
                pairDistance = drawR_gbd( rnd, radius12, self.dt, D12 )

                unitVector = randomUnitVector()
                vector = unitVector * pairDistance # * (1.0 + 1e-10) # safety
            
                # place particles according to the ratio D1:D2
                # this way, species with D=0 doesn't move.
                # FIXME: what if D1 == D2 == 0?
                newpos1 = oldpos + vector * ( D1 / D12 )
                newpos2 = oldpos - vector * ( D2 / D12 )
                #FIXME: check surfaces here
            
                self.applyBoundary( newpos1 )
                self.applyBoundary( newpos2 )

                # accept the new positions if there is enough space.
                if self.checkOverlap( newpos1, radius1,
                                      ignore = [ particle, ]) and \
                       self.checkOverlap( newpos2, radius2,
                                          ignore = [ particle, ]):
                    break
            else:
                print 'no space for product particles.'
                raise NoSpace()

            # move accepted
            self.removeParticle( particle )

            self.createParticle( productSpecies1, newpos1 )
            self.createParticle( productSpecies2, newpos2 )

        else:
            raise RuntimeError, 'num products >= 3 not supported.'

        self.reactionEvents += 1
        self.populationChanged = True



    def fireReaction2( self, particle1, particle2, rt ):

        pos1 = particle1.pos.copy()
        pos2 = particle2.pos.copy()

        if len( rt.products ) == 1:
                
            productSpecies = rt.products[0]

            D1 = particle1.species.D
            D2 = particle2.species.D

            pos2t = cyclicTranspose( pos2, pos1, self.main.worldSize )
            newPos = ( D2 * pos1 + D1 * pos2t ) / ( D1 + D2 )
            self.applyBoundary( newPos )

            self.clearVolume( newPos, productSpecies.radius,
                              ignore=[ particle1, particle2 ] )

            self.removeParticle( particle1 )
            self.removeParticle( particle2 )

            try:
                self.createParticle( productSpecies, newPos )
            except NoSpace:
                assert False, "this shouldn't happen"


            #print particle1, particle2, '->', newparticle


            try:
                self.particlesToStep.remove( particle2 )
            except ValueError:  
                pass     # particle2 already stepped, which is fine.

            self.reactionEvents += 1
            self.populationChanged = True
            return
        
        else:
            raise NotImplementedError,\
                'num products >= 2 not supported.'


    def check( self ):

        # particles don't overlap

        for particle in self.particleList:
            assert self.checkOverlap( particle.pos, particle.radius,
                                      ignore=[particle,] )



class BDSimulatorCore( BDSimulatorCoreBase ):
    

    def __init__( self, main ):

        BDSimulatorCoreBase.__init__( self, main )

        self.checkOverlap = self.main.checkOverlap

        self.moveParticle = self.main.moveParticle

        self.getNeighborParticles = main.getNeighborParticles
        self.getClosestParticle = main.getClosestParticle

        
    def initialize( self ):
        BDSimulatorCoreBase.initialize( self )

        self.updateParticleList()

    def updateParticleList( self ):

        self.clearParticleList()

        for s in self.speciesList.values():
            for i in range( s.pool.size ):
                self.addToParticleList( Particle( s, index = i ) )

    def addParticle( self, particle ):
        self.main.addParticle( particle )
        self.addToParticleList( particle )

    def removeParticle( self, particle ):
        self.main.removeParticle( particle )
        self.removeFromParticleList( particle )

    def createParticle( self, species, pos ):
        particle = self.main.createParticle( species, pos )
        self.addToParticleList( particle )



class BDSimulator( GFRDSimulatorBase ):
    
    def __init__( self ):

        GFRDSimulatorBase.__init__( self )

        self.core = BDSimulatorCore( self )
        self.isDirty = True

    def gett( self ):
        return self.core.t

    def sett( self, t ):
        self.core.t = t

    def getDt( self ):
        return self.core.dt

    def getStepCounter( self ):
        return self.core.stepCounter

    t = property( gett, sett )
    dt = property( getDt )
    stepCounter = property( getStepCounter )

    def initialize( self ):

        self.setAllRepulsive()

        self.core.initialize()

        self.isDirty = False



    def getNextTime( self ):
        return self.core.t + self.core.dt

    def stop( self, t ):
        # dummy
        self.core.stop( t )

    def step( self ):

        self.populationChanged = False

        if self.isDirty:
            self.initialize()

        if self.stepCounter % 10000 == 0:
            self.check()

        self.core.step()

        print self.stepCounter, ': t = ', self.t,\
            'dt = ', self.dt,\
            'reactions:', self.reactionEvents,\
            'rejected moves:', self.rejectedMoves
        #print ''



    def check( self ):
        pass
