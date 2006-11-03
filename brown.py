#!/usr/env python


import math
import bisect
import random
import sys
#from Numeric import *

import gfrdfunctions

#from scipy import gplt

import numpy
import scipy
import scipy.optimize


Pi = scipy.pi
Pi2 = scipy.pi * 2.0
PiSqrt = math.sqrt( scipy.pi )

N_A = 6.0221367e23
INF = numpy.Inf

NOWHERE = numpy.array( ( INF, INF, INF ) )


def MsTom3s( rate ):
    return rate / ( 1000 * N_A )

def distanceSq_Simple( position1, position2, fsize=None ):

    return ( ( position1 - position2 ) **2 ).sum()

def distanceSqArray_Simple( position1, positions, fsize=None ):
    
    return ( ( position1 - positions ) **2 ).sum(1)
    #return ( ( positions - position1 ) **2 ).sum(1)



def distanceSq_Cyclic( position1, position2, fsize ):

    halfsize = fsize * 0.5
    location = numpy.less( position1, halfsize ) * 2.0 - 1.0
    xtransposes = ( 0.0, location[0] * fsize )
    ytransposes = ( 0.0, location[1] * fsize )
    ztransposes = ( 0.0, location[2] * fsize )

    array = numpy.zeros( ( 8, 3 ), numpy.floating )

    i = 0
    for xtranspose in xtransposes:
        for ytranspose in ytransposes:
            for ztranspose in ztransposes:
                array[i] = ( ( position1 +\
                               ( xtranspose, ytranspose, ztranspose ) )\
                             - position2 ) **2
                i += 1

    return array.sum(1).min()


def distanceSqArray_Cyclic( position1, positions, fsize ):

    halfsize = fsize * 0.5

    location = numpy.less( position1, halfsize ) * 2.0 - 1.0
    xtransposes = ( 0.0, location[0] * fsize )
    ytransposes = ( 0.0, location[1] * fsize )
    ztransposes = ( 0.0, location[2] * fsize )

    array = numpy.zeros( ( 8, len(positions), 3 ), numpy.floating )

    i = 0
    for xtranspose in xtransposes:
        for ytranspose in ytransposes:
            for ztranspose in ztransposes:
                array[i] = ( ( position1 +\
                               ( xtranspose, ytranspose, ztranspose ) )\
                             - positions ) **2
                i += 1
                
    return array.sum(2).min(0)



def randomUnitVector():
    theta = numpy.random.uniform( 0, Pi2 )
    phi = numpy.random.uniform( 0, Pi2 )

    vector = numpy.array( [ math.cos( theta ) * math.sin( phi ),
                            math.sin( theta ) * math.sin( phi ),
                            math.cos( phi ) ] )
    return vector


class Species:
    
    def __init__( self, id, D, radius ):
        self.id = id
        self.D = D
        self.radius = radius
        self.pool = ParticlePool()

    def newParticle( self, position ):

        serial = self.pool.newParticle( position )


    def removeParticleByIndex( self, index ):
        self.pool.removeByIndex( index )

    def removeParticleBySerial( self, serial ):
        self.pool.removeBySerial( serial )


class ReactionType:

    def order( self ):
        return len( self.reactants )

    def dump( self ):
        s = ''
        for i in self.reactants:
            s += i.id
            s += ' '

        s += ' -> '

        for i in self.products:
            s += i.id
            s += ' '

        return s

class UnimolecularReactionType( ReactionType ):

    def __init__( self, s1, p1, k ):
        self.reactants = ( s1, )
        self.products = ( p1, )
        self.k  = k

class BindingReactionType( ReactionType ):

    def __init__( self, s1, s2, p1, k ):
        self.reactants = ( s1, s2 )
        self.products = ( p1, )
        self.k  = k

class UnbindingReactionType( ReactionType ):

    def __init__( self, s1, p1, p2, k ):
        self.reactants = ( s1, )
        self.products = ( p1, p2 )
        self.k  = k
    

class ParticlePool:

    def __init__( self ):

        self.indexMap = {}
        self.serialCounter = 0

        self.serials = numpy.array( [0,], numpy.integer )

        self.positions = numpy.array( [], numpy.floating )
        self.positions.shape = ( 0, 3 )

        self.size = 0

    def newParticle( self, position ):

        newindex = self.size
        newserial = self.serialCounter
        
        self.size += 1
        self.serialCounter += 1

        self.__resizeArrays( self.size )

        self.indexMap[ newserial ] = newindex
        self.serials[ newindex ] = newserial
        self.positions[ newindex ] = position

        return newserial
    

    def __resizeArrays( self, newsize ):

        self.serials = numpy.resize( self.serials, newsize )
        self.positions = numpy.resize( self.positions, ( newsize, 3 ) )

    def removeBySerial( self, serial ):

        index = self.getIndex( serial )
        self.removeByIndex( index )


    def removeByIndex( self, index ):

        self.size -= 1

        serialOfRemovedItem = self.serials[ index ]
        serialOfMovedItem = self.serials[ self.size ]

        # swap items at index and the end, and discard the item
        # to be removed, which is now at the end, by shrinking
        # the arrays by one.
        self.serials[ index ] = self.serials[ self.size ]
        self.positions[ index ] = self.positions[ self.size ]
        self.__resizeArrays( self.size )

        # book keeping
        del self.indexMap[ serialOfRemovedItem ]
        self.indexMap[ serialOfMovedItem ] = index

    def getSerialByIndex( self, index ):
        
        return self.serials[index]

    def getIndex( self, serial ):

        index = self.indexMap[ serial ]

        return index


class Reaction:

    def __init__( self ):
        pass

    


'''
class Particle:
    
    def __init__( self, pos, species ):
        self.species = species

        pool = species.pool
        
        self.serial = pool.newLot()

        index = pool.getIndexBySerial( self.serial )


        print 'idx', index
        print pool.positions
        pool.positions[ index ] = pos

        self.partner = None

        # temporary area, should be gotten rid of
        self.dt = array( (scipy.Inf,scipy.Inf) )

    def getPos( self ):
        pool = self.species.pool
        return pool.positions[ pool.getIndexBySerial( self.serial ) ]

    def setPos( self ):
        pool = self.species.pool
        return pool.positions[ pool.getIndexBySerial( self.serial ) ]
'''        



class Simulator:
    
    def __init__( self ):
        self.speciesList = {}
        self.reactionTypeList1 = {}
        self.reactionTypeList2 = {}

        self.dt = 1e-7
        self.t = 0.0

        self.H = 4.0
        
        self.dtLimit = 1e-3
        self.dtMax = self.dtLimit

        self.nextReaction = None

        self.fsize = 0.0

        self._newparticles = []

        #self._distanceSq = distanceSq_Simple
        #self._distanceSqArray = distanceSqArray_Simple
        self._distanceSq = distanceSq_Cyclic
        self._distanceSqArray = distanceSqArray_Cyclic

        self.resampledMoves = 0

    def simpleDiffusion( self, speciesIndex, particleIndex ):

        species = self.speciesList.values()[speciesIndex]

        limitSq = self.H * numpy.sqrt( 6.0 * species.D * self.dt ) - \
                  species.radius
        limitSq *= limitSq

        while True:
            displacement = gfrdfunctions.p1( species.D, self.dt )
            
            distanceSq = ( displacement * displacement ).sum()

            if distanceSq <= limitSq:
                break

            self.resampledMoves += 1

        pos = species.pool.positions[particleIndex]
        pos += displacement
        pos %= self.fsize


    def distanceSq( self, position1, position2 ):
        return self._distanceSq( position1, position2, self.fsize )

    def distance( self, position1, position2 ):
        return math.sqrt( self.distanceSq( position1, position2 ) )
        
    def distanceSqArray( self, position1, positions ):
        return self._distanceSqArray( position1, positions, self.fsize )

    def setSize( self, size ):
        self.fsize = size


    def addSpecies( self, species ):
        self.speciesList[ species.id ] = species

    def addReactionType( self, r ):

        numReactants = len( r.reactants )


        if numReactants == 1:
            species1 = r.reactants[0]
            self.reactionTypeList1[species1] = r
        elif numReactants == 2:
            species1 = r.reactants[0]
            species2 = r.reactants[1]
            self.reactionTypeList2[ (species1,species2) ] = r
            if species1 != species2:
                self.reactionTypeList2[ (species2,species1) ] = r
        else:
            raise 'unexpected'

    def throwInParticles( self, id, n ):
        print 'throwing in %s %s particles' % ( n, id )

        species = self.speciesList[ id ]
        
        for i in range( n ):

            while True:

                position= numpy.array( [ random.uniform( 0, self.fsize ),\
                                         random.uniform( 0, self.fsize ),\
                                         random.uniform( 0, self.fsize ) ] )
                if self.checkOverlap( position, species.radius ):
                    break
            
            species.newParticle( position )



    def placeParticle( self, id, position ):

        position = numpy.array(position)
        species = self.speciesList[ id ]

        if not self.checkOverlap( position, species.radius ):
            raise 'placeParticle: overlap check failed'
            
        species.newParticle( position )
        
        
    def checkOverlap( self, position, radius ):
        
        for species2 in self.speciesList.values():

            if species2.pool.size == 0:
                continue

            positions2 = species2.pool.positions

            radius12 = radius + species2.radius
            radius12sq = radius12 * radius12

            closestdistsq = self.distanceSqArray( position, positions2 ).min()

#            print distsq

            if closestdistsq <= radius12sq:
                print 'reject:', math.sqrt(closestdistsq), math.sqrt( radius12sq )
                return False

        return True
    

    
    def step( self ):
    
        self.clear()

        self.formPairs()


        # proceed slightly even if dtMax is 0.
        if self.dtMax <= 1e-18:
            self.dtMax = 1e-18
        
        self.determineNextReaction()

        print 'maxdt', self.dtMax, 'dt', self.dt,\
              'resampled moves', self.resampledMoves
        
        self.propagateParticles()

        self.t += self.dt

        self.newParticles()


    def isPopulationChanged( self ):
        return self.nextReaction != None


    def determineNextReaction( self ):

        self.dt = self.dtMax

        # first order reactions
        for rt in self.reactionTypeList1.values():
            reactantSpecies = rt.reactants[0]

            pool = reactantSpecies.pool

            dt, i = self.nextReactionTime1( rt, pool )
            
            reaction1 = ( dt, reactantSpecies, i, rt )

            if i != None:
            
                dt1 = reaction1[0]

                if self.dt >= dt1:
                    self.dt = dt1
                    self.nextReaction = reaction1

                

        # second order reactions

        self.pairs = [ ( self.nextReactionTime2( r ),\
                         r[1], r[2], r[3], r[4], r[5], r[6] )\
                       for r in self.pairs ]

        if len( self.pairs ) != 0:
            
            self.pairs.sort()
            reaction2 = self.pairs[0]
            dt2 = reaction2[0]

            if self.dt >= dt2:

                self.dt = dt2
                self.nextReaction = reaction2

        print 'next reaction = ', self.nextReaction


    def nextReactionTime1( self, rt, pool ):

        if pool.size == 0:
            return None, None

        dts = numpy.array( [ random.random() for i in range( pool.size ) ] )
        dts = - numpy.log( dts ) / rt.k

        i = numpy.argmin( dts )

        return dts[i], i 
        


    def nextReactionTime2( self, pair ):

        ( dt, ndt, s1, i1, s2, i2, rt ) = pair

        if rt == None:
            return INF

        species1 = self.speciesList.values()[s1]
        species2 = self.speciesList.values()[s2]
        pos1 = species1.pool.positions[i1]
        pos2 = species2.pool.positions[i2]

        radius = species1.radius + species2.radius
        r0 = math.sqrt( ( ( pos1 - pos2 ) ** 2 ).sum() )
        D = species1.D + species2.D
        k = rt.k
        

        if radius > r0:
            print 'CRITICAL: radius > r0', str(radius), str(r0)
            #            return scipy.Inf
            print pair
            print rt.dump()
            #sys.exit(-1)

        u = random.random()

#        infp = gfrdreaction.p_survival2_inf( radius, r0, D, k )
#        print 'sur', infp
#        if u <= infp:
#            return scipy.Inf
        
        maxtp = gfrdfunctions.p_survival2( self.dtMax, radius, r0, D, k )
#        print 'maxtp', maxtp
#        if k == 0.0:
#            print 'kmaxp', maxtp

        if u <= maxtp:
            return INF

        def __s( f, u, radius, r0, D, k ):
            print '__s', ( f,u,radius, r0, D, k ),\
                  u - gfrdfunctions.p_survival2( f, radius, r0, D, k )
            return u - gfrdfunctions.p_survival2( f, radius, r0, D, k )


        res = scipy.optimize.brenth( __s, 1e-100, self.dtMax,\
                                     args = ( u, radius, r0, D, k) )

        return res
        

    
    def fireReaction( self ):

        if self.nextReaction == None:
            return

        if len( self.nextReaction.reactants ) == 1:
            self.fireReaction1( self.nextReaction )
        else:
            self.fireReaction2( self.nextReaction )

    def fireReaction1( self, reaction ):

        dt, reactantSpecies, index, rt = reaction
        
        pos = reactantSpecies.pool.positions[index].copy()

        if len( rt.products ) == 1:
            productSpecies = rt.products[0]

            reactantSpecies.removeParticleByIndex( index )

            if reactantSpecies.radius < productSpecies.radius and \
                   not self.checkOverlap( pos, productSpecies1.radius ):
                raise 'fireReaction1: overlap check failed'

            productSpecies.newParticle( pos )

        elif len( rt.products ) == 2:
            
            productSpecies1 = rt.products[0]
            productSpecies2 = rt.products[1]

            D1 = productSpecies1.D
            D2 = productSpecies2.D

            reactantSpecies.removeParticleByIndex( index )

            unitVector = randomUnitVector()

            #print 'unit', self.distance( unitVector, numpy.array([0,0,0]) )
            distance = productSpecies1.radius + productSpecies2.radius
            vector = unitVector * distance * ( 1.0 + 1e-18 ) # safety
            
            newpos1 = pos + vector * ( D1 / (D1 + D2) )
            newpos2 = pos - vector * ( D2 / (D1 + D2) )
            
            newpos1 %= self.fsize
            newpos2 %= self.fsize
            
            productSpecies1.newParticle( newpos1 )
            productSpecies2.newParticle( newpos2 )

        elif len( rt.products ) == 0:

            reactantSpecies.removeParticleByIndex( index )

        else:
            raise "num products >= 3 not supported."

            
    def fireReaction2( self, pair ):

        print 'fire:', pair

        dt, ndt, speciesIndex1, index1, speciesIndex2, index2, rt = pair

        species1 = self.speciesList.values()[speciesIndex1]
        species2 = self.speciesList.values()[speciesIndex2]

        pos1 = species1.pool.positions[ index1 ].copy()
        pos2 = species2.pool.positions[ index2 ].copy()

        D1 = species1.D
        D2 = species2.D

        if len( rt.products ) == 1:

            species3 = rt.products[0]

            serial1 = species1.pool.getSerialByIndex( index1 )
            serial2 = species2.pool.getSerialByIndex( index2 )

            if D1 == 0.0:
                newpos = pos1
            elif D2 == 0.0:
                newpos = pos2
            else:
                
                D2D1Sqrt = math.sqrt( D2 / D1 ) 
                D1D2Sqrt = math.sqrt( D1 / D2 )
                
                R0 = D2D1Sqrt * pos1 + D1D2Sqrt * pos2
                
                dR = gfrdfunctions.p2_R( D1, D2, self.dt )
                newpos = ( R0 + dR ) / ( D1D2Sqrt + D2D1Sqrt )
                
                
            newpos %= self.fsize
                
            species1.removeParticleBySerial( serial1 )
            species2.removeParticleBySerial( serial2 )
            species3.newParticle( newpos )


            '''
            else: 
                species1.pool.positions[ index1 ] = pos1
                species2.pool.positions[ index2 ] = pos2
                try:  # no space for the product molecule, just propagate.
                    print 'fireReaction2: no space for product, just propagate'
                    self.propagatePairs( [ pair, ] )
                except:  #FIXME: catch a class
                    print 'fireReaction2: fail, restoring pair positions.'
                    # restore pair position (already done above)
'''

        
    def clear( self ):

        self.dtMax = self.dtLimit
        self.dt = self.dtLimit

        self.nextReaction = None

        self.pairs = []
        self.singles = []

#        for species in self.speciesList.values():
#            for p in species.pool.keys():
#                p.dt[0] = scipy.Inf
#                p.dt[1] = scipy.Inf
#                p.partner = None
                
    def debug( self ):

        for species in self.speciesList.values():
            #            print 'species: %s' % species.id

            for i in species.pool.keys():
                print i.species.id, i.getPos(), i.dt


    def propagateParticles( self ):

        self.propagateSingles()


        # fireReaction should come last in any case because
        # it can change particle identities.

        if self.nextReaction == None:  # no reaction
            self.propagatePairs( self.pairs )
        elif len( self.pairs ) != 0 and\
                 self.nextReaction == self.pairs[0]: # binary reaction
            if len( self.pairs ) > 1: 
                self.propagatePairs( self.pairs[1:] )
            self.fireReaction2( self.nextReaction )
        else:                           # unary reaction
            self.propagatePairs( self.pairs )
            self.fireReaction1( self.nextReaction ) # after propagatePairs()


    def propagatePairs( self, pairs ):

        if len( pairs ) == 0:
            return


        # dummy procedure, assumes no correlation.
        for pair in pairs:

            print pair
            dt, ndt, speciesIndex1, i1, speciesIndex2, i2, rt = pair

            species1 = self.speciesList.values()[speciesIndex1]
            species2 = self.speciesList.values()[speciesIndex2]
            

            if species1.D != 0.0:
                dts1 = self.dtCache[speciesIndex1][i1]
                pairradii1 = self.H * numpy.sqrt( 6.0 * species1.D * dts1 )
                self.simpleDiffusion( speciesIndex1, i1 )


            if species2.D != 0.0:
                dts2 = self.dtCache[speciesIndex2][i2]
                pairradii2 = self.H * numpy.sqrt( 6.0 * species2.D * dts2 )
                self.simpleDiffusion( speciesIndex2, i2 )

            
    def propagateSingles( self ):


        for i in range( len( self.singles ) ):
            species = self.speciesList.values()[i]

            if species.D != 0.0:
                for j in self.singles[i]:

                    self.simpleDiffusion( i, j )



    def newParticles( self ):

        for i in self._newparticles:
            i.species.insert( i )

        self._newparticles = []
        
    def formPairs( self ):

        # 1. form pairs in self.pairs
        # 2. list singles in self.singles

        speciesList = self.speciesList.values()

        # list up pair candidates

        # partner -> nearest particle
        # neighbor -> second nearest particle

        # dtTable[ speciesIndex ][ particleIndex ][ 0 .. 1 ]
        self.dtCache = []
        self.neighborCache = []
        #self.drCache = []
        for speciesIndex in range( len( speciesList ) ):
            size = speciesList[speciesIndex].pool.size
            self.dtCache.append( numpy.zeros( ( size, 2), numpy.floating ) )
            self.neighborCache.append( [[[ -1, -1 ],[-1,-1]]] * size )
            #self.drCache.append( numpy.zeros( size, numpy.floating ) )

        for speciesIndex in range( len( speciesList ) ):
            size = speciesList[speciesIndex].pool.size
            for particleIndex in range( size ):
                ret = self.checkPairs( speciesIndex, particleIndex )
                self.dtCache[ speciesIndex ][ particleIndex ] = ret[0]
                self.neighborCache[ speciesIndex ][ particleIndex ] = ret[1]
                #self.drCache[ speciesIndex ][ particleIndex ] = ret[2]


        self.pairs = []
        for speciesIndex1 in range( len( speciesList ) ):

            species1 = speciesList[speciesIndex1]
            positions1 = species1.pool.positions

            for particleIndex in range( len( positions1 ) ):

                # A partner: the other of the pair.
                # A neighbor of a pair: closer of the second closest of
                #                       the particles in the pair.
                #                       This is different from neighbors of
                #                       a particle.
                
                # (1) Find the closest particle (partner).
                partner = self.neighborCache[ speciesIndex1 ][ particleIndex ][0]

                # (skip the last particle)
                #if partner[0] == -1:
                # continue

                #print partner
                
                dts = self.dtCache[ speciesIndex1 ][ particleIndex ]
                ( partnerSpeciesIndex, partnerIndex ) = partner

                partnersPartner = self.neighborCache\
                                  [ partnerSpeciesIndex ][ partnerIndex ][0]
                partnerDts = self.dtCache\
                             [ partnerSpeciesIndex ][ partnerIndex ]

                # (2) If the partner's partner has to be this, else
                #     this combination isn't a pair.
                # (3) 'Neighbor' of this pair is the closer of
                #     this and the partner's second closest.
                #     We take the particle that has this neighbor.
                if ( partnersPartner[0] != speciesIndex1 or \
                     partnersPartner[1] != particleIndex  ) or \
                       partnerDts[1] < dts[1]:
                    continue
                
                # (4) Now we have a candidate pair.
                species2 = speciesList[partnerSpeciesIndex]
                rt = self.reactionTypeList2.get( ( species1, species2 ), None )
                pair = ( dts[0], dts[1], speciesIndex1, particleIndex,\
                         partnerSpeciesIndex, partnerIndex, rt )
                self.pairs.append( pair )

                # dtMax = the minimum neighbor dt of all pairs.
                self.dtMax = min( self.dtMax, dts[1] )



        # screening pairs

        self.pairs.sort()
        #print pairs
        #for i in pairs:
        #print i[0]

        checklist = []
        for i in range( len( speciesList ) ):
            checklist.append( numpy.ones( speciesList[i].pool.size ) )

        # screening 1
        for i in range( len( self.pairs ) ):
            ( dt, ndt, si1, i1, si2, i2, rt ) = self.pairs[i]


            # Don't take pairs with partner dt greater than dtMax.
            if dt > self.dtMax:
                self.pairs = self.pairs[:i]
                break

            # debug
            if checklist[si1][i1] == 0 or checklist[si2][i2] == 0:
                print self.pairs[:i+1]
                print self.dtCache[si1][i1], self.dtCache[si2][i2]
                print self.neighborCache[si1][i1], self.neighborCache[si2][i2]
                raise 'pairs not mutually exclusive. critical.'

            checklist[si1][i1] = 0
            checklist[si2][i2] = 0


        #print self.dtMax, self.pairs

        # now we have the final list of pairs

        #for i in range(len(checklist)):
        #    print checklist[i]
            

        # make the list of singles.
        # a single is a particle that is not in the list
        # of pairs.
        self.singles = []

        for i in range( len( checklist ) ):
            singleIndices = numpy.flatnonzero( checklist[i] )
            self.singles.append( singleIndices )

        #debug
        numSingles = 0
        for i in self.singles:
            numSingles += len(i)

        print '# pairs = ', len(self.pairs),\
              ', # singles = ', numSingles



    def checkPairs( self, speciesIndex1, particleIndex ):

        speciesList = self.speciesList.values()

        neighbordts = [ INF, INF ]
        #neighbordts = []
        # [ ( index, reaction type ), ]
        neighbors = [ ( -1, -1 ), ( -1, -1 ) ]
        #neighbors = []
        drSqs = []

        species1 = speciesList[ speciesIndex1 ]
        positions = species1.pool.positions
        position1 = positions[ particleIndex ].copy()

        if len( position1 ) >= 2 and species1.D != 0.0:

            # temporarily displace the particle
            positions[particleIndex] = NOWHERE
            
            topDts, topIndices = self.checkDistance( position1, positions,
                                                     species1, species1 )

            # restore the particle.
            positions[particleIndex] = position1

            neighbordts.extend( topDts )

            if len( topIndices ) == 2:
                neighbors.extend( ( ( speciesIndex1, topIndices[0] ),\
                                    ( speciesIndex1, topIndices[1] ) ) )
            else:
                neighbors.extend( ( ( speciesIndex1, topIndices[0] ), ) )

        #for speciesIndex2 in range( speciesIndex1 + 1, len( speciesList ) ):
        for speciesIndex2 in range( speciesIndex1 )\
                + range( speciesIndex1 + 1, len( speciesList ) ):
            species2 = speciesList[speciesIndex2]

            if species2.pool.size == 0:
                continue

            if species1.D + species2.D == 0.0:
                continue
                    
            positions = species2.pool.positions

            if species2.pool.size == 1:  # insert a dummy
                positions = numpy.concatenate( ( positions, [NOWHERE,] ) )
                
            topDts, topIndices = self.checkDistance( position1, positions,
                                                     species1, species2 )
            neighbordts.extend( topDts )
            neighbors.extend( ( ( speciesIndex2, topIndices[0] ),\
                                ( speciesIndex2, topIndices[1] ) ) )

        topargs = numpy.argsort( neighbordts )[:2]
        topNeighborDts = numpy.take( neighbordts, topargs )
        topNeighbors = numpy.take( neighbors, topargs, 0  )

        #topDrSq = numpy.amin( drSqs )
        return topNeighborDts, topNeighbors #, topDrSq


    def checkDistance( self, position1, positions2, species1, species2 ):

        #positions2 = species2.pool.positions

        distanceSq = self.distanceSqArray( position1, positions2 )
        sortedindices = distanceSq.argsort()

        #print sortedindices
        #debug
        radius12 = species1.radius + species2.radius
        radius12sq = radius12 * radius12

        # check if particles overlap
        if distanceSq[ sortedindices[0] ] < radius12sq - 1e-20 and \
               distanceSq[ sortedindices[0] ] != 0.0:
            print position1, positions2[ sortedindices[0] ]
            print 'dr<radius', math.sqrt(distanceSq[sortedindices[0]]), radius12
            print species1.id, species2.id, sortedindices[0]
            raise "critical"

        factor = ( math.sqrt( species1.D ) + math.sqrt( species2.D ) ) * self.H
        factor *= factor * 6.0

        distanceSqSorted = distanceSq.take( sortedindices[:2] )

        # instead of just
        dts = distanceSqSorted / factor
        # below takes into account of particle radii.
        #dts = numpy.sqrt( distanceSqSorted ) - radius12
        #dts *= dts
        #dts /= factor

        #if min( dts ) < 0.0:
        #print 'negative dts occured, clipping.', dts
        #dts = numpy.clip( dts, 0.0, INF )
        # raise 'stop'

        indices = sortedindices[:2]

        return dts, indices # , distanceSqSorted[0]
