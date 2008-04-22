#!/usr/env python

import weakref

import math

import numpy
#import scipy
#import scipy.optimize


from utils import *
from surface import *

from ObjectMatrix import *

from gfrdbase import *
from bd import *


SAFETY = 1.0 + 1e-5


class Delegate( object ):

    def __init__( self, obj, method ):
        self.obj = weakref.proxy( obj )
        self.method = method


    def __call__( self, arg ):
        return self.method( self.obj, arg )


class Shell( object ):
    def __init__( self, pos, radius ):
        self.pos = pos.copy()
        self.radius = radius


'''
Used internally by Multi.
'''

class MultiBDCore( BDSimulatorCoreBase ):
    
    def __init__( self, main, multi ):

        BDSimulatorCoreBase.__init__( self, main )

        self.multi = multi

        self.particleMatrix = SimpleObjectMatrix()
        self.particleMatrix.setWorldSize( self.main.worldSize )

        self.shellMatrix = SimpleObjectMatrix()
        self.shellMatrix.setWorldSize( self.main.worldSize )

        self.escaped = False

        #self.checkOverlap is implemented below.
        #self.getNeighborParticles() is implemented below.
        #self.getClosestParticle() is implemented below.

        
    def updateParticle( self, particle, pos ):
        self.particleMatrix.update( ( particle.species, particle.serial ), pos,
                                    particle.radius )

    def initialize( self ):
        BDSimulatorCoreBase.initialize( self )

        self.updateShellMatrix()


    def step( self ):
        self.escaped = False
        
        BDSimulatorCoreBase.step( self )


    def updateShellMatrix( self ):

        self.shellMatrix.clear()
        for shell in self.multi.shellList:
            self.shellMatrix.add( shell, shell.pos, shell.radius )

    def _addParticle( self, particle ):
        self.addToParticleList( particle )
        self.particleMatrix.add( ( particle.species, particle.serial ),
                                 particle.pos, particle.radius )

    def addParticle( self, particle ):
        self._addParticle( particle )

    def removeParticle( self, particle ):
        self.main.removeParticle( particle )
        self.removeFromParticleList( particle )
        self.particleMatrix.remove( ( particle.species, particle.serial ) )


    def createParticle( self, species, pos ):

        if self.withinShell( pos, species.radius ):
            if not self.checkOverlap(pos, species.radius, ignore ):
                raise NoSpace()
        else:
            self.escaped = True
            self.clearOuterVolume( pos, species.radius )

            
        particle = self.main.createParticle( species, pos )
        self._addParticle( particle )

        return particle


    def moveParticle( self, particle, pos ):

        if self.withinShell( pos, particle.radius ):
            if not self.checkOverlap( pos, particle.radius, 
                                      ignore=[particle] ):
                raise NoSpace()
        else:
            self.escaped = True
            self.clearOuterVolume( pos, particle.radius, ignore=[particle] )

        self.main.moveParticle( particle, pos )
        self.updateParticle( particle, pos )



        
    def clearVolume( self, pos, radius, ignore=[] ):

        if not self.checkOverlap(pos, radius, ignore ):
            raise NoSpace()

        if not self.withinShell( pos, radius ):
            self.clearOuterVolume( pos, radius, ignore )

    def clearOuterVolume( self, pos, radius, ignore=[] ):

        self.main.clearVolume( pos, radius, ignore=[self.multi,] )
        if not self.main.checkOverlap( pos, radius, ignore ):
            raise NoSpace()

 
    def withinShell( self, pos, radius ):
        n, _ = self.shellMatrix.getNeighborsWithinRadius( pos, - radius )

        if n:
            return True
        else:
            return False

        
    def checkOverlap( self, pos, radius, ignore=[] ):
        
        n, _ = self.particleMatrix.getNeighborsWithinRadius( pos, radius )

        for particle in ignore:
            n.remove( ( particle.species, particle.serial ) )

        if n:
            return False
        else:
            return True


    def getNeighborParticles( self, pos, n=None, dummy=None ):

        n, d = self.particleMatrix.getNeighbors( pos, n, dummy )
        neighbors = [ Particle( i[0], i[1] ) for i in n ]
        return neighbors, d


    def getClosestParticle( self, pos, ignore=[] ):

        neighbors, distances =\
            self.getNeighborParticles( pos, 
                                       len( ignore ) + 1,
                                       dummy = ( None, -1 ) )

        for i, neighbor in enumerate( neighbors ):
            if neighbor not in ignore:
                return neighbor, distances[i]

        # default case: none left.
        return None, INF
        #return DummyParticle(), INF


    def check( self ):

        BDSimulatorCoreBase.check( self )

        # shellMatrix consistency
        for shell in self.multi.shellList:
            pos, radius = self.shellMatrix.get( shell )
            assert not ( pos - shell.pos ).any()
            assert not ( radius - shell.radius ).any()


        # shells are contiguous
        for shell in self.multi.shellList:
            n, d = self.shellMatrix.getNeighbors( shell.pos )
            assert d[1] - shell.radius < 0.0, 'shells are not contiguous.'

        # all particles within the shell.
                
        for p in self.particleList:
            assert self.withinShell( p.pos, p.species.radius ),\
                'not all particles within the shell.'



class Single( object ):

    def __init__( self, particle, reactiontypes ):

        self.multiplicity = 1

        self.particle = particle
        self.reactiontypes = reactiontypes

        self.k_tot = 0

        self.lastTime = 0.0
        self.dt = 0.0
        self.eventType = None

        self.shellList = [ Shell( self.particle.pos, self.getMinRadius() ), ]

        self.eventID = None

        self.gf = FirstPassageGreensFunction( particle.species.D )

        self.updatek_tot()


    def getD( self ):

        return self.particle.species.D

    def getPos( self ):

        return self.shellList[0].pos

    def setPos( self, pos ):
        self.shellList[0].pos = pos
        self.particle.pos = pos

    pos = property( getPos, setPos )


    def setRadius( self, radius ):

        assert radius >= self.getMinRadius()

        self.shellList[0].radius = radius


    '''
    A radius of a Single is the distance from the current position
    of the particle to the shell of the Single.   The shell defines the
    farthest point in space that it can occupy when it makes the maximum
    displacement.
    '''

    def getRadius( self ):

        return self.shellList[0].radius

    radius = property( getRadius, setRadius )


    def getMinRadius( self ):

        return self.particle.species.radius



    '''
    Initialize this Single.

    The radius (shell size) is shrunken to the radius of the particle
    it represents.   
    self.lastTime is reset to the current time, and self.dt
    is set to zero.
    '''

    def initialize( self, t ):

        self.reset()
        self.lastTime = t



    def calculateShellSize( self, closest, distance, shellDistance ):

        minRadius1 = self.getMinRadius()
        D1 = self.getD()

        if D1 == 0:
            return minRadius1

        D2 = closest.getD()
        minRadius2 = closest.getMinRadius()
        minRadius12 = minRadius1 + minRadius2
        sqrtD1 = math.sqrt( D1 )
            
        shellSize = min( sqrtD1 / ( sqrtD1 + math.sqrt( D2 ) )
                         * ( distance - minRadius12 ) + minRadius1,
                         shellDistance )
        shellSize /= SAFETY
        shellSize = max( shellSize, minRadius1 ) # not smaller than the radius

        return shellSize

        


    '''
    A mobility radius indicates the maximum displacement this single
    particle can make.

    Mobility radius of a particle is calculated as follows;

    mobility radius = Single radius - Particle radius.

    '''
    
    def getMobilityRadius( self ):

        #return self.radius - ( self.getMinRadius() * 2 )
        return self.radius - self.getMinRadius()




    def drawDisplacement( self, r ):

        rnd = numpy.random.uniform( size=2 )

        displacementS = [ r, rnd[0] * Pi, rnd[1] * 2 * Pi ]
        displacement = sphericalToCartesian( displacementS )

        return displacement


    '''
    Reset the Single.

    Radius (shell size) is shrunken to the actual radius of the particle.
    self.dt is reset to 0.0.  Do not forget to reschedule this Single
    after calling this method.
    '''

    def reset( self ):

        self.setRadius( self.getMinRadius() )
        self.dt = 0.0
        self.eventType = EventType.ESCAPE


    def isReset( self ):

        return self.radius == self.getMinRadius() and self.dt == 0.0\
               and self.eventType == EventType.ESCAPE
        

    def drawR( self, dt ):
        assert dt >= 0

        self.gf.seta( self.getMobilityRadius() )
        rnd = numpy.random.uniform()
        r = self.gf.drawR( rnd , dt )
        return r



    def determineNextEvent( self ):
        if self.getD() == 0:
            firstPassageTime = INF
        else:
            firstPassageTime = self.drawEscapeTime()
            
        reactionTime = self.drawReactionTime()

        if firstPassageTime <= reactionTime:
            self.dt = firstPassageTime
            self.eventType = EventType.ESCAPE
        else:
            self.dt = reactionTime
            self.eventType = EventType.REACTION


    def drawReactionTime( self ):
        
        if self.k_tot == 0:
            return INF

        rnd = numpy.random.uniform()
        dt = ( 1.0 / self.k_tot ) * math.log( 1.0 / rnd )

        return dt


    def drawEscapeTime( self ):
        
        rnd = numpy.random.uniform()
        self.gf.seta( self.getMobilityRadius() )
        dt = self.gf.drawTime( rnd )
        return dt


    def updatek_tot( self ):

        self.k_tot = 0

        if not self.reactiontypes:
            return

        for rt in self.reactiontypes:
            self.k_tot += rt.k


    def drawReactionType( self ):

        k_array = [ rt.k for rt in self.reactiontypes ]
        k_array = numpy.add.accumulate( k_array )
        k_max = k_array[-1]

        rnd = numpy.random.uniform()
        i = numpy.searchsorted( k_array, rnd * k_max )

        return self.reactiontypes[i]


    def check( self ):
        pass

    def __str__( self ):
        return 'Single' + str( self.particle )



'''
Just a free func ver of Pair.getCoM().
'''

def calculatePairCoM( pos1, pos2, D1, D2, worldSize ):

    pos2t = cyclicTranspose( pos2, pos1, worldSize )

    return ( ( D2 * pos1 + D1 * pos2t ) / ( D1 + D2 ) ) % worldSize


class Pair( object ):
    
    # CUTOFF_FACTOR is a threshold to choose between the real and approximate
    # Green's functions.
    # H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7,
    # 5.6: ~1e-8, 6.0: ~1e-9
    CUTOFF_FACTOR = 5.6

    def __init__( self, single1, single2, rt, distFunc, worldSize ):

        self.multiplicity = 2

        # Order single1 and single2 so that D1 < D2.
        if single1.particle.species.D <= single2.particle.species.D:
            self.single1, self.single2 = single1, single2 
        else:
            self.single1, self.single2 = single2, single1 

        self.rt = rt

        self.distance = distFunc
        self.worldSize = worldSize
        
        particle1 = self.single1.particle
        particle2 = self.single2.particle

        self.D1, self.D2 = particle1.species.D, particle2.species.D

        self.D_tot = self.D1 + self.D2
        self.D_geom = math.sqrt( self.D1 * self.D2 )  # geometric mean

        #self.minRadius = max( particle1.species.radius,
        #                      particle2.species.radius )
        self.sigma = particle1.species.radius + particle2.species.radius

        self.sgf = FirstPassageGreensFunction( self.D_geom )
        self.pgf = FirstPassagePairGreensFunction( self.D_tot, 
                                                   rt.k, self.sigma )
        self.pgf_free = FreePairGreensFunction( self.D_tot )
        self.pgf_basic = BasicPairGreensFunction( self.D_tot, rt.k, self.sigma )
        self.pgf_nocol = FirstPassageNoCollisionPairGreensFunction( self.D_tot )

        self.eventID = None

        self.lastTime = 0.0
        self.dt = 0.0
        self.eventType = None

        self.shellList = [ Shell( self.getCoM(), self.getMinRadius() ), ]



    #def __del__( self ):
    #pass
    #        print 'del', str( self )

    def initialize( self, t ):

        self.lastTime = t
        self.radius = self.getMinRadius()
        self.dt = 0
        self.eventType = None

    def setPos( self, pos ):
        self.shellList[0].pos = pos

    def getPos( self ):

        return self.shellList[0].pos

    pos = property( getPos, setPos )

    def getD( self ):

        return self.D_tot #FIXME: is this correct?

    def setRadius( self, radius ):

        assert radius >= self.getMinRadius()
        self.shellList[0].radius = radius


    def getRadius( self ):

        return self.shellList[0].radius

    radius = property( getRadius, setRadius )


    '''
    This method returns the radius from its CoM that this Pair must reserve
    to remain mobile.
    '''

    def getMinRadius( self ):

        pairDistance = self.distance( self.single1.pos,
                                      self.single2.pos )
        minRadius = max( pairDistance * self.D1 /
                         self.D_tot + self.single1.getMinRadius(),
                         pairDistance * self.D2 /
                         self.D_tot + self.single2.getMinRadius() )
        return minRadius


    '''
    Calculate and return the "Center of Mass" (== CoM) of this pair.
    '''

    def getCoM( self ):

        particle1 = self.single1.particle
        particle2 = self.single2.particle
        
        pos1 = particle1.pos
        pos2 = particle2.pos

        pos2t = cyclicTranspose( pos2, pos1, self.worldSize ) #FIXME:
        
        com = ( pos1 * self.D2 + pos2t * self.D1 ) / self.D_tot

        return com % self.worldSize


    def choosePairGreensFunction( self, r0, t ):

        distanceFromSigma = r0 - self.sigma
        distanceFromShell = self.a_r - r0;

        thresholdDistance = Pair.CUTOFF_FACTOR * \
            math.sqrt( 6.0 * self.D_tot * t )

        if distanceFromSigma < thresholdDistance:
        
            if distanceFromShell < thresholdDistance:
                # near both a and sigma;
                # use FirstPassagePairGreensFunction
                print 'normal'
                return self.pgf
            else:
                # near sigma; use BasicPairGreensFunction
                print 'near only sigma'
                return self.pgf_basic
                #return self.pgf
        else:
            if distanceFromShell < thresholdDistance:
                # near a;
                print 'near only a'
                return self.pgf_nocol
                
            else:
                # distant from both a and sigma; 
                print 'free'
                return self.pgf_free


    '''
    Calculate new positions of the pair particles using
    a new center-of-mass, a new inter-particle vector, and
    an old inter-particle vector.

    '''

    def newPositions( self, CoM, newInterParticle, oldInterParticle ):

        #FIXME: needs better handling of angles near zero and pi.

        # I rotate the new interparticle vector along the
        # rotation axis that is perpendicular to both the
        # z-axis and the original interparticle vector for
        # the angle between these.
        
        # the rotation axis is a normalized cross product of
        # the z-axis and the original vector.
        # rotationAxis = crossproduct( [ 0,0,1 ], interParticle )

        angle = vectorAngleAgainstZAxis( oldInterParticle )
        if angle % numpy.pi != 0.0:
            rotationAxis = crossproductAgainstZAxis( oldInterParticle )
            rotationAxis = normalize( rotationAxis )
            rotated = rotateVector( newInterParticle,
                                    rotationAxis,
                                    angle )
        elif angle == 0.0:
            rotated = newInterParticle
        else:
            rotated = numpy.array( [ newInterParticle[0], newInterParticle[1],
                                     - newInterParticle[2] ] )
            #rotated = newInterParticle * numpy.array( [ 1, 1, -1 ] )

        newpos1 = CoM - rotated * ( self.D1 / self.D_tot )
        newpos2 = CoM + rotated * ( self.D2 / self.D_tot )

        return newpos1, newpos2
        

    def determineNextEvent( self ):

        particle1 = self.single1.particle
        particle2 = self.single2.particle

        species1 = particle1.species
        species2 = particle2.species
        radius1 = species1.radius
        radius2 = species2.radius

        pos1 = particle1.pos
        pos2 = particle2.pos

        D1 = self.D1
        D2 = self.D2

        D1_factor = D1 / self.D_tot
        D2_factor = D2 / self.D_tot

        shellSize = self.radius

        sqrtD_tot = math.sqrt( self.D_tot )
        sqrtD_geom = math.sqrt( self.D_geom )

        r0 = self.distance( pos1, pos2 )

        assert r0 >= self.sigma, \
            '%s;  r0 %g < sigma %g' % ( self, r0, self.sigma )

        # equalize expected mean t_r and t_R.

        r0_1 = r0 * D1_factor
        r0_2 = r0 * D2_factor

        D_factor = sqrtD_tot + sqrtD_geom

        qrrtD1D25 = ( D1    * D2**5 ) ** 0.25
        qrrtD15D2 = ( D1**5 * D2 ) ** 0.25

        if qrrtD15D2 * r0 + ( qrrtD15D2 + qrrtD1D25 ) * radius1 \
                + D1 * ( sqrtD_tot * ( shellSize - radius2 ) 
                         - sqrtD_geom * radius2 )\
                - D2 * ( sqrtD_geom * r0 + sqrtD_tot * 
                         ( shellSize - radius1 ) )\
                         - qrrtD1D25 * radius2 >= 0:

            den1 = qrrtD1D25 + D1 * ( sqrtD_geom + sqrtD_tot )

            a_R_1 = sqrtD_geom * ( D2 * ( shellSize - radius1) + 
                                   D1 * ( shellSize - r0 - radius1 ) ) / den1

            a_r_1 = self.D_tot * ( sqrtD_geom * r0 + sqrtD_tot * 
                                   ( shellSize - radius1 ) ) / den1

            assert a_R_1 + a_r_1 * D1_factor + radius1 >= \
                a_R_1 + a_r_1 * D2_factor + radius2

            assert abs( a_R_1 + a_r_1 * D1_factor + radius1 - shellSize ) \
                < 1e-12 * shellSize

            self.a_r = a_r_1
            self.a_R = a_R_1
        else:
            den2 = qrrtD15D2 + D2 * ( sqrtD_geom + sqrtD_tot )

            a_R_2 = sqrtD_geom * ( D1 * ( shellSize - radius2 ) + 
                                   D2 * ( shellSize - r0 - radius2 ) ) / den2

            a_r_2 = self.D_tot * ( sqrtD_geom * r0 + sqrtD_tot * 
                                   ( shellSize - radius2 ) ) / den2

            assert a_R_2 + a_r_2 * D2_factor + radius2 >= \
                a_R_2 + a_r_2 * D1_factor + radius1

            assert abs( a_R_2 + a_r_2 * D2_factor + radius2 - shellSize ) \
                < 1e-12 * shellSize

            self.a_r = a_r_2
            self.a_R = a_R_2

        #print 'r R', self.a_r, self.a_R
        #print 'tr, tR', (( self.a_r - r0 ) / math.sqrt(6 * self.D_tot))**2,\
        #      (self.a_R / math.sqrt( 6*self.D_geom ))**2

        #print 'a a_r a_R', shellSize, self.a_r, self.a_R
        assert self.a_r > 0
        assert self.a_R > 0 or ( self.a_R == 0 and ( D1 == 0 or D2 == 0 ) )
        assert self.a_r > r0, '%g %g' % ( self.a_r, r0 )

        rnd = numpy.random.uniform( size=3 )

        # draw t_R
        self.sgf.seta( self.a_R )
        self.t_R = self.sgf.drawTime( rnd[0] )

        # draw t_r
        self.pgf.seta( self.a_r )
        self.t_r = self.pgf.drawTime( rnd[1], r0 )

        # draw t_reaction
        t_reaction1 = self.single1.drawReactionTime()
        t_reaction2 = self.single2.drawReactionTime()

        if t_reaction1 < t_reaction2:
            self.t_single_reaction = t_reaction1
            self.reactingsingle = self.single1
        else:
            self.t_single_reaction = t_reaction2
            self.reactingsingle = self.single2

        #print ( self.t_R, self.t_r, self.t_single_reaction )
        self.dt = min( self.t_R, self.t_r, self.t_single_reaction )

        assert self.dt >= 0
        print 'dt ', self.dt, 't_R', self.t_R, 't_r', self.t_r
        #print self.pgf.dump()
        if self.dt == self.t_r:  # type = 0 (REACTION) or 1 (ESCAPE_r)
            self.eventType = self.pgf.drawEventType( rnd[2],
                                                     r0, self.t_r )
        elif self.dt == self.t_R: # type = ESCAPE_R (2)
            self.eventType = 2
        elif self.dt == self.t_single_reaction:  # type = single reaction (3)
            self.eventType = 3 
        else:
            raise 'never get here'

        #assert False




    def drawR_single( self, t, a ):

        self.sgf.seta( a )
        r = self.sgf.drawR( numpy.random.uniform(), t )
        while r > self.a_R: # redraw; shouldn't happen often
            print 'drawR_single: redraw'
            r = self.sgf.drawR( numpy.random.uniform(), t )

        return r


    '''
    Draw r for the pair inter-particle vector.
    '''
    def drawR_pair( self, r0, t, a ):

        gf = self.choosePairGreensFunction( r0, t )

        if hasattr( gf, 'seta' ):  # FIXME: not clean
            gf.seta( a )

        r = gf.drawR( numpy.random.uniform(), r0, t )

        while r >= self.a_r or r <= self.sigma: # redraw; shouldn't happen often
            print 'drawR_pair: redraw'
            self.sim.rejectedMoves += 1
            r = gf.drawR( numpy.random.uniform(), r0, t )

        return r


    '''
    Draw theta for the pair inter-particle vector.
    '''
    def drawTheta_pair( self, rnd, r, r0, t ):

        gf = self.choosePairGreensFunction( r0, t )
        theta = gf.drawTheta( rnd, r, r0, t )

        return theta


    def checkNewpos( self, pos1, pos2, com ):

        species1 = self.single1.particle.species
        species2 = self.single2.particle.species

        oldCoM = com
        
        # debug: check if the new positions are valid:
        newDistance = distance( pos1, pos2 )
        particleRadius12 = species1.radius + species2.radius

        # check 1: particles don't overlap.
        if newDistance <= particleRadius12:
            print 'rejected move: ', 'radii, interp',\
                  species1.radius + species2.radius, newDistance
            print 'DEBUG: dt, pos1, pos2, pos1, pos2',\
                  self.dt, pos1, pos2, pos1, pos2
            raise RuntimeError, 'New particles overlap'

        # check 2: particles within mobility radius.
        if self.distance( oldCoM, pos1 ) + species1.radius \
               > self.radius or \
               self.distance( oldCoM, pos2 ) + species2.radius \
               > self.radius:
            raise RuntimeError, 'New particle(s) out of protective sphere.'


    def check( self ):
        pass

    def __str__( self ):
        buf = 'Pair( ' + str(self.single1.particle) +\
              ', ' + str(self.single2.particle) + ' )'

        return buf


class Multi( object ):

    def __init__( self, main ):

        self.shellList = []
        self.eventID = None

        self.sim = MultiBDCore( main, self )


    def initialize( self, t ):

        self.lastTime = t
        self.startTime = t

        self.sim.initialize() # ??


    def getDt( self ):
        return self.sim.dt

    dt = property( getDt )


    def getMultiplicity( self ):
        return len( self.sim.particleList )

    multiplicity = property( getMultiplicity )

    def addParticle( self, particle ):
        self.sim.addParticle( particle )

    def addShell( self, pos, radius ):
        self.shellList.append( Shell( pos, radius ) )


    def check( self ):

        self.sim.check()


    def __str__( self ):

        if len( self.sim.particleList ) == 0:
            return 'Multi()'

        buf = 'Multi( '
        buf += str( self.sim.particleList[0] )
        for particle in self.sim.particleList[1:]:
            buf += ', ' + str( particle )
        buf += ' )'

        return buf


class DummySingle( object ):
    def __init__( self ):
        self.multiplicity = 1

        self.radius = 0.0

    def getMinRadius( self ):
        return 0.0

    def getD( self ):
        return 0.0

    def getPos( self ):
        return NOWHERE

    pos = property( getPos )




class EGFRDSimulator( GFRDSimulatorBase ):
    
    def __init__( self ):


        self.MULTI_SHELL_FACTOR = 0.1
        self.SINGLE_SHELL_FACTOR = 0.1

        #self.shellMatrix = ObjectMatrix()
        self.shellMatrix = SimpleObjectMatrix()

        GFRDSimulatorBase.__init__( self )

        self.isDirty = True
        self.scheduler = EventScheduler()

        self.smallT = 1e-8  # FIXME: is this ok?

        self.maxShellSize = INF

        self.reset()

    def setWorldSize( self, size ):

        if isinstance( size, list ) or isinstance( size, tuple ):
            size = numpy.array( size )

        GFRDSimulatorBase.setWorldSize( self, size )
        self.shellMatrix.setWorldSize( size )

    def setMatrixSize( self, size ):
        GFRDSimulatorBase.setMatrixSize( self, size )
        self.shellMatrix.setMatrixSize( size )

    def getMatrixCellSize( self ):

        return self.shellMatrix.cellSize

    def getNextTime( self ):
        return self.scheduler.getTopTime()

    def setMaxShellSize( self, maxShellSize ):

        self.maxShellSize = maxShellSize

    def getMaxShellSize( self ):

        return self.maxShellSize

    def reset( self ):

        self.t = 0.0
        self.dt = INF
        self.stepCounter = 0
        self.zeroSteps = 0
        self.rejectedMoves = 0
        self.reactionEvents = 0
        self.lastEvent = None
        self.populationChanged = False

        self.isDirty = True
        #self.initialize()
        

    def initialize( self ):

        GFRDSimulatorBase.initialize( self )

        self.setAllRepulsive()

        self.scheduler.clear()
        self.shellMatrix.clear()

        for species in self.speciesList.values():
            for i in range( species.pool.size ):
                particle = Particle( species, index=i )
                single = self.createSingle( particle )
                self.addToShellMatrix( single )
                self.addSingleEvent( single )

        self.isDirty = False


    def stop( self, t ):

        print 'stop at', t

        if self.t == t:
            return

        if t >= self.scheduler.getTopEvent().getTime():
            raise RuntimeError, 'Stop time <= next event time.'

        self.t = t
        
        scheduler = self.scheduler
        
        nonSingleList = []

        # first burst all Singles.
        for i in range( scheduler.getSize() ):
            obj = scheduler.getEventByIndex(i).getArg()
            if isinstance( obj, Pair ) or isinstance( obj, Multi ):
                nonSingleList.append( obj )
            elif isinstance( obj, Single ):
                self.burstSingle( obj )
            else:
                assert False, 'do not reach here'


        # then burst all Pairs.
        for obj in nonSingleList:
            self.burstObj( obj )

        self.dt = 0.0


    def step( self ):

        self.populationChanged = False

        if self.isDirty:
            self.initialize()

        #if self.stepCounter % 100 == 0:
        #self.check()

        self.stepCounter += 1

        event = self.scheduler.getTopEvent()
        self.t, self.lastEvent = event.getTime(), event.getArg()

        print self.stepCounter, ': t = ', self.t, ': event = ', self.lastEvent
        
        self.scheduler.step()

        nextEvent = self.scheduler.getTopEvent()
        nextTime, nextEventObject = nextEvent.getTime(), nextEvent.getArg()
        self.dt = nextTime - self.t


        # assert if not too many successive dt=0 steps occur.
        if self.dt == 0:
            self.zeroSteps += 1
        else:
            self.zeroSteps = 0
        assert self.zeroSteps < self.scheduler.getSize() * 2,\
            'too many dt=zero steps.  simulator halted?'


        assert self.scheduler.getSize() != 0

        print 'next dt = ', self.dt, 'reactions', self.reactionEvents,\
              'rejected moves', self.rejectedMoves
        assert self.scheduler.check()
        print ''



    def createSingle( self, particle ):

        rt = self.getReactionType1( particle.species )
        single = Single( particle, rt )
        single.initialize( self.t )

        return single


    def createPair( self, single1, single2 ):

        assert single1.dt == 0.0
        assert single2.dt == 0.0
        assert single1.getMobilityRadius() == 0.0
        assert single2.getMobilityRadius() == 0.0

        species1 = single1.particle.species
        species2 = single2.particle.species
        rt = self.reactionTypeMap2.get( ( species1, species2 ) )

        pair = Pair( single1, single2, rt, self.distance, self.getWorldSize() )
        pair.initialize( self.t )

        return pair

    def createMulti( self ):

        multi = Multi( self )

        #multi.initialize( self.t )

        return multi


    def moveSingle( self, single, pos ):
        single.pos = pos
        self.updateOnParticleMatrix( single.particle )

    def movePair( self, pair, pos ):
        pair.pos = pos
        self.updateOnParticleMatrix( single.particle )


    def addToShellMatrix( self, obj ):
        for i, shell in enumerate( obj.shellList ):
            self.shellMatrix.add( ( obj, i ), shell.pos, shell.radius )

    def removeFromShellMatrix( self, obj ):
        for i in range( len( obj.shellList ) ):
            self.shellMatrix.remove( ( obj, i ) )

    def updateShellMatrix( self, obj ):
        for i, shell in enumerate( obj.shellList ):
            self.shellMatrix.update( ( obj, i ), shell.pos, shell.radius )



    def addEvent( self, t, func, arg ):

        return self.scheduler.addEvent( t, func, arg )

    def addSingleEvent( self, single ):

        eventID = self.addEvent( self.t + single.dt, 
                                 Delegate( self, EGFRDSimulator.fireSingle ), 
                                 single )
        single.eventID = eventID


    def addPairEvent( self, pair ):

        eventID = self.addEvent( self.t + pair.dt, 
                                 Delegate( self, EGFRDSimulator.firePair ), 
                                 pair )
        pair.eventID = eventID

    def addMultiEvent( self, multi ):

        eventID = self.addEvent( self.t + multi.dt, 
                                 Delegate( self, EGFRDSimulator.fireMulti ), 
                                 multi )
        multi.eventID = eventID


    def removeEvent( self, event ):

        self.scheduler.removeEvent( event.eventID )


    def updateEvent( self, t, event ):
        self.scheduler.updateEventTime( event.eventID, t )


    def burstObj( self, obj ):
        
        if isinstance( obj, Single ):
            print 'bursting', obj
            self.burstSingle( obj )
            return [obj,]
        elif isinstance( obj, Pair ):  # Pair
            print 'bursting', obj
            single1, single2 = self.burstPair( obj )
            self.removeEvent( obj )
            self.addSingleEvent( single1 )
            self.addSingleEvent( single2 )
            return [ single1, single2 ]
        else:  # Multi
            print 'bursting', obj
            bursted = self.burstMulti( obj )
            self.removeEvent( obj )
            return bursted

    def burstObjs( self, objs ):

        bursted = []
        for obj in objs:
            b = self.burstObj( obj )
            bursted.extend( b )

        return bursted

    def clearVolume( self, pos, radius, ignore=[] ):

        neighbors = self.getNeighborsWithinRadius( pos, radius, ignore )

        return self.burstObjs( neighbors )


    def burstNonMultis( self, neighbors ):

        bursted = []

        for obj in neighbors:
            if not isinstance( obj, Multi ):
                b = self.burstObj( obj )
                bursted.extend( b )
            else:
                bursted.append( obj )

        return bursted



    def fireSingleReaction( self, single ):

        reactantSpecies = single.particle.species
        oldpos = single.particle.pos.copy()
        
        rt = single.drawReactionType()

        if len( rt.products ) == 0:
            
            self.removeParticle( single.particle )
            
        elif len( rt.products ) == 1:
            
            productSpecies = rt.products[0]

            if reactantSpecies.radius < productSpecies.radius:
                self.clearVolume( oldpos, productSpecies.radius )

            if not self.checkOverlap( oldpos, productSpecies.radius,
                                      ignore = [ single.particle, ] ):
                print 'no space for product particle.'
                raise NoSpace()

            self.removeParticle( single.particle )
            newparticle = self.createParticle( productSpecies, oldpos )
            newsingle = self.createSingle( newparticle )
            self.addToShellMatrix( newsingle )
            self.addSingleEvent( newsingle )
            print 'product;', newsingle

            
        elif len( rt.products ) == 2:
            
            productSpecies1 = rt.products[0]
            productSpecies2 = rt.products[1]
            
            D1 = productSpecies1.D
            D2 = productSpecies2.D
            D12 = D1 + D2
            
            particleRadius1 = productSpecies1.radius
            particleRadius2 = productSpecies2.radius
            particleRadius12 = particleRadius1 + particleRadius2

            # clean up space.
            rad = max( particleRadius12 * ( D1 / D12 ) + particleRadius1,
                       particleRadius12 * ( D2 / D12 ) + particleRadius2 )

            self.clearVolume( oldpos, rad )

            for i in range( 100 ):
                unitVector = randomUnitVector()
                vector = unitVector * particleRadius12 * SAFETY
            
                # place particles according to the ratio D1:D2
                # this way, species with D=0 doesn't move.
                # FIXME: what if D1 == D2 == 0?
                newpos1 = oldpos + vector * ( D1 / D12 )
                newpos2 = oldpos - vector * ( D2 / D12 )

                newpos1 = self.applyBoundary( newpos1 )
                newpos2 = self.applyBoundary( newpos2 )

                # accept the new positions if there is enough space.
                if ( self.checkOverlap( newpos1, particleRadius1,
                                        ignore = [ single.particle, ] ) and
                     self.checkOverlap( newpos2, particleRadius2,
                                        ignore = [ single.particle, ] ) ):
                    break
            else:
                print 'no space for product particles.'
                raise NoSpace()

            self.removeParticle( single.particle )

            particle1 = self.createParticle( productSpecies1, newpos1 )
            particle2 = self.createParticle( productSpecies2, newpos2 )
            newsingle1 = self.createSingle( particle1 )
            newsingle2 = self.createSingle( particle2 )

            self.addToShellMatrix( newsingle1 )
            self.addToShellMatrix( newsingle2 )
            self.addSingleEvent( newsingle1 )
            self.addSingleEvent( newsingle2 )

            print 'products;', newsingle1, newsingle2

        else:
            raise RuntimeError, 'num products >= 3 not supported.'

        self.reactionEvents += 1
        self.populationChanged = True


    def propagateSingle( self, single, r ):
        
        closest, distanceToClosestShell =\
                 self.getClosestObj( single.particle.pos, 
                                     ignore = [ single, ] )

        particleRadius = single.getMinRadius()
        oldpos = single.particle.pos.copy()
        
        displacement = single.drawDisplacement( r )
            
        newpos = oldpos + displacement
        newpos = self.applyBoundary( newpos )
            
        assert self.checkOverlap( newpos, particleRadius,
                                  ignore = [ single.particle, ] )

        self.moveSingle( single, newpos )

        single.initialize( self.t )

        self.updateShellMatrix( single )



    def fireSingle( self, single ):

        # Reaction.
        if single.eventType == EventType.REACTION:

            print 'single reaction', single
            r = single.drawR( single.dt )

            self.propagateSingle( single, r )

            try:
                self.removeFromShellMatrix( single )
                self.fireSingleReaction( single )
            except NoSpace:
                print 'single reaction; placing product failed.'
                self.addToShellMatrix( single )
                self.rejectedMoves += 1
                single.reset()
                return single.dt

            single.dt = -INF  # remove this Single from the Scheduler
            return single.dt

        # Propagate, if not reaction.

        D0 = single.getD()

        # Handle immobile case first.
        if D0 == 0:
            # no propagation, just calculate reaction times.
            single.determineNextEvent() 
            return single.dt
        
        # Propagate this particle to the exit point on the shell.
        
        self.propagateSingle( single, single.getMobilityRadius() )

        # (2) Clear volume.

        minShell = single.getMinRadius() * ( 1.0 + self.SINGLE_SHELL_FACTOR )

        neighbors = self.getNeighborsWithinRadius( single.pos, radius=minShell,
                                                   ignore=[single,] )

        bursted = self.burstNonMultis( neighbors )

        neighborDists = self.objDistanceArray( single.pos, bursted )
        neighbors = [ bursted[i] for i in 
                      ( neighborDists <= minShell ).nonzero()[0] ]

        if neighbors:
            obj, b = self.formPairOrMulti( single, neighbors )
            bursted.extend( b )

            if obj:
                single.dt = -INF # remove by rescheduling to past.
                return single.dt


        # If this Single bursted something,
        # Recheck closest and closest shell distance.

        closest, closestShellDistance =\
            self.getClosestObj( single.pos, ignore = [ single, ] )

        self.updateSingle( single, closest, closestShellDistance )

        bursted = list( set( bursted ) ) # uniq()

        for s in bursted:

            if isinstance( s, Single ):

                c, d = self.getClosestObj( s.pos, ignore = [s,] )
                self.updateSingle( s, c, d )
                self.updateEvent( self.t + s.dt, s )
                #print 'restore shell', s, s.radius, c, d
            



        # (3) If a new Pair was not formed, this Single continues.
        #     Determine a new shell size and dt.

        # recheck the closest and distance to it.


        #print 'single shell', single.radius, 'dt', single.dt

        return single.dt


    def updateSingle( self, single, closest, distanceToShell ): 

        if isinstance( closest, Single ):
            distanceToClosest = self.distance( single.pos, closest.pos )
            shellSize = single.calculateShellSize( closest, distanceToClosest,
                                                   distanceToShell )
        else:  # Pair or Multi
            shellSize = distanceToShell / SAFETY

        shellSize = min( shellSize, self.getMatrixCellSize(),
                         self.maxShellSize )

        single.setRadius( shellSize )
        single.determineNextEvent()
        self.updateShellMatrix( single )



    def firePair( self, pair ):

        print 'fire:', pair, pair.eventType

        particle1 = pair.single1.particle
        particle2 = pair.single2.particle
        species1 = particle1.species
        species2 = particle2.species
        particleRadius1 = species1.radius
        particleRadius2 = species2.radius
        
        oldInterParticle = particle2.pos - particle1.pos

        oldCoM = self.applyBoundary( pair.getCoM() )

        # Three cases:
        #  0. Reaction
        #  1. Escaping through a_r.
        #  2. Escaping through a_R.
        #  3. Single reaction 

        # First handle single reaction case.
        if pair.eventType == 3:

            reactingsingle = pair.reactingsingle

            print 'pair: single reaction', reactingsingle

            if reactingsingle == pair.single1:
                theothersingle = pair.single2
            else:
                theothersingle = pair.single1

            self.burstPair( pair )

            self.addSingleEvent( theothersingle )

            try:
                self.removeFromShellMatrix( reactingsingle )
                self.fireSingleReaction( reactingsingle )
            except NoSpace:
                self.addToShellMatrix( reactingsingle )
                self.rejectedMoves += 1
                reactingsingle.dt = 0
                self.addSingleEvent( reactingsingle )

            pair.dt = -INF
            return pair.dt
        


        #
        # 0. Reaction
        #
        if pair.eventType == EventType.REACTION:

            print 'reaction'

            if len( pair.rt.products ) == 1:
                
                species3 = pair.rt.products[0]

                rnd = numpy.random.uniform( size=2 )

                # calculate new R
            
                r_R = pair.drawR_single( pair.dt, pair.a_R )
            
                displacement_R_S = [ r_R, rnd[0] * Pi, rnd[1] * 2 * Pi ]
                displacement_R = sphericalToCartesian( displacement_R_S )
                newCoM = oldCoM + displacement_R
                
                assert self.distance( oldCoM, newCoM ) + species3.radius <\
                    pair.radius

                #FIXME: SURFACE
                newPos = self.applyBoundary( newCoM )

                self.removeParticle( particle1 )
                self.removeParticle( particle2 )

                particle = self.createParticle( species3, newPos )
                newsingle = self.createSingle( particle )
                self.addToShellMatrix( newsingle )
                self.addSingleEvent( newsingle )

                self.reactionEvents += 1
                self.populationChanged = True
                
            else:
                raise NotImplementedError,\
                      'num products >= 2 not supported.'

            self.removeFromShellMatrix( pair )

            pair.dt = -INF
            return pair.dt


        #
        # Escape 
        #

        r0 = self.distance( particle1.pos, particle2.pos )

        # 1 Escaping through a_r.
        if pair.eventType == EventType.ESCAPE:

            print 'escape r'

            print 'r0 = ', r0, 'dt = ', pair.dt, pair.pgf.dump()
            
            rnd = numpy.random.uniform( size=4 )

            # calculate new R
            
            r_R = pair.drawR_single( pair.dt, pair.a_R )
                
            displacement_R_S = [ r_R, rnd[0] * Pi, rnd[1] * 2 * Pi ]
            displacement_R = sphericalToCartesian( displacement_R_S )
            newCoM = oldCoM + displacement_R

            # calculate new r
            theta_r = pair.drawTheta_pair( rnd[2], pair.a_r, r0, pair.dt )
            phi_r = rnd[3] * 2 * Pi
            newInterParticleS = numpy.array( [ pair.a_r, theta_r, phi_r ] )
            newInterParticle = sphericalToCartesian( newInterParticleS )
                
            newpos1, newpos2 = pair.newPositions( newCoM, newInterParticle,
                                                  oldInterParticle )
            newpos1 = self.applyBoundary( newpos1 )
            newpos2 = self.applyBoundary( newpos2 )


        # 2 escaping through a_R.
        elif pair.eventType == 2:

            print 'escape R'

            rnd = numpy.random.uniform( size = 4 )

            # calculate new r
            print 'r0 = ', r0, 'dt = ', pair.dt, pair.pgf.dump()
            r = pair.drawR_pair( r0, pair.dt, pair.a_r )
            print 'new r = ', r
            #assert r >= pair.sigma
            
            theta_r = pair.drawTheta_pair( rnd[0], r, r0, pair.dt )
            phi_r = rnd[1] * 2*Pi
            newInterParticleS = numpy.array( [ r, theta_r, phi_r ] )
            newInterParticle = sphericalToCartesian( newInterParticleS )
                
            # calculate new R
            displacement_R_S = [ pair.a_R, rnd[2] * Pi, rnd[3] * 2 * Pi ]
            displacement_R = sphericalToCartesian( displacement_R_S )
            
            newCoM = oldCoM + displacement_R
                
            #print 'COM', oldCoM, newCoM

            newpos1, newpos2 = pair.newPositions( newCoM, newInterParticle,
                                                      oldInterParticle )
            newpos1 = self.applyBoundary( newpos1 )
            newpos2 = self.applyBoundary( newpos2 )

        else:
            raise SystemError, 'Bug: invalid eventType.'

        # this has to be done before the following clearVolume()

        self.removeFromShellMatrix( pair )

        assert self.distance( newpos1, newpos2 ) >= pair.sigma

        assert self.checkOverlap( newpos1, particleRadius1,
                                  ignore = [ particle1, particle2 ] )
        assert self.checkOverlap( newpos2, particleRadius2,
                                  ignore = [ particle1, particle2 ] )

        single1, single2 = pair.single1, pair.single2

        self.moveSingle( single1, newpos1 )
        self.moveSingle( single2, newpos2 )

        single1.initialize( self.t )
        single2.initialize( self.t )
            
        self.addSingleEvent( single1 )
        self.addSingleEvent( single2 )

        self.addToShellMatrix( single1 )
        self.addToShellMatrix( single2 )

        pair.dt = -INF
        return pair.dt


    def fireMulti( self, multi ):
        
        self.updateEvent( INF, multi )
        nextObjTime = self.scheduler.getTopTime()

        sim = multi.sim

        startCount = sim.stepCounter

        while nextObjTime >= self.t:
            sim.step()

            if sim.populationChanged:
                print 'bd reaction'
                self.reactionEvents += 1
                
                self.breakUpMulti( multi )

                dt = -INF
                break

            if sim.escaped:
                print 'multi particle escaped.'
                self.breakUpMulti( multi )
                
                dt = -INF
                break
        else:
            dt = multi.dt

        steps = sim.stepCounter - startCount
        assert steps >= 1
        self.stepCounter += steps-1

        print 'multi stepped', steps, 'steps, duration', sim.t

        return dt


    def breakUpMulti( self, multi ):

        self.removeFromShellMatrix( multi )

        singles = []
        for particle in multi.sim.particleList:
            single = self.createSingle( particle )
            self.addToShellMatrix( single )
            self.addSingleEvent( single )
            singles.append( single )

        return singles


    def burstMulti( self, multi ):
        
        singles = self.breakUpMulti( multi )

        return singles


    def burstSingle( self, single ):

        #print 'b', self, 't ', t, 'last ', self.lastTime, 'dt ', self.dt
        assert self.t >= single.lastTime
        assert self.t <= single.lastTime + single.dt
        assert single.radius >= single.getMinRadius()

        dt = self.t - single.lastTime

        particleRadius = single.particle.species.radius
        oldpos = single.particle.pos.copy()

        r = single.drawR( dt )

        displacement = single.drawDisplacement( r )
            
        newpos = oldpos + displacement

        newpos = self.applyBoundary( newpos )

        assert self.distance( newpos, oldpos ) <= single.getMobilityRadius()
        assert self.checkOverlap( newpos, particleRadius,
                                  ignore = [ single.particle, ] )

        single.initialize( self.t )
        self.moveSingle( single, newpos )

        self.updateShellMatrix( single )
        self.updateEvent( self.t, single )


    def breakUpPair( self, pair ):

        assert self.t >= pair.lastTime

        dt = self.t - pair.lastTime 

        if dt != 0.0:

            particle1 = pair.single1.particle
            particle2 = pair.single2.particle

            particleRadius1 = particle1.species.radius
            particleRadius2 = particle2.species.radius
            
            oldpos1 = particle1.pos.copy()
            oldpos2 = particle2.pos.copy()

            oldInterParticle = particle2.pos - particle1.pos
            oldCoM = pair.getCoM()
            r0 = pair.distance( particle1.pos, particle2.pos )
            
            rnd = numpy.random.uniform( size = 4 )

            # calculate new CoM
            r_R = pair.drawR_single( dt, pair.a_R )
            
            displacement_R_S = [ r_R, rnd[0] * Pi, rnd[1] * 2 * Pi ]
            displacement_R = sphericalToCartesian( displacement_R_S )
            newCoM = oldCoM + displacement_R
            
            # calculate new interparticle
            r_r = pair.drawR_pair( r0, dt, pair.a_r )
            theta_r = pair.drawTheta_pair( rnd[2], r_r, r0, dt )
            phi_r = rnd[3] * 2 * Pi
            newInterParticleS = numpy.array( [ r_r, theta_r, phi_r ] )
            newInterParticle = sphericalToCartesian( newInterParticleS )

            newpos1, newpos2 = pair.newPositions( newCoM, newInterParticle,
                                                  oldInterParticle )

        return ( pair.single1, pair.single2 )


    def burstPair( self, pair ):
        print 'burst', pair

        single1, single2 = self.breakUpPair( pair )
        single1.initialize( self.t )
        single2.initialize( self.t )
        
        self.removeFromShellMatrix( pair )
        self.addToShellMatrix( single1 )
        self.addToShellMatrix( single2 )

        return single1, single2


    def formPairOrMulti( self, single, neighbors ):

        assert neighbors

        bursted = []

        # Try forming a Pair.
        if len( neighbors ) == 1 and isinstance( neighbors[0], Single ):
            obj, bursted = self.formPair( single, neighbors[0] )
            if obj:
                return obj, bursted

        # Then, a Multi.
        closest = neighbors[0]


        if isinstance( closest, Single ):

            multi = self.createMulti()
            self.addToMulti( single, multi )
            self.removeFromShellMatrix( single )
            for neighbor in neighbors:
                self.addToMultiRecursive( neighbor, multi )

            multi.initialize( self.t )
            
            self.addToShellMatrix( multi )
            self.addMultiEvent( multi )

            return multi, bursted

        elif isinstance( closest, Multi ):

            multi = closest
            print 'multi merge', single, multi#, closestShellDistance

            self.removeFromShellMatrix( multi )

            self.addToMulti( single, multi )
            self.removeFromShellMatrix( single )
            for neighbor in neighbors[1:]:
                self.addToMultiRecursive( neighbor, multi )

            multi.initialize( self.t )

            self.addToShellMatrix( multi )
            self.updateEvent( self.t + multi.dt, multi )

            return multi, bursted


        assert False, 'do not reach here'


    def formPair( self, single, closest ):

        bursted=[]
        if not closest.isReset():
            self.burstSingle( closest )

        bursted.append( closest )

        obj = self.tryPair( single, closest )
        return obj, bursted



    '''
    Determine if given couple singles meet the criteria of forming
    a Pair.

    '''

    def tryPair( self, single1, single2 ):

        assert single1.isReset()
        assert single2.isReset()

        species1 = single1.particle.species
        species2 = single2.particle.species

        radius1 = species1.radius
        radius2 = species2.radius
        sigma1 = radius1 * 2
        sigma2 = radius2 * 2

        D1, D2 = species1.D, species2.D
        D12 = D1 + D2

        radius12 = radius1 + radius2


        com = calculatePairCoM( single1.pos, single2.pos,\
                                single1.getD(), single2.getD(),\
                                self.getWorldSize() )
        com = self.applyBoundary( com )
        closest, closestShellDistance =\
            self.getClosestObj( com, ignore = ( single1, single2 ) )
        pairDistance = self.distance( single1.pos, single2.pos )

        minShellSize = max( pairDistance * D1 / D12 + radius1,
                            pairDistance * D2 / D12 + radius2 )
        pairShellSizeMargin = min( sigma1, sigma2 ) * self.SINGLE_SHELL_FACTOR

        
        # Squeezing check:
        #     Pair shell size cannot be smaller than minPairShellSize.
        if closestShellDistance <= minShellSize + pairShellSizeMargin:
            print 'Pair not formed because squeezed by ', closest
            return None

        # pairGap = real distance including radii
        pairGap = pairDistance - radius12
        assert pairGap >= 0, 'pairGap between %s and %s = %g < 0' \
            % ( single1, single2, pairGap )


        maxShellSize = min( self.getMatrixCellSize(), self.maxShellSize )


        # 1. Shell cannot be larger than max shell size or sim cell size.
        if minShellSize > maxShellSize:
            print 'tryParing(): minShellSize > maxShellSize: none created'
            return None

        # 2. Check if a Pair is better than two Singles

        # pairDistance < min( particles' distances to the shell )
        if isinstance( closest, Single ):
            shellSize = minShellSize \
                + ( closestShellDistance - minShellSize ) * .5
        else:
            shellSize = closestShellDistance / SAFETY
        
        if shellSize < 1.5 * minShellSize:  #FIXME: factor ok?
            print 'Singles are better than Pair'
            return None

        # 3. Ok, Pair makes sense.  Create one.

        shellSize = max( shellSize, minShellSize )
        shellSize = min( shellSize, maxShellSize )

        pair = self.createPair( single1, single2 )
        pair.setRadius( shellSize )

        self.removeFromShellMatrix( single1 )
        self.removeFromShellMatrix( single2 )
        self.addToShellMatrix( pair )


        pair.determineNextEvent()

        print pair, 'dt=', pair.dt, 'type=', pair.eventType
                
        self.addPairEvent( pair )

        self.removeEvent( single2 )


        print 'Pair formed: ', pair, 'pair distance', pairDistance,\
            'shell size=', pair.radius,\
            ', closest = ', closest,\
            ', distance to shell = ', closestShellDistance


        return pair
    


    def addToMultiRecursive( self, obj, multi ):
        
        if isinstance( obj, Single ):
            if obj.particle in multi.sim.particleList:  # Already in the Multi.
                return
            assert obj.isReset()
            
            self.addToMulti( obj, multi )
            self.removeFromShellMatrix( obj )
            self.removeEvent( obj )

            radius = obj.particle.species.radius *\
                ( 1.0 + self.MULTI_SHELL_FACTOR )
            neighbors = self.getNeighborsWithinRadius( obj.pos, radius,
                                                       ignore=[obj,] )
            bursted = self.burstNonMultis( neighbors )
            neighborDists = self.objDistanceArray( obj.pos, bursted )
            neighbors = [ bursted[i] for i in 
                          ( neighborDists <= radius ).nonzero()[0] ]

            for obj in neighbors:
                self.addToMultiRecursive( obj, multi )

        elif isinstance( obj, Multi ):
            try:
                self.mergeMultis( obj, multi )
                self.removeFromShellMatrix( obj )
                self.removeEvent( obj )
            except AlreadyIn:
                print obj, ' already added. skipping.'
        else:
            assert False, 'do not reach here.'  # Pairs are bursted




    def addToMulti( self, single, multi ):
        print 'adding', single, 'to', multi

        shellSize = single.particle.species.radius * \
            ( 1.0 + self.MULTI_SHELL_FACTOR )
        multi.addParticle( single.particle )
        multi.addShell( single.pos, shellSize )


    '''
        merge multi1 into multi2
    '''
    def mergeMultis( self, multi1, multi2 ):

        print 'merging', multi1, 'to', multi2

        if multi1.sim.particleList[0] in multi2.sim.particleList:
            raise AlreadyIn

        for i, particle in enumerate( multi1.sim.particleList ):
            
            # FIXME: shells should be renewed

            multi2.addParticle( particle )
            shell = multi1.shellList[i]
            multi2.addShell( shell.pos, shell.radius )

        multi2.initialize( self.t )


    '''
    Find closest n shells.

    This method returns a tuple ( neighbors, distances ).
    '''

    def getNeighborShells( self, pos, n, dummy=DummySingle() ):

        return self.shellMatrix.getNeighbors( pos, n, dummy )

    def getNeighborShellsWithinRadius( self, pos, radius ):

        return self.shellMatrix.getNeighborsWithinRadius( pos, radius )


    def getNeighborsWithinRadius( self, pos, radius, ignore=[] ):

        shells, distances = self.getNeighborShellsWithinRadius( pos, radius )

        neighbors = [ s[0] for s in shells if s[0] not in ignore ]
        neighbors = list( set( neighbors ) )  # uniq()

        return neighbors

    '''
    Find the closest shell from the position pos.

    When the sequence parameter ignore is given, this method finds the
    closest ignoring the objects in it.

    This method returns a tuple ( neighbors, distances ).
    '''

    def getClosestShell( self, pos, ignore=[] ):

        neighbors, distances = self.getNeighborShells( pos, len( ignore ) + 1 )

        for i, neighbor in enumerate( neighbors ):
            if neighbor not in ignore:
                return neighbor, distances[i]

        # default case: none left.
        return ( DummySingle(), 0 ), INF

    def getClosestObj( self, pos, ignore=[] ):

        neighbors, distances = self.getNeighborShells( pos, len( ignore ) + 1 )

        for i, neighbor in enumerate( neighbors ):
            neighborObj = neighbor[0]
            if neighborObj not in ignore:
                return neighborObj, distances[i]

        # default case: none left.
        return DummySingle(), INF


    def objDistance( self, pos, obj ):
        
        dists = numpy.zeros( len( obj.shellList ) )
        for i, shell in enumerate( obj.shellList ):
            dists[i] = self.distance( pos, shell.pos ) - shell.radius

        return min( dists )

    def objDistanceArray( self, pos, objs ):

        dists = numpy.array( [ self.objDistance( pos, obj ) for obj in objs ] )
        return dists
            


    '''
    Get neighbors simply by distance.

    This method does not take into account of particle radius or shell size.
    This method picks top n neighbors simply by the distance between given
    position pos to positions of objects (either Singles or Pairs) around.

    This method returns a tuple ( neighbors, distances ).

    def getNeighbors( self, pos, n=None ):
        return self.shellMatrix.getNeighbors( pos, n )

    def getNeighbors( self, pos, n=None ):

        objMatrix = self.matrix

        size = objMatrix.size
        if not n:
            n = size 

        distances = self.distanceSqArray( objMatrix.positions, pos )

        topargs = distances.argsort()[:n]
        distances = distances.take( topargs )
        distances = numpy.sqrt( distances )
        neighbors = [ objMatrix.objList[arg] for arg in topargs ]

        return neighbors, distances

    def getClosestNeighbor( self, pos, ignore=[] ):

        neighbors, distances = self.getNeighbors( pos, len( ignore ) + 1 )

        for i in range( len( neighbors ) ): 
            if neighbors[i] not in ignore:
                closest, distance = neighbors[i], distances[i]

                #assert not closest in ignore
                return closest, distance

        # default case: none left.
        return DummySingle(), INF
    '''


    #
    # consistency checkers
    #

    
    def checkObj( self, obj ):

        obj.check()


        allshells = [ ( obj, i ) for i in range( len( obj.shellList ) ) ]
        for i, shell in enumerate( obj.shellList ):
            shellid = (obj,i)
            closest, distance = self.getClosestShell( shell.pos,
                                                      ignore = allshells )
            radius = shell.radius

            if radius > self.getMatrixCellSize():
                raise RuntimeError,\
                    '%s shell size larger than simulator cell size' % \
                    str( shellid )

            if radius > self.maxShellSize:
                raise RuntimeError,\
                    '%s shell size larger than maxShellSize' % \
                    str( shellid )

            if distance - radius < 0.0:
                raise RuntimeError,\
                    '%s overlaps with %s. (shell: %g, dist: %g, diff: %g.' \
                    % ( str( obj ), str( closest[0] ), radius, distance,\
                            distance - radius )


    def checkObjForAll( self ):

        for i in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex(i).getArg()
            self.checkObj( obj )


    def checkEventStoichiometry( self ):

        population = 0
        for species in self.speciesList.values():
            population += species.pool.size

        eventPopulation = 0
        for i in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex(i).getArg()
            eventPopulation += obj.multiplicity

        if population != eventPopulation:
            raise RuntimeError, 'population %d != eventPopulation %d' %\
                  ( population, eventPopulation )

    def checkShellMatrix( self ):

        if self.worldSize != self.shellMatrix.worldSize:
            raise RuntimeError,\
                'self.worldSize != self.shellMatrix.worldSize'

        shellPopulation = 0
        for i in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex(i).getArg()
            shellPopulation += len( obj.shellList )

        if shellPopulation != self.shellMatrix.size:
            raise RuntimeError,\
                'num shells != self.shellMatrix.size'
        
        self.shellMatrix.check()

        for key in self.shellMatrix.keyList:
            pos, radius = self.shellMatrix.get( key )
            obj, i = key

            if ( obj.shellList[i].pos - pos ).any():
                raise RuntimeError, \
                    '%s shellMatrix positions consistency broken' % str( key )

            if obj.shellList[i].radius != radius:
                raise RuntimeError, \
                    '%s shellMatrix radii consistency broken' % str( key )



    def check( self ):

        GFRDSimulatorBase.check( self )

        assert self.t >= 0.0
        assert self.dt >= 0.0

        self.checkShellMatrix()

        self.checkEventStoichiometry()
        
        self.checkObjForAll()




    #
    # methods for debugging.
    #

    def dumpScheduler( self ):
        scheduler = self.scheduler
        for i in range( scheduler.getSize() ):
            event = scheduler.getEventByIndex(i)
            print i, event.getTime(), event.getArg()

    def dump( self ):
        scheduler = self.scheduler
        for i in range( scheduler.getSize() ):
            event = scheduler.getEventByIndex(i)
            print i, event.getTime(), event.getArg(), event.getArg().pos



