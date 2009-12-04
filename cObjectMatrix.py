#!/usr/bin/env python


import numpy

import _gfrd

from utils import *

__all__ = (
    'ObjectMatrix',
    )


class ObjectMatrix( object ):
    def __init__( self ):
        self.worldSize = 1.0
        self.setMatrixSize( 3 )
        self.initialize()

    def initialize( self ):
        self.impl = _gfrd.ObjectContainer( self.worldSize, 
                                           self.matrixSize )

    def __iter__(self):
        return self.impl.__iter__()

    def setWorldSize( self, size ):
        self.worldSize = size
        self.initialize()

    def setMatrixSize( self, size ):
        if size < 3:
            raise RuntimeError,\
                'Size of distance cell matrix must be at least 3'
        self.matrixSize = size
        self.initialize()

    def getSize( self ):
        return self.impl.size()
    size = property( getSize )

    def getCellSize( self ):
        return self.impl.cell_size
    cellSize = property( getCellSize )

    def clear( self ):
        self.initialize()

    def remove( self, key ):
        self.impl.erase( key )

    def update( self, key, pos, radius ):
        assert radius < self.cellSize * .5
        self.impl.update( key, (pos, radius) )

    def get( self, key ):
        return self.impl.get( key )

    def getNeighborsCyclicNoSort( self, pos ):
        return self.impl.get_neighbors_cyclic( pos )

    def getNeighborsWithinRadiusNoSort( self, pos, radius ):
        return self.impl.get_neighbors_within_radius( pos, radius )

    def getNeighborsCyclic( self, pos, n=None ):
        neighbors, distances = self.impl.get_neighbors_cyclic( pos )
        topargs = distances.argsort()[:n]
        return neighbors.take( topargs ), distances.take( topargs )

    def getNeighborsWithinRadius( self, pos, radius ):
        neighbors, distances = \
            self.impl.get_neighbors_within_radius( pos, radius )
        topargs = distances.argsort()
        return neighbors.take( topargs ), distances.take( topargs )

    def getNeighbors( self, pos, n=None ):
        return self.getNeighborsCyclic( pos, n )

    def getNeighborsNoSort( self, pos ):
        return self.getNeighborsCyclicNoSort( pos )


    def check( self ):

        pass

