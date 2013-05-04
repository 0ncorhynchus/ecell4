

from cython.operator cimport dereference as deref

from libcpp.vector cimport vector
from libcpp cimport set
from libcpp cimport bool

include "types.pxi"

from PyEcell4 cimport *
from PyEcell4 import PySpecies


cdef extern from "ecell4/gillespie/GillespieWorld.hpp" namespace "ecell4::gillespie":
    cdef cppclass GillespieWorld:
        GillespieWorld(Real) except +
        void set_t(Real)
        Real t()
        Real volume()
        Integer num_species()
        bool has_species(Species &)
        Integer num_molecules(Species &)
        
        void add_species(Species &)
        void remove_species(Species &)
        void add_molecules(Species &sp, Integer &num)
        void remove_molecules(Species &sp, Integer &num)

cdef extern from "boost/shared_ptr.hpp" namespace "boost":
    cdef cppclass shared_ptr[T]:
        shared_ptr(T *ptr)
        T* get()

cdef class PyGillespieWorld:
    #cdef GillespieWorld *thisptr
    cdef shared_ptr[GillespieWorld] *thisptr
    # XXX
    # If you don't use shared_ptr, please remove calling 'get()'.
    def __cinit__(self, Real vol):
        #self.thisptr = new GillespieWorld(vol)
        self.thisptr = new shared_ptr[GillespieWorld](new GillespieWorld(vol))
    def __dealloc__(self):
        #XXX Here, we release shared pointer, and if reference count to the GillespieWorld object,
        # it will be released automatically.
        del self.thisptr
    
    def set_t(self, Real t):
        self.thisptr.get().set_t(t)
    def t(self):
        return self.thisptr.get().t()
    def volume(self):
        return self.thisptr.get().volume()
    def num_species(self):
        return self.thisptr.get().num_species()
    def has_species(self, PySpecies sp):
        return self.thisptr.get().has_species( deref(sp.thisptr) )
    def num_molecules(self, PySpecies sp):
        return self.thisptr.get().num_molecules( deref(sp.thisptr) )

    def add_species(self, PySpecies sp):
        self.thisptr.get().add_species(deref(sp.thisptr) )
    def remove_species(self, PySpecies sp):
        self.thisptr.get().remove_species(deref(sp.thisptr))
    def add_molecules(self, PySpecies sp, Integer num):
        self.thisptr.get().add_molecules(deref(sp.thisptr), num)
    def remove_species(self, PySpecies sp, Integer num):
        self.thisptr.get().remove_molecules(deref(sp.thisptr), num)
