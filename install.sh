#!/bin/bash -x

# PREFIX=/usr/local
# PREFIX=${HOME}/local
# PREFIX=

make clean; rm -rf ${PREFIX}; rm -rf python/build; rm CMakeCache.txt
cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} .
# cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} -DWITH_LATTICE=OFF .
make
make test
make install
cd python
python -c 'import sys; print sys.path'
PYTHONPATH=${PREFIX}/lib/python2.7/site-packages${PYTHONPATH} python -c 'import sys; print sys.path'
LD_LIBRARY_PATH=${PREFIX}/lib python setup.py build_ext -L${PREFIX}/lib -I${PREFIX}/include install --prefix=${PREFIX}
PYTHONPATH=${PREFIX}/lib/python2.7/site-packages${PYTHONPATH} LD_LIBRARY_PATH=${PREFIX}/lib python setup.py test
