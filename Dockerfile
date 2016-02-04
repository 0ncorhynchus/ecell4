FROM ubuntu:14.04

# E-Cell4 wheel
RUN apt-get update; apt-get install -y python python-dev cmake gcc g++ libboost-dev libgsl0-dev libhdf5-dev wget; wget https://bootstrap.pypa.io/get-pip.py; python get-pip.py; pip install cython
ADD . /usr/src/ecell4
RUN cd /usr/src/ecell4; cmake .; make BesselTables; cd python; python setup.py build_ext; python setup.py bdist_wheel; ls dist

# matplotlib and jupyter
RUN apt-get install -y libfreetype6-dev libpng-dev pkg-config python-numpy pandoc
RUN pip install matplotlib jupyter

# ffmpeg and avconv
RUN apt-get install -y software-properties-common libav-tools
RUN add-apt-repository ppa:mc3man/trusty-media -y; apt-get update; apt-get install -y ffmpeg

# install with install.sh
RUN cd /usr/src/ecell4; export PREFIX=/usr/local; export PYTHONPATH=/usr/local/lib/python2.7/site-packages:$PYTHONPATH; ./install.sh --python2 --hdf5

EXPOSE 8888
CMD LD_LIBRARY_PATH=/usr/local/lib jupyter-notebook --notebook-dir='/usr/src/ecell4/ipynb' --no-browser --ip='*' --port 8888
