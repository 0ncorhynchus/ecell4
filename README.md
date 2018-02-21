E-Cell System version 4
=======================

[![Build Status](https://travis-ci.org/ecell/ecell4.svg?branch=develop)](https://travis-ci.org/ecell/ecell4)
[![CircleCI](https://circleci.com/gh/ecell/ecell4.svg?style=svg)](https://circleci.com/gh/ecell/ecell4)
[![Build status](https://ci.appveyor.com/api/projects/status/github/ecell/ecell4?svg=true)](https://ci.appveyor.com/project/kaizu/ecell4)
[![Documentation Status](https://readthedocs.org/projects/ecell4/badge/?version=latest)](http://ecell4.readthedocs.org/en/latest/?badge=latest)
[![PyPI](https://img.shields.io/pypi/v/ecell.svg)](https://pypi.python.org/pypi/ecell)
[![License: GPL v2](https://img.shields.io/badge/license-GPL%20v2-blue.svg)](https://github.com/ecell/ecell4/blob/master/licenses/LICENSE)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/87e076986e354b508f66af0a0ca3373d)](https://www.codacy.com/app/ecell/ecell4?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ecell/ecell4&amp;utm_campaign=Badge_Grade)
[![Slack Status](https://img.shields.io/badge/chat-on%20slack-50baa6.svg)](https://ecell-project.herokuapp.com/)
<!---[![Slack Status](https://ecell-project.herokuapp.com/badge.svg)](https://ecell-project.herokuapp.com/)--->

What is E-Cell System?
----------------------

E-Cell System is, a software platform for modeling, simulation and analysis of complex, heterogeneous and multi-scale systems like the cell.
E-Cell has multi-algorithm, multi-timescale and multi-spatial-representation as its central feature.

Features
--------

- Single particle simulations ([The enhanced Green's Function Reaction Dynamics (eGFRD) method](http://gfrd.org), [Spatiocyte](http://spatiocyte.org), and the Reaction Brownian Dynamics (RBD) method)
- Multi-algorithm support
- Rule-based modeling
- Python programmable

Try online
----------

You can try this package online from the following links:

<a href="https://notebooks.azure.com/import/gh/ecell/ecell4"><img src="https://notebooks.azure.com/launch.png" /></a>
[![Binder](http://mybinder.org/badge.svg)](http://beta.mybinder.org/v2/gh/ecell/ecell4-notebooks/master)

- Code fragments that depend on ffmpeg will not work with Azure Notebooks. If you use movie export, please try binder instead.
- If you use Azure Notebooks, please move `readthedocs/tutorials` or `readthedocs/examples`. There are several Jupyter Notebooks.

Installation
-------------

Please see [INSTALL.md](https://github.com/ecell/ecell4/blob/master/INSTALL.md).

How to try E-Cell4 examples
---------------------------
Here we download example notebooks from https://github.com/ecell/ecell4-notebooks and open it with Jupyter Notebook.

### Windows
Open powershell and run these commands.
Here we assume that you installed Miniconda(Python3.5) to C:¥Miniconda3

```shell
cd C:¥Miniconda3¥Scripts
wget https://github.com/ecell/ecell4-notebooks/archive/master.zip -OutFile master.zip
Expand-Archive master.zip
.¥jupyter-notebook.exe .¥master¥ecell4-notebooks-master¥
```

### Mac or Linux
Here we assume that you installed Miniconda(Python3.5) to ~/miniconda3

```shell
wget https://github.com/ecell/ecell4-notebooks/archive/master.zip
unzip master.zip
cd ecell4-notebooks-master
~/miniconda3/bin/jupyter-notebook
```

Simple examples
---------------

Here are two extremely simple examples, See http://ecell4.readthedocs.org for more details on running E-Cell4.

```
Python 2.7.6 (default, Mar 22 2014, 22:59:56)
[GCC 4.8.2] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> from ecell4.core import *
>>> sp = Species("A.B.C")
>>> print sp.serial()
A.B.C
>>>
```

### Binding and unbinding reactions

```python
%matplotlib inline
from ecell4 import *

with reaction_rules():
    A + B == C | (0.01, 0.3)

run_simulation(100, {'A': 60, 'B': 60})
```

![png](./readthedocs/images/output_7_0.png)

### Diffusion on a spherical surface

```python
%matplotlib inline
from ecell4 import *

with species_attributes():
    A | {'D': '1', 'location': 'M'}

surface = Sphere(ones() * 0.5, 0.5).surface()
obs = FixedIntervalTrajectoryObserver(1e-4)
run_simulation(
    0.4, y0={'A': 10}, structures={'M': surface},
    solver='spatiocyte', observers=obs, return_type=None)

viz.plot_trajectory(obs, interactive=False)
```

![png](./readthedocs/images/hairball.png)

Citation
========

If this package contributes to a project which leads to a scientific publication, I would appreciate a citation.

[![DOI](https://zenodo.org/badge/6348303.svg)](https://zenodo.org/badge/latestdoi/6348303)

Licensing terms
===============

This product is licensed under the terms of the [GNU General Public License v2](https://github.com/ecell/ecell4/blob/master/licenses/LICENSE),
See [NOTICE](https://github.com/ecell/ecell4/blob/master/licenses/NOTICE.txt) for the software included in this product.

- Copyright (c) 2010-, RIKEN

All rights reserved.
