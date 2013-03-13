ff
====================

Introdução
------------

Dependências
------------
Build Dependencies
^^^^^^^^^^^^^^^^^^

- Python_
- CMake_
- C++ compiler (g++ for instance)
- Eigen3_
- Boost_ 
sudo apt-get install libboost-all-dev

- Cython_


GTK:
----------------
- Instalar GTK e informações de como compilar um programa.
- Small Change
- Another small change
::

	sudo apt-get install libgtk-3-*
	cc source_file.c $(pkg-config --cflags --libs gtk+-3.0)

..

Instalar Eigen
------------------
- Descompacte o arquivo do site da Eigen

::

  mkdir build
  cd build
  cmake .. -DCMAKE_INSTALL_PREFIX=/usr/
  make
  sudo make install
  
..

Test Dependencies
^^^^^^^^^^^^^^^^^

-Nose_

Build Instructions
------------------

::

  mkdir build
  cd build
  cmake ..
  make
  
.. warning::

  bla bla

To run the tests::

  nosetests


Documentação
^^^^^^^^^^^^^^^^^

PYTHON
------------------
 - CythonDoc_ 
 - PythonDoc_
 - IPython_
 - EPD8_
 - Mayavi_
 - Sphinx_
 - Scipy_
 - Reportlab_
 - IntroNumericPython_

Artigos matemáticos
------------------
 - PoissonEquationWithPython_
 - Fipy_
 - NavierStokesPython_
 - kmkns_
GIT
-------------------
 - Para remover arquivos deletados da pasta:
 - git rm $(git ls-files --deleted)
 - adicionar uma pasta
 - git add pasta
 - criar commit
 - git commit -m "comentário"
 - atualizar o repositório no site
 - git push origin master

.. _Eigen3: http://eigen.tuxfamily.org/
.. _Boost:  http://www.boost.org/
.. _Cython: http://cython.org/
.. _CMake:  http://cmake.org/
.. _Nose:   http://pypi.python.org/pypi/nose/
.. _Python: http://python.org/
.. _CythonDoc: http://docs.cython.org/
.. _PythonDoc: http://docs.python.org/
.. _IPython: http://ipython.org/
.. _EPD8: http://docs.enthought.com/EPD_8/
.. _Sphinx: http://sphinx-doc.org/
.. _Mayavi: http://docs.enthought.com/mayavi
.. _Scipy: http://docs.scipy.org/
.. _Reportlab: http://www.reportlab.com/
.. _PoissonEquationWithPython: http://www.scientificpython.net/1/post/2012/05/poisson-equation-on-the-square.html
.. _Fipy: http://www.hasenkopf2000.net/wiki/page/fipy-solving-pdes-python/
.. _NavierStokesPython: http://fenicsproject.org/documentation/dolfin/1.0.0/python/demo/pde/navier-stokes/python/documentation.html
.. _kmkns: http://code.google.com/p/kmkns/
.. _IntroNumericPython: http://math.jacobs-university.de/oliver/teaching/numpy-intro/numpy-intro/index.html
