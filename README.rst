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
- Cython_

GTK
sudo apt-get install libgtk-3-*

para compilar
cc source_file.c $(pkg-config --cflags --libs gtk+-3.0)

Ou:
baixe a boost no site
http://www.boost.org/
Descompacte o arquivo
Dê o comando:
$ sudo ./bootstrap.sh
$ sudo ./bjam install
#Baixe a Eigen em
http://eigen.tuxfamily.org/
Execute o cmake como
Entre na pasta descompactada do eigen
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr/
sudo make install

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
.. _Eigen3: http://eigen.tuxfamily.org/
.. _Boost:  http://www.boost.org/
.. _Cython: http://cython.org/
.. _CMake:  http://cmake.org/
.. _Nose:   http://pypi.python.org/pypi/nose/
.. _Python: http://python.org/

PYTHON
http://docs.cython.org/
http://docs.python.org/
http://ipython.org/
http://docs.enthought.com/EPD_8/
http://docs.enthought.com/mayavi
http://sphinx-doc.org/
http://docs.scipy.org/
http://www.reportlab.com/

GIT
Para remover arquivos deletados da pasta:
git rm $(git ls-files --deleted)
adicionar uma pasta
git add pasta
atualizar o repositório no site
git push origin master
