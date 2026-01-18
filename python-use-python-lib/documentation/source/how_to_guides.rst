Installing the application
**************************

There are multiple ways of installing the application libs. The recommended way is to use `poetry` package manager, which will create virtual environment for the project and install all required dependencies.

Alternative way is to install the application libs using `pip` package manager within already existing virtual environment.

Installation the application libs with `poetry` from provided source code
=========================================================================


1. Installing  and initializing `poetry` through `pipx`

Firstly, `pipx` should be installed.

.. code-block::

        > python -m pip install --user pipx 
        > python -m pipx ensurepath


When `pipx` is installed, the installation of `poetry` can be executed,

.. code-block::

        > pipx install poetry

2. Check if `poetry` is successfully installed

.. code-block::

        > poetry --version

3. Install project's packets and documentation builder packets with `poetry` 

.. code-block::

        > poetry env use 3.13

    The previous command will create virtual environment based on `python3.13` in subdirectory `/.venv`. After that, `poetry` will installed specified libs, as well as libs that are used for documentation.

.. code-block::

        > poetry install --with docs


4. Install `universal-optimizer-lib` in editable form  

Library `universal-optimizer-lib`, described `here <https://matf-r-i.github.io/universal-optimizer-lib-python/>`_ should be installed. Source code of that library id within folder `lib`. 

.. code-block::

        > pip install -e libs


Installation the application libs with `pip` from provided source code
======================================================================

1. Create and activate virtual environment

.. code-block::

        > python -m venv venv
        > source venv/bin/activate      # On Linux
        > .\venv\Scripts\activate       # On Windows

2. Install required dependencies from `requirements.txt` file

.. code-block::

        > pip install -r requirements.txt

3. Install `universal-optimizer-lib` in editable form  

Library `universal-optimizer-lib`, described `here <https://matf-r-i.github.io/universal-optimizer-lib-python/>`_ should be installed. Source code of that library id within folder `lib`. 

.. code-block::

        > pip install -e libs


Running of all the unit tests within application
************************************************

Unit tests are located within `tests` subdirectory . To run all the unit tests, the following commands can be used.

.. code-block::

        > python -m unittest

Obtaining coverage analysis of all unit tests within library can be done with `coverage` package. The following commands can be used:


.. code-block::

        > python -m coverage run -m unittest
        > python -m coverage report


Building documentation for the application
*******************************************

1. Build documentation sources into `/documentation/source` folder from `python` source files 

.. code-block::

        > sphinx-apidoc -o documentation/source/ opt


2. Change current directory to `/documentation` 

.. code-block::

        > cd documentation

3. Clean previously builded HTML documentation 

.. code-block::

        /documentation> ./make clean html

4. Build HTML documentation from `/documentation/source` directory. Created documentation is within `/documentation/build/html` directory. 

.. code-block::

        /documentation> ./make html

5. Generated documentation, that is in folder `/documentation/build/html` should be then copied into folder `/docs`.

.. code-block::

        /documentation> cp build/html/*.* ../docs


