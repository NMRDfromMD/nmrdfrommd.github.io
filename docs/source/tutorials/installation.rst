.. include:: ../additional/links.rst

Installation
============

..
    Install the latest published version
    ------------------------------------

    To install the latest stable release, run the following command in a terminal:

    .. code-block:: bash

        pip install nmrdfrommd

Install the development version
-------------------------------

To install the latest development version of NMRDfromMD, clone the repository,
|NMRDfromMD-code|, from GitHub, and use ``pip`` from the main directory:

.. code-block:: bash

    git clone https://github.com/NMRDfromMD/nmrdfrommd.git

    cd nmrdfrommd/

    pip install .

Run the tests
-------------

To run the test suite, install the testing dependencies:

.. code-block:: bash

    pip install pytest coverage

Then run from the tests folder:

.. code-block:: bash

    pytest
