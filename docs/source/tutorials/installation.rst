.. include:: ../additional/links.rst

Installation
============

Install the latest published version
------------------------------------

To install the latest stable release, run the following command in a terminal:

.. code-block:: bash

    pip install nmrdfrommd

Install the development version
-------------------------------

To install the latest development version of NMRDfromMD, clone the repository
from |NMRDfromMD-code| and use ``pip`` from the main directory:

.. code-block:: bash

    git clone https://github.com/NMRDfromMD/nmrdfrommd.git

    cd nmrdfrommd/

    pip install .

You can then run the test suite using ``pytest``:

.. code-block:: bash

    cd tests
    pytest .
