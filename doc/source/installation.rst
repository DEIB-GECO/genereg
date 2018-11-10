Installation
============================================

-------------
Prerequisites
-------------

Python (version >= 3.6) is required.

The library is currently tested for running on Windows. Portability on other OS (iOS and Linux) is being tested.

**PyGMQL** package is required. Follow the instructions at this `link <https://github.com/DEIB-GECO/PyGMQL>`_
to install it and to known more about it.
PyGMQL can be used in two different ways: either locally as any other computational library, or through the connection to a remote server. In order to use the execution through a remote server, loggin in to the remote service is required. You need to configure PyGMQL by specifying the remote server address (the Web service offered by the GeCo group at Politecnico di Milano can be found `here <http://www.gmql.eu/gmql-rest/>`_ and then loggin in into the system.
Here, the links to follow for the installation and configuration of the library: `installation <https://pygmql.readthedocs.io/en/latest/installation.html>`_ and `remote data management <https://pygmql.readthedocs.io/en/latest/remote.html>`_.


----------------------
Installation using PIP
----------------------
The package (available on `PyPI <https://pypi.org/project/genereg/>`_) can be downloaded and installed directly in your Python distribution using::

    pip install genereg
