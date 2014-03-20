    PyRosetta installation.


1. Download appropriate PyRosetta <your-architecture>.tar.bz2 file.
2. Unpack archive using 'tar -vjxf PyRosetta-<version>.tar.bz2' command. Please note, there is no
special install procedure required, after unpacking PyRosetta is ready to use. So unpack it to the
location from where you want to execute it.


    Using PyRosetta.

To use PyRosetta you have two options:
  1. Run it from it file system location.
or
  2. Set up environment variables so you can use it from any file system location.

To set up environment variables we provide shell script. To use it execute
'source <path_to_your_PyRosetta_install>/SetPyRosettaEnvironment.sh' each time you planning to use
PyRosetta. This line could also be added to your .bashrc or global shell configuration file.


    iPython.
The PyRosetta distribution contains the iPython package.  iPython is an advanced shell for
Python available under a BSD-like license.  For more information see http://ipython.scipy.org.
iPython provides additional features like tab completion and instant help, so we *strongly*
encourage you to use iPython when working with PyRosetta.  We also provide an iPyRosetta shell
script that sets the PyRosetta environment variables and then executes ipython. Simply call the script
'<path_to_PyRosetta_install>/iPyRosetta [my_script.py]'
