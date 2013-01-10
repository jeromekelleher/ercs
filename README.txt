
==============================================
Extinction/recolonisation coalescent simulator
==============================================

Simulates the coalescent for populations evolving in a spatial 
continuum under the extinction/recolonisation model. The simulations 
support:
        
- A sample of ``n`` individuals with ``m`` loci at arbitrary locations on a 
  torus of diameter ``L``.
- Arbitrary recombination rates between adjacent loci.
- An arbitrary number of classes of event occuring at fixed
  rates. 

Ercs supports both Python 2 and 3.


-------------
Documentation
-------------

Here's a quick example for the impatient::

        import ercs
        sim = ercs.Simulator(10)
        sim.sample = [None, (3, 2), (6, 4), (7, 0)]
        sim.event_classes = [ercs.DiscEventClass(u=0.5, r=1)]
        pi, tau = sim.run(1)

Full documentation for ``ercs`` is available at `<http://jeromekelleher.github.com/ercs>`_.

------------
Installation
------------

Ercs depends on the `GNU Scientific Library <http://www.gnu.org/software/gsl/>`_,
which must be installed before the ``ercs`` module can be built.
Fortunately, this is straightforward on most platforms. For example, 
on Debian or Ubuntu use::

        $ sudo apt-get install libgsl0-dev

or on Fedora::

        $ sudo yum install gsl-devel

GSL is available on most packaging systems; if it is not available on your
platform, it can be installed from source.

Once GSL has been installed we can build the ``ercs`` module using the 
standard Python `methods <http://docs.python.org/install/index.html>`_. For 
example, using pip we have ::
        
        $ sudo pip install ercs

Or, we can manually download the package, unpack it and then run::
        
        $ python setup.py build
        $ sudo python setup.py install

Most of the time this will compile and install the module without difficulty.

It is also possible to download the latest development version of 
``ercs`` from `github <https://github.com/jeromekelleher/ercs>`_. 

******************
Potential problems
******************

On platforms that GSL is not available as part of the native packaging 
system (or GSL was installed locally because of non-root access)
there can be issues with finding the correct headers and libraries
when compiling ``ercs``. For example, on FreeBSD we get something 
like this::

        $ python setup.py build
        ... [Messages cut for brevity] ...
        In file included from _ercsmodule.c:27:
        lib/ercs.h:26:25: error: gsl/gsl_rng.h: No such file or directory
        In file included from _ercsmodule.c:27:
        lib/ercs.h:73: error: expected declaration specifiers or '...' before 'gsl_rng'
        lib/ercs.h:94: error: expected specifier-qualifier-list before 'gsl_rng'
        _ercsmodule.c: In function 'pyercs_simulate':
        _ercsmodule.c:351: error: 'ercs_t' has no member named 'pi'
        _ercsmodule.c:356: error: 'ercs_t' has no member named 'tau'
        error: command 'cc' failed with exit status 1

This can be remedied by using the ``gsl-config`` program to set the 
the ``LDFLAGS`` and ``CFLAGS`` environment variables to 
their correct values::
        
         $ CFLAGS=`gsl-config --cflags` LDFLAGS=`gsl-config --libs` python setup.py build

*****
Tests
*****

Ercs provides some test cases to ensure that the installation has gone smoothly.
It is a good idea to run these immediately after installation::

        $ python tests.py

It is also possible to run some tests on the low-level C library.
To do this, ``cd`` to the ``lib`` directory and run::

        $ make check 

This will build the C library test cases and run them. If compilation fails, it 
may be necessary set ``LDFLAGS`` and ``CFLAGS`` as before:: 

        $ make CFLAGS="`gsl-config --cflags`" LDFLAGS="`gsl-config --libs`" check 


****************
Tested platforms
****************

Ercs has been successfully built and tested on the following platforms:

================        ========        ======          ========
Operating system        Platform        Python          Compiler
================        ========        ======          ========
Ubuntu 8.04             i386            2.5.2           gcc 4.2.3 
NetBSD 5.0              i386            2.7.3           gcc 4.1.3
Fedora 17               i386            2.7.3           gcc 4.7.2
Fedora 17               i386            3.2.3           gcc 4.7.2
Cygwin                  i386            2.6.8           gcc 4.5.3
Ubuntu 12.04            x86-64          2.7.3           gcc 4.6.3
Ubuntu 12.04            x86-64          3.2.3           gcc 4.6.3
FreeBSD 9.0             i386            3.2.2           gcc 4.2.2        
FreeBSD 9.0             i386            2.7.2           gcc 4.2.2        
FreeBSD 9.0             i386            3.1.4           clang 3.0 
Solaris 11              x86-64          2.6.4           Sun C 5.12
Mac OSX 10.6.8          x86-64          2.6.1           gcc 4.2.1
Mac OSX 10.6.8          x86-64          3.2.3           gcc 4.2.1
Mac OS X 10.4.11        ppc             3.2.3           gcc 4.0.1
Mac OS X 10.4.11        ppc             2.7.3           gcc 4.0.1
Debian wheezy           armv6l          2.7.3           gcc 4.6.3
Debian squeeze          ppc64           2.6.6           gcc 4.4.5	
Debian squeeze          ppc64           3.1.3           gcc 4.4.5	
================        ========        ======          ========

The C library has additionally been compiled and tested with the 
following processors and compilers:

==========================================        ========
Platform                                          Compiler
==========================================        ========
Linux 2.6.18-164.15.1.el5 x86_64 GNU/Linux        Intel(R) C 11.0
Linux 2.6.32-5-amd64 x86_64 GNU/Linux             gcc 4.4.5
Linux 2.6.32-5-amd64 x86_64 GNU/Linux             clang 1.1 
Linux 3.2.0-32-generic x86_64 GNU/Linux           gcc 4.6.3
Linux 3.2.0-32-generic x86_64 GNU/Linux           clang 3.0 
Linux 3.2.27+ armv6l GNU/Linux                    gcc 4.6.3
Linux 2.6.32-5-powerpc64 ppc64 GNU/Linux          gcc 4.4.5
Linux 2.6.32-5-powerpc64 ppc64 GNU/Linux          clang 1.1 
SunOS 5.11 11.0 i86pc i386 i86pc                  Sun C 5.12
SunOS 5.10 sun4u sparc SUNW,Ultra-4               Sun C 5.8
==========================================        ========

