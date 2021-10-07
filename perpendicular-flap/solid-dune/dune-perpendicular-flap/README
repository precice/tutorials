Preparing the Sources
=========================

Additional to the software mentioned in README you'll need the
following programs installed on your system:

  cmake >= 3.13

Getting started
---------------

If these preliminaries are met, you should run

  dunecontrol all

which will find all installed dune modules as well as all dune modules
(not installed) which sources reside in a subdirectory of the current
directory. Note that if dune is not installed properly you will either
have to add the directory where the dunecontrol script resides (probably
./dune-common/bin) to your path or specify the relative path of the script.

Most probably you'll have to provide additional information to dunecontrol
(e. g. compilers, configure options) and/or make options.

The most convenient way is to use options files in this case. The files
define four variables:

CMAKE_FLAGS      flags passed to cmake (during configure)

An example options file might look like this:

#use this options to configure and make if no other options are given
CMAKE_FLAGS=" \
-DCMAKE_CXX_COMPILER=g++-5 \
-DCMAKE_CXX_FLAGS='-Wall -pedantic' \
-DCMAKE_INSTALL_PREFIX=/install/path" #Force g++-5 and set compiler flags

If you save this information into example.opts you can pass the opts file to
dunecontrol via the --opts option, e. g.

  dunecontrol --opts=example.opts all

More info
---------

See

     dunecontrol --help

for further options.


The full build system is described in the dune-common/doc/buildsystem (Git version) or under share/doc/dune-common/buildsystem if you installed DUNE!
