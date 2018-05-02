Mappability
===========

How to build the project
------------------------

```sh
   $ git clone git@github.com:cpockrandt/mappability.git --recursive
   $ mkdir mappability-build && cd mappability-build
   $ cmake ../mappability -DCMAKE_BUILD_TYPE=Release
   $ make -j
```

Documentation
-------------

All binaries document possible flags that you can get with ``--help``, e.g. ``./create_index --help``

The index is built using secondary memory. If you're getting a runtime error, you're most likely runing out of disk space or quota. You can change the TMPRDIR environment variable TMPRDIR on UNIX systems (and TEMP on Windows).

```sh
   $ export TMPDIR=/somewhere/else/with/more/space
```
