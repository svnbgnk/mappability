Mappability
===========

The program is under development and in an alpha-stage. Below you can find examples on how to run the application. If something is not clear, just send me an e-mail: christopher.pockrandt AT fu-berlin DOT de.

How to build the project
------------------------

```sh
   $ git clone git@github.com:cpockrandt/mappability.git --recursive
   $ mkdir mappability-build && cd mappability-build
   $ cmake ../mappability -DCMAKE_BUILD_TYPE=Release
   $ make -j mappability create_index mappability_dump
```

Documentation
-------------

All binaries document possible flags that you can get with ``--help``, e.g. ``./create_index --help``

### Building the index ###

```sh
   $ ./create_index --genome /path/to/fasta.fa(sta) --index /path/to/index/index_name
```

The index is built using secondary memory. If you're getting a runtime error, you're most likely runing out of disk space or quota. You can change the TMPRDIR environment variable TMPRDIR on UNIX systems (and TEMP on Windows).

```sh
   $ export TMPDIR=/somewhere/else/with/more/space
```

### Computing mappability ###

```sh
   $ ./mappability --index /path/to/index/index_name --output /output/of/mappability --errors 2 --length 50 --overlap 10 --mmap --threads 4
```

* overlap: overlap is a parameter that does not have an effect on the result, only on the running time. The value has to be set s.t. `length - overlap >= errors + 2` and `overlap <= length - 1`. This parameter will be removed soon. A good choice (in terms of performance) is: `min(length, 65) * 0.85^(errors - 1)`.

* The maximum mappability value is 255 by default, i.e. if a k-mer occurs more than 255 times, it is still only reported as 255 hits. To increase this threshold, you can add the flag `--high` which increases the upper bound to 65535.

* output: the outputfile will be named `/output/of/mappability_2_50.gmapp8` or `/output/of/mappability_2_50.gmapp16` if `--high` has been set.

* mmap: activates memory mapping. This might speed up the computation a bit but has no effect on the result.

* threads: please not that there is no locking implemented. If you test trivial inputs (i.e. tiny fasta files), turn off parallelization by using only one thread, otherwise threads might interfere and produce wrong results.

### Transforming the result into a readable format ###

The output from the mappability program is written from RAM to disk, i.e. it is the raw std::vector format. To make it more readable in a text editor (and significantly increasing the file size), you can run:

```sh
   $ ./mappability_dump --input /output/of/mappability_2_50.gmapp8 --output /output/of/mappability_2_50.txt
```
