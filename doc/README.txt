#! README.TXT

-[BUILDING]-----------------------------------------------------------
To build execute:
$ make clean
$ make -j
this will produce the binary file with name "tsi" in the current directory.
----------------------------------------------------------------------
-[BUILDING WITH MPI]--------------------------------------------------
To build execute:
$ make clean
$ CC=mpicc make -j
this will produce the binary file with name "tsi-mpi" in the current directory.
----------------------------------------------------------------------

-[RUNNING]------------------------------------------------------------
to run the application execute it with one argument, being the configuration file, 
   eg:
 $ ./tsi conf/bigTest.par
or
 $ ./tsi
 if executed without arguments, tsi will look for the configuration file in "conf/tsi_config.par"
note: tsi needs more files besides the configurarion file, those should be indicated in the configuration file.
----------------------------------------------------------------------

--RUNNING WITH MPI----------------------------------------------------
to run the application using MPI, you must have MPI environment allready set.
After that you should run it like this:

$ mpiexec -np NCPUS tsi-mpi conf/tsi_config.par
----------------------------------------------------------------------



-[RUNNING REQUIREMENTS]-----------------------------------------------
To Run the aplication the following files are required:
- executable file
- configuration file
 - Seismic file 
 - wavelet file
 - wells file
 The last three files depend on the reservoir and are not supplied with the application.
 File formats:
 	Configuration file:
		please see configuration files on "conf" directory for examples.
 	Seismic file:
		gslib ASCII format, one value per line.
	 	values are atributed in the folowing order ( grid[x,y,z] )
	-----example--------
	grid[1,0,0]
	...
	grid[XNUMBER,0,0]
	grid[1,1,0]
	...
	grid[XNUMBER,1,0]
	grid[1,2,0]
	...
	grid[XNUMBER,YNUMBER,ZNUMBER]
	--------------------
	wavelet file:
		gslib ASCII format, two values per line: "x value", expecting x in [-A,A] integer ranges
	Wells file:
		gslib ASCII format, 4 values per line: "x y z value".

 The output of the application, is "bestAiCube.out", a grid of acoustic impedances in the same
 gslib format has the SEISMIC grid. the file is saved in the directory indicated by the config file.
----------------------------------------------------------------------



-[REQUIREMENTS]-------------------------------------------------------
Hardware:
 - architecture:
   PC x86 computer (intel, amd ou compatible), preferably SMP, or cc-NUMA.
 - RAM memory
  minimum advised: +/- 2Gbytes for grid with 40 milion points
  (depends of data size and number of parallel simulations)
 - Disk space:
  application itself is fairly small, but datasets are usually quite big.
  minimum advised for testing environment: ~ 20Gb
  recomended for testing and debugging: ~ 100Gb
  minimum production environment: 2 x seismic file size (for final result)

more information about memory requirements & scalability in PERFORMANCE.txt

Software:
 - Operating System:
    Linux (tested e developed in Linux-2.6)
    glibc 2.3.5
 - Compiler
	gcc-3 ou icc (intel compiler)
 - Aditional Tools
	GNU make
	file editor
    openmpi for mpi environment.
note: code is fairly POSIX, should build and work 
on most UNIX systems, and even on cygwin..
although this was not tested yet.
----------------------------------------------------------------------

-[ISSUES]-------------------------------------------------------------
 see BUGS.txt
----------------------------------------------------------------------


---------------------------------------------- info about versions ---
tsi-v5.6 - 2007.02.XX
	- TSI with MPI support.
	- print layers size in log file.
	- print seed in log file.
	- supports windows (single processor mode)
	- new pseudo random number generator (mersenne twister)
	- better error messages when problems are found in input files.

tsi-5.4 - 2006.07.18
	- TSI with MPI support, tested and with great measured performance.
	- tested with grids with up to 120 million points. (32bit mode, with mpi)
	- memory usage is configurable (proporcional to the grid size, can be 9x grid size to 13x grid size)
	- better 64bit support (faster than 32bits, good/valid results up to 34 million points)
	- considerable speed optimizations.
	- reworked random path generation (faster)
	- lots of bugs fixed on DSS

tsi-5.1 - 2006.04.21
	- First release with MPI support, MPI has some performance issues..
	- virtual grid manager support, will allow to configure the number of grids to keep in memory.
	- better logging.


tsi-v4.0 - 2005.08.25
	- TSI multi-threaded, execution time is inversely proporcional to the number of 
	CPUs/THREADS if number of simulations is a multiple of the number of simulation THREADS
	- predicted peak memory use is indicated in the log file,
   		at the beginning of execution.
	- memory consumption was reduced
	- memory consumption depends of number of simulation threads
	WARNING:
	- compiler optimizations bigger than:
			-O3 -march=pentium2 -ffast-math -fomit-frame-pointer
		are unsafe and unreliable.
	- code is not 64bits clean, only x86 is suported.
	- its not amd64 safe.

tsi-v3.2 - 2005.08.10  
	- DSS ported to ansi C.
	- configuration file can be passed has argument
	- memory consumption was greatly reduced: 30% lower
	- execution time was reduced to 85% of tsi-v3.0
	WARNING:
	- compiler optimizations bigger than:
			-O3 -march=pentium2 -ffast-math -fomit-frame-pointer
		are unsafe and unreliable.
	- code is not 64bits clean, only x86 is suported.
	- its not amd64 safe.

tsi-v3.1 - 2005.06.29 
	- automatic translation of DSS from fortran to C
	- asserted quality of results
	- internal release
	
tsi-v3.0 
	- single and friendlier configuration file
	- ported to Linux
	- DSS in fortran77 and bits in fortran90
	- SI in C++

----------------------------------------------------------------------
