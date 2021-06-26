// dummy proxy for ParSplice driver

// Syntax: mpirun -np P parsplice switch args ...
// Switches:
//   -mdi MDI
//      MDI = "-method LINK -name d -role DRIVER -plugin_path $HOME/LATTE"
//      plugin_path = full path to dir where engine shlibs or scripts are
//   -plugs LMPname LATTEname FITname
//      LMPname = lammps_plug if liblammps_plug.so is in plugin_path dir
//      LATTEname = latte_plug if liblatte_plug.so is in plugin_path dir
//      FITname = fitSNAP if fitSNAP.py is in plugin_path dir
//      NOTE: any of these can be NULL, Nlammps or Nlatte or Nfit must also then be 0
//   -p Nps PSdrivers Nlammps Plammps Nlatte Platte Nfit Pfit (def = P 1 0 0 0 0 0 0)
//      Nps = # of ParSplice procs
//      PSdrivers = # of Nps procs used to manage engine tasks
//      Nlammps = # of LAMMPS instances
//      Plammps = # of MPI tasks per LAMMPS instance
//      Nlatte = # of LATTE instances
//      Platte = # of MPI tasks per LATTE instance
//      Nfit = # of fitSNAP instances
//      Pfit = # of MPI tasks per fitSNAP instance
//      NOTE: P = Nps + Nlammps*Plammps + Nlatte*Platte + Nfit*Pfit is required
//   -w PS Tlammps Dlammps Tlatte Dlatte Tfit Dfit (def = 1 0 0.0 0 0.0 0 0.0)
//      Tlammps = # of LAMMPS tasks to perform
//      Dlammps = time in seconds to delay for a LAMMPS task to run
//      Tlatte = # of LATTE tasks to perform
//      Dlatte = time in seconds to delay for a LATTE task to run
//      Tfit = # of fitSNAP tasks to perform
//      Dfit = time in seconds to delay for a fitSNAP task to run
//   -e Nsize (def = 9)
//      Nsize = # of doubles exchanged at start/end of task
//   -o 0/1
//      0 = minimal output
//      1 = full ouptut (3 lines per work task)
//   -f <latteFile> Default is "latte.in" 

#include <mpi.h>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "mdi.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <sstream>
#include <cstdio>
// --------------------------------------------
// data
// --------------------------------------------

int me, nprocs;			// all of ParSplice
MPI_Comm world;			// MPI_COMM_WORLD for ParSplice

// Data passed to the engine
double dlammps, dlatte, dfit;
int natoms;			// Number of atoms     
int ntypes;			// Number of atom types (elements)
double *coords;			// Coordinates (3*natoms)
double *cell;			// Cell parameters 
char *symbols;			// Symbols for every atom type 
char *fname;			// Latte file name (default in latte.in)
int *types;			// Array that indexes the element type (type(1) = 1 <=> Atom 1 is of type 1)

// Data passed from the enegine
double *venerg;			// Potential energy (eV in LATTE)
double *forces;			// Forces (eV/Ang in LATTE) 

int engineflag;
enum
{ NONE, LAMMPS, LATTE, FITSNAP };

// --------------------------------------------
// functions
// --------------------------------------------

void
error (const char *str)
{
  if (me == 0)
    printf ("ERROR: %s\n", str);
  MPI_Abort (MPI_COMM_WORLD, 1);
}

// MDI driver to communicate with one LAMMPS or LATTE or fitSNAP engine instance
// the 3 engine wrappers are identical for this test, could each be different
//   except for sending of engine-specific time delay
// called by procs within each engine communicator after engine is instantiated
// mpicomm_ptr = ptr to MPI comm the engine instance is running on
// mdicomm = used by all calls to MDI
// class_object = arg passed to MDI_Launch_Plugin()

int
engine_wrapper (void *mpicomm_ptr, MDI_Comm mdicomm, void *class_object)
{
  int me_engine, nprocs_engine;
  MPI_Comm world_engine = *(MPI_Comm *) mpicomm_ptr;
  MPI_Comm_rank (world_engine, &me_engine);
  MPI_Comm_size (world_engine, &nprocs_engine);

  // query name of engine

  int nlen;
  char *engine_name = new char[MDI_NAME_LENGTH];
  MDI_Send_Command ("<NAME", mdicomm);
  MDI_Recv (engine_name, MDI_NAME_LENGTH, MDI_CHAR, mdicomm);
  MPI_Bcast (engine_name, MDI_NAME_LENGTH, MPI_CHAR, 0, world_engine);

  // send setup info to engine
  // data size, time delay
  // test on engineflag for which time delay to send

  // Send the name of the latte file 
  MDI_Send_Command (">FNAME", mdicomm);
  MDI_Send (fname, 20, MDI_CHAR, mdicomm);

  // Send the number of atoms
  natoms = 2;
  MDI_Send_Command (">NATOMS", mdicomm);
  MDI_Send (&natoms, 1, MDI_INT, mdicomm);

  // Send number of types
  ntypes = 2;
  MDI_Send_Command (">NTYPES", mdicomm);
  MDI_Send (&ntypes, 1, MDI_INT, mdicomm);

  // Send atoms symbols/names
  symbols = "C O";
  MDI_Send_Command (">SYMBOLS", mdicomm);
  MDI_Send (symbols, 100, MDI_CHAR, mdicomm);

  // Send a vector of lenght natms with the types
  types = new int[natoms];
  types[0] = 1;
  types[1] = 2;
  MDI_Send_Command (">TYPES", mdicomm);
  MDI_Send (types, natoms, MDI_INT, mdicomm);

  // Send the coordinates 
  coords = new double[3 * natoms];
  coords[0] = 0.0;
  coords[1] = 0.0;
  coords[2] = 0.0;
  coords[3] = 0.8;
  coords[4] = 0.0;
  coords[5] = 0.0;
  MDI_Send_Command (">COORDS", mdicomm);
  MDI_Send (coords, 3 * natoms, MDI_DOUBLE, mdicomm);

  // Send cell parameters
  cell = new double[9];
  cell[0] = 0.0;		//xlo(1) 
  cell[1] = 5.0;		//xhi(1) 
  cell[2] = 0.0;		//xlo(2) 
  cell[3] = 5.0;		//xhi(2)
  cell[4] = 0.0;		//xlo(3)
  cell[5] = 5.0;		//xhi(3) 
  cell[6] = 0.0;		//xy
  cell[7] = 0.0;		//xz 
  cell[8] = 0.0;		//yz 
  MDI_Send_Command (">CELL", mdicomm);
  MDI_Send (cell, 9, MDI_DOUBLE, mdicomm);

  // Print info on engine
  MDI_Send_Command ("ECHO", mdicomm);

  // Allocate forces
  forces = new double[3 * natoms];
  venerg = new double[1];

  // Run LATTE 
  for (int i = 0; i < 5; i++)
    {
      coords[3] = coords[3] + 0.1;	// Stretch the CO molecule
      MDI_Send_Command (">COORDS", mdicomm);
      MDI_Send (coords, 3 * natoms, MDI_DOUBLE, mdicomm);
      MDI_Send_Command ("RUN", mdicomm);

      // Receive forces
      MDI_Send_Command ("<FORCES", mdicomm);
      MDI_Recv (forces, 3 * natoms, MDI_DOUBLE, mdicomm);
      MPI_Bcast (forces, 3 * natoms, MPI_DOUBLE, 0, world_engine);
      printf ("Forces %f %f\n", forces[0], forces[3]);

      // Receive potential energy
      MDI_Send_Command ("<PE", mdicomm);
      MDI_Recv (venerg, 1, MDI_DOUBLE, mdicomm);
      MPI_Bcast (venerg, 1, MPI_DOUBLE, 0, world_engine);
      printf ("Potential Energy %f\n", venerg[0]);
    }

  MDI_Send_Command ("EXIT", mdicomm);

  // clean up
  delete[]engine_name;
  delete[]coords;

  return 0;
}

// --------------------------------------------
// main program
// --------------------------------------------

int
main (int narg, char **args)
{
  // initialize MPI

  MPI_Init (&narg, &args);
  world = MPI_COMM_WORLD;
  MPI_Comm_rank (world, &me);
  MPI_Comm_size (world, &nprocs);

  // parse command line args
  // defaults as listed

  bool initialized_mdi = false;
  char *pluglmp = NULL;
  char *pluglatte = NULL;
  char *plugfit = NULL;
//  char *fname = "latte.in";

  int nparsplice, psdrivers, nlammps, plammps, nlatte, platte, nfit, pfit;
  int tlammps, tlatte, tfit;
  int outflag;

  nparsplice = nprocs;
  nlammps = plammps = 0;
  nlatte = platte = 0;
  nfit = pfit = 0;
  psdrivers = 1;
  tlammps = tlatte = tfit = 0;
  dlammps = dlatte = dfit = 0.0;
  natoms = 100;
  outflag = 1;

  int iarg = 1;
  while (iarg < narg)
    {
      if (strcmp (args[iarg], "-mdi") == 0)
	{
	  if (iarg + 2 > narg)
	    error ("Invalid command line syntax");
	  int err = MDI_Init (&narg, &args);	// this removes -mdi arg
	  if (err)
	    error ("Driver MDI lib not initialized correctly");
	  err = MDI_MPI_get_world_comm (&world);
	  if (err)
	    error ("MDI_MPI_get_world_comm returned error");
	  initialized_mdi = true;
	  //iarg += 2;     // not needed b/c -mdi arg removed
	}
      else if (strcmp (args[iarg], "-plugs") == 0)
	{
	  if (iarg + 4 > narg)
	    error ("Invalid command line syntax");
	  if (strcmp (args[iarg + 1], "NULL") == 0)
	    pluglmp = NULL;
	  else
	    pluglmp = args[iarg + 1];
	  if (strcmp (args[iarg + 2], "NULL") == 0)
	    pluglatte = NULL;
	  else
	    pluglatte = args[iarg + 2];
	  if (strcmp (args[iarg + 3], "NULL") == 0)
	    plugfit = NULL;
	  else
	    plugfit = args[iarg + 3];
	  iarg += 4;
	}
      else if (strcmp (args[iarg], "-p") == 0)
	{
	  if (iarg + 9 > narg)
	    error ("Invalid command line syntax");
	  nparsplice = atoi (args[iarg + 1]);
	  psdrivers = atoi (args[iarg + 2]);
	  nlammps = atoi (args[iarg + 3]);
	  plammps = atoi (args[iarg + 4]);
	  nlatte = atoi (args[iarg + 5]);
	  platte = atoi (args[iarg + 6]);
	  nfit = atoi (args[iarg + 7]);
	  pfit = atoi (args[iarg + 8]);
	  iarg += 9;
	}
      else if (strcmp (args[iarg], "-w") == 0)
	{
	  if (iarg + 7 > narg)
	    error ("Invalid command line syntax");
	  tlammps = atoi (args[iarg + 1]);
	  dlammps = atof (args[iarg + 2]);
	  tlatte = atoi (args[iarg + 3]);
	  dlatte = atof (args[iarg + 4]);
	  tfit = atoi (args[iarg + 5]);
	  dfit = atof (args[iarg + 6]);
	  iarg += 7;
	}
      else if (strcmp (args[iarg], "-e") == 0)
	{
	  if (iarg + 2 > narg)
	    error ("Invalid command line syntax");
	  natoms = atoi (args[iarg + 1]);
	  iarg += 2;
	}
      else if (strcmp (args[iarg], "-o") == 0)
	{
	  if (iarg + 2 > narg)
	    error ("Invalid command line syntax");
	  outflag = atoi (args[iarg + 1]);
	  iarg += 2;
	}
      else if (strcmp (args[iarg], "-f") == 0)
	{
	  if (strcmp (args[iarg + 1], "NULL") == 0)
	    fname = "latte.in";
	  else
	    fname = args[iarg + 1];
	  iarg += 2;
	}
      else
	error ("Unrecognized command line option");
    }

  // error checks

  if (not initialized_mdi)
    error ("-mdi command line option must be used");
  if (!pluglmp && nlammps)
    error ("no LAMMPS plugin, nlammps must be 0");
  if (!pluglatte && nlatte)
    error ("no LATTE plugin, nlatte must be 0");
  if (!plugfit && nfit)
    error ("no FitSNAP plugin, nfit must be 0");

  if (nparsplice < 1)
    error ("nparsplice must be > 0");
  if (psdrivers < 1 || psdrivers > nparsplice)
    error ("psdrivers must be > 0 and <= nparsplice");
  if (nlammps < 0)
    error ("nlammps must be >= 0");
  if (nlatte < 0)
    error ("nlatte must be >= 0");
  if (nfit < 0)
    error ("nfit must be >= 0");
  if (nlammps && plammps <= 0)
    error ("plammps must be >= 1");
  if (nlatte && platte <= 0)
    error ("platte must be >= 1");
  if (nfit && pfit <= 0)
    error ("pfit must be >= 1");

  if (nlammps && dlammps < 0.0)
    error ("dlmp must be >= 0.0");
  if (nlatte && dlatte < 0.0)
    error ("dlatte must be >= 0.0");
  if (nfit && dfit < 0.0)
    error ("dfit must be >= 0.0");
  if (nlammps == 0 && tlammps > 0)
    error ("no LAMMPS instances to run LAMMPS tasks");
  if (nlatte == 0 && tlatte > 0)
    error ("no LATTE instances to run LATTE tasks");
  if (nfit == 0 && tfit > 0)
    error ("no FitSNAP instances to run FitSNAP tasks");
  if (outflag < 0 || outflag > 1)
    error ("invalid outflag setting");

  if ((nparsplice + nlammps * plammps + nlatte * platte + nfit * pfit) !=
      nprocs)
    error ("Nprocs = nparsplice + nlammps*plammps + nlatte*platte + "
	   "nfit*pfit is required");

  // split MPI world comm into sub-comms for the driver and each engine instance
  // each sub-comm is a color index

  int color, me_intra, nprocs_intra;
  MPI_Comm intra_comm;

  if (me < nparsplice)
    color = 0;
  else if (me < nparsplice + nlammps * plammps)
    color = (me - nparsplice) / plammps + 1;
  else if (me < nparsplice + nlammps * plammps + nlatte * platte)
    color = (me - nparsplice - nlammps * plammps) / platte + nlammps + 1;
  else if (me <
	   nparsplice + nlammps * plammps + nlatte * platte + nfit * pfit)
    color =
      (me - nparsplice - nlammps * plammps - nlatte * platte) / pfit +
      nlammps + nlatte + 1;

  MPI_Comm_split (world, color, me, &intra_comm);
  MPI_Comm_rank (intra_comm, &me_intra);
  MPI_Comm_size (intra_comm, &nprocs_intra);

  // print run info

  if (me == 0)
    {
      printf ("ParSplice procs: %d\n", nparsplice);
      printf ("LAMMPS instances and procs: %d %d\n", nlammps,
	      nlammps * plammps);
      printf ("LATTE instances and procs %d %d\n", nlatte, nlatte * platte);
      printf ("FitSNAP instances and procs %d %d\n", nfit, nfit * pfit);
    }

  // start timer

  MPI_Barrier (world);
  double time_start = MPI_Wtime ();

  // driver procs

  if (color == 0)
    {
      engineflag = NONE;

      if (me_intra < psdrivers)
	{

	  // divvy up engine instances and tasks across psdrivers

	  // nlammps_mine = # of LAMMPS instances this proc controls
	  // nlammps_first = first instance this proc controls (0 to Nlammps-1)
	  // tlammps_mine = # of LAMMPS tasks this proc runs on its instances
	  // tlammps_first = first LAMMPS task this proc runs on its instances

	  // nlatte_mine = # of LATTE instances this proc controls
	  // nlatte_first = first instance this proc controls (0 to Nlatte-1)
	  // tlatte_mine = # of LATTE tasks this proc runs on its instances
	  // tlatte_first = first LATTE task this proc runs on its instances

	  // nfit_mine = # of fitSNAP instances this proc controls
	  // nfit_first = first instance this proc controls (0 to Nfit-1)
	  // tfit_mine = # of fitSNAP tasks this proc runs on its instances
	  // tfit_first = first fitSNAP task this proc runs on its instances

	  // needed to avoid round-off issues in t lammps/latte/fit mine

	  //double eps = 0.0;
	  double eps = 1.0e-10;

	  int nlammps_first =
	    static_cast < int >(1.0 * me_intra / psdrivers * nlammps);
	  int nlammps_next =
	    static_cast < int >(1.0 * (me_intra + 1) / psdrivers * nlammps);
	  int nlammps_mine = nlammps_next - nlammps_first;

	  int tlammps_first =
	    static_cast <
	    int >(1.0 * nlammps_first / nlammps * tlammps + eps);
	  int tlammps_next =
	    static_cast < int >(1.0 * nlammps_next / nlammps * tlammps + eps);
	  int tlammps_mine = tlammps_next - tlammps_first;

	  int nlatte_first =
	    static_cast < int >(1.0 * me_intra / psdrivers * nlatte);
	  int nlatte_next =
	    static_cast < int >(1.0 * (me_intra + 1) / psdrivers * nlatte);
	  int nlatte_mine = nlatte_next - nlatte_first;
	  int tlatte_first =
	    static_cast < int >(1.0 * nlatte_first / nlatte * tlatte + eps);
	  int tlatte_next =
	    static_cast < int >(1.0 * nlatte_next / nlatte * tlatte + eps);
	  int tlatte_mine = tlatte_next - tlatte_first;

	  int nfit_first =
	    static_cast < int >(1.0 * me_intra / psdrivers * nfit);
	  int nfit_next =
	    static_cast < int >(1.0 * (me_intra + 1) / psdrivers * nfit);
	  int nfit_mine = nfit_next - nfit_first;
	  int tfit_first =
	    static_cast < int >(1.0 * nfit_first / nfit * tfit + eps);
	  int tfit_next =
	    static_cast < int >(1.0 * nfit_next / nfit * tfit + eps);
	  int tfit_mine = tfit_next - tfit_first;

	  // each PS proc triggers its engines to perform tasks
	  // loops until all this proc's tasks in each category are done
	  // driver/engine comm is done within world comm

	  int itask, iproc, flag;
	  int complete_lammps, next_lammps;
	  int complete_latte, next_latte;
	  int complete_fit, next_fit;
	  MPI_Status status;

	  // assign initial LAMMPS tasks to each of my LAMMPS instances
	  // itask = 0 if no tasks left to perform

	  for (int i = 0; i < nlammps_mine; i++)
	    {
	      if (i < tlammps_mine)
		itask = tlammps_first + i + 1;
	      else
		itask = 0;
	      iproc = nparsplice + (nlammps_first + i) * plammps;
	      MPI_Send (&itask, 1, MPI_INT, iproc, 0, world);
	    }

	  next_lammps =
	    nlammps_mine < tlammps_mine ? nlammps_mine + 1 : tlammps_mine + 1;
	  complete_lammps = 0;

	  // assign initial LATTE tasks to each of my LATTE instances
	  // itask = 0 if no tasks left to perform

	  for (int i = 0; i < nlatte_mine; i++)
	    {
	      if (i < tlatte_mine)
		itask = tlatte_first + i + 1;
	      else
		itask = 0;
	      iproc =
		nparsplice + nlammps * plammps + (nlatte_first + i) * platte;
	      MPI_Send (&itask, 1, MPI_INT, iproc, 0, world);
	    }

	  next_latte =
	    nlatte_mine < tlatte_mine ? nlatte_mine + 1 : tlatte_mine + 1;
	  complete_latte = 0;

	  // assign initial fitSNAP tasks to each of my fitSNAP instances
	  // itask = 0 if no tasks left to perform

	  for (int i = 0; i < nfit_mine; i++)
	    {
	      if (i < tfit_mine)
		itask = tfit_first + i + 1;
	      else
		itask = 0;
	      iproc = nparsplice + nlammps * plammps + nlatte * platte +
		(nfit_first + i) * pfit;
	      MPI_Send (&itask, 1, MPI_INT, iproc, 0, world);
	    }

	  next_fit = nfit_mine < tfit_mine ? nfit_mine + 1 : tfit_mine + 1;
	  complete_fit = 0;

	  // loop until all LAMMPS and LATTE and fitSNAP tasks are completed
	  // recv a done message
	  // assign a new task to that instance
	  // itask = 0 if no tasks left to perform

	  while (complete_lammps < tlammps_mine
		 || complete_latte < tlatte_mine || complete_fit < tfit_mine)
	    {
	      MPI_Recv (&flag, 1, MPI_INT, MPI_ANY_SOURCE, 0, world, &status);
	      iproc = status.MPI_SOURCE;

	      if (iproc < nparsplice + nlammps * plammps)
		{
		  complete_lammps++;
		  itask = next_lammps++;
		  if (itask > tlammps_mine)
		    itask = 0;
		  else
		    itask += tlammps_first;
		}
	      else if (iproc <
		       nparsplice + nlammps * plammps + nlatte * platte)
		{
		  complete_latte++;
		  itask = next_latte++;
		  if (itask > tlatte_mine)
		    itask = 0;
		  else
		    itask += tlatte_first;
		}
	      else
		{
		  complete_fit++;
		  itask = next_fit++;
		  if (itask > tfit_mine)
		    itask = 0;
		  else
		    itask += tfit_first;
		}

	      MPI_Send (&itask, 1, MPI_INT, iproc, 0, world);
	    }
	}

      // LAMMPS engine procs

    }
  else if (color <= nlammps)
    {
      engineflag = LAMMPS;
      int instance = (me - nparsplice) / plammps + 1;

      //if (outflag && me_intra == 0)
      //printf("LAMMPS instance %d of size %d begins on proc %d\n",
      //       instance,nprocs_intra,me);

      int itask, flag, driver;
      MPI_Status status;
      char cmdline_args[64] =
	"-mdi \"-name LMP -role ENGINE -method LINK\" -foo bar";

      while (1)
	{

	  // engine proc 0 waits for message from driver to perform a task or exit
	  // if itask == 0, done

	  if (me_intra == 0)
	    {
	      MPI_Recv (&itask, 1, MPI_INT, MPI_ANY_TAG, 0, world, &status);
	      driver = status.MPI_SOURCE;
	    }
	  MPI_Bcast (&itask, 1, MPI_INT, 0, intra_comm);
	  if (itask == 0)
	    break;

	  // instantiate the LAMMPS lib
	  // plugin_name = "lammps_plug" if shlib is liblammps_plug.so
	  // cmdline_args = NULL for now, eventually args for LAMMPS
	  // intra_comm = MPI comm this instance will run on
	  // engine_wrapper = callback func, invoked after LAMMPS is running
	  //   it contains logic to be the MDI driver for LAMMPS engine
	  // final arg = optional class pointer if callback func is a class method

	  if (outflag && me_intra == 0)
	    printf ("LAMMPS task %d starting on LAMMPS instance %d "
		    "on engine proc %d %d from driver proc %d\n",
		    itask, instance, me_intra, me, driver);

	  int err = MDI_Launch_plugin (pluglmp, cmdline_args,
				       &intra_comm, engine_wrapper, NULL);
	  if (err)
	    error ("MDI_Launch_plugin failed");

	  if (me_intra == 0)
	    {
	      if (outflag)
		printf ("LAMMPS task %d completed on LAMMPS instance %d\n",
			itask, instance);
	      flag = 1;
	      MPI_Send (&flag, 1, MPI_INT, driver, 0, world);
	    }
	}

      // LATTE engine procs

    }
  else if (color <= nlammps + nlatte)
    {
      engineflag = LATTE;
      int instance = (me - nparsplice - nlammps * plammps) / platte + 1;

      if (outflag && me_intra == 0)
	printf ("LATTE instance %d of size %d begins on proc %d\n",
		instance, nprocs_intra, me);

      int itask, flag, driver;
      MPI_Status status;
      char cmdline_args[64] =
	"-mdi \"-name LATTE -role ENGINE -method LINK\" -foo bar";

      while (1)
	{

	  // engine proc 0 waits for message from driver to perform a task or exit
	  // if itask == 0, done

	  if (me_intra == 0)
	    {
	      MPI_Recv (&itask, 1, MPI_INT, MPI_ANY_TAG, 0, world, &status);
	      driver = status.MPI_SOURCE;
	    }
	  MPI_Bcast (&itask, 1, MPI_INT, 0, intra_comm);
	  if (itask == 0)
	    break;

	  // instantiate the LATTE lib
	  // plugin_name = "latte_plug" if shlib is liblatte_plug.so
	  // cmdline_args = NULL for now, eventually args for LATTE
	  // intra_comm = MPI comm this instance will run on
	  // engine_wrapper = callback func, invoked after LATTE is running
	  //   it contains logic to be the MDI driver for LATTE engine
	  // final arg = optional class pointer if callback func is a class method

	  if (outflag && me_intra == 0)
	    printf ("LATTE task %d starting on LATTE instance %d "
		    "on engine proc %d %d from driver proc %d\n",
		    itask, instance, me_intra, me, driver);

	  int err = MDI_Launch_plugin (pluglatte, cmdline_args,
				       &intra_comm, engine_wrapper, NULL);
	  if (err)
	    error ("MDI_Launch_plugin failed");

	  if (me_intra == 0)
	    {
	      if (outflag)
		printf ("LATTE task %d completed on LATTE instance %d\n",
			itask, instance);
	      flag = 1;
	      MPI_Send (&flag, 1, MPI_INT, driver, 0, world);
	    }
	}

      // fitSNAP engine procs

    }
  else if (color <= nlammps + nlatte + nfit)
    {
      engineflag = FITSNAP;
      int instance =
	(me - nparsplice - nlammps * plammps - nlatte * platte) / pfit + 1;

      if (outflag && me_intra == 0)
	printf ("FitSNAP instance %d of size %d begins on proc %d\n",
		instance, nprocs_intra, me);

      int itask, flag, driver;
      MPI_Status status;
      char cmdline_args[64] =
	"-mdi \"-name FitSNAP -role ENGINE -method LINK\" -foo bar";

      while (1)
	{

	  // engine proc 0 waits for message from driver to perform a task or exit
	  // if itask == 0, done

	  if (me_intra == 0)
	    {
	      MPI_Recv (&itask, 1, MPI_INT, MPI_ANY_TAG, 0, world, &status);
	      driver = status.MPI_SOURCE;
	    }
	  MPI_Bcast (&itask, 1, MPI_INT, 0, intra_comm);
	  if (itask == 0)
	    break;

	  // instantiate the fitSNAP Python script
	  // plugin_name = "fitSNAP" if script is fitSNAP.py
	  // cmdline_args = NULL for now, eventually args for fitSNAP
	  // intra_comm = MPI comm this instance will run on
	  // engine_wrapper = callback func, invoked after fitSNAP is running
	  //   it contains logic to be the MDI driver for fitSNAP engine
	  // final arg = optional class pointer if callback func is a class method

	  if (outflag && me_intra == 0)
	    printf ("FitSNAP task %d starting on FitSNAP instance %d "
		    "on engine proc %d %d from driver proc %d\n",
		    itask, instance, me_intra, me, driver);

	  int err = MDI_Launch_plugin (plugfit, cmdline_args,
				       &intra_comm, engine_wrapper, NULL);
	  if (err)
	    error ("MDI_Launch_plugin failed");

	  if (me_intra == 0)
	    {
	      if (outflag)
		printf ("FitSNAP task %d completed on fitSNAP instance %d\n",
			itask, instance);
	      flag = 1;
	      MPI_Send (&flag, 1, MPI_INT, driver, 0, world);
	    }
	}
    }

  // stop timer

  MPI_Barrier (world);
  double time_stop = MPI_Wtime ();

  if (me == 0)
    printf ("Total run time (secs): %g\n", time_stop - time_start);

  // shut down MPI

  MPI_Barrier (world);
  MPI_Finalize ();
}
