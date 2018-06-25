/*
*******************************************************************************

EPANET.C -- Hydraulic & Water Quality Simulator for Water Distribution Networks

AUTHOR:     M. Lopez-Ibanez (M. L\'opez-Ib\'a\~nez)
            Napier University, Edinburgh, UK.

$LastChangedBy: manu $ $Revision: 174 $
$LastChangedDate: 2008-09-10 16:17:40 +0200 (Wed, 10 Sep 2008) $

---------------------------------------------------------------------

 Copyright (c) 2005, 2006, 2007 Manuel Lopez-Ibanez
 TeX: \copyright 2005, 2006, 2007 Manuel L\'opez-Ib\'a\~nez

 This program is free software (software libre); you can redistribute
 it and/or modify it under the terms of version 2 of the GNU General
 Public License version as published by the Free Software Foundation.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, you can obtain a copy of the GNU
 General Public License at:
                http://www.gnu.org/copyleft/gpl.html
 or by writing to:
          Free Software Foundation, Inc., 59 Temple Place,
                Suite 330, Boston, MA 02111-1307 USA

----------------------------------------------------------------------
----------------------------------------------------------------------

Original public domain sources (no license) from:

VERSION:    2.00
DATE:       5/30/00
            9/7/00
            10/25/00
            3/1/01
            11/19/01
AUTHOR:     L. Rossman
            US EPA - NRMRL

EPANET performs extended period hydraulic and water quality analysis of
looped, pressurized piping networks. The program consists of the
following code modules:

    EPANET.C  -- main module providing supervisory control
    INPUT1.C  -- controls processing of input data
    INPUT2.C  -- reads data from input file
    INPUT3.C  -- parses individual lines of input data
    INPFILE.C -- saves modified input data to a text file
    RULES.C   -- implements rule-based control of piping system
    HYDRAUL.C -- computes extended period hydraulic behavior
    QUALITY.C -- tracks transport & fate of water quality
    OUTPUT.C  -- handles transfer of data to and from binary files
    REPORT.C  -- handles reporting of results to text file
    SMATRIX.C -- sparse matrix linear equation solver routines
    MEMPOOL.C -- memory allocation routines
    HASH.C    -- hash table routines

The program can be compiled as either a stand-alone console application
or as a dynamic link library (DLL) of function calls depending on whether
the macro identifier 'DLL' is defined or not.

See TOOLKIT.H for function prototypes of exported DLL functions
See FUNCS.H for prototypes of all other functions
See TYPES.H for declaration of global constants and data structures
See VARS.H for declaration of global variables
See TEXT.H for declaration of all string constants
See ENUMSTXT.H for assignment of string constants to enumerated types

The following naming conventions are used in all modules of this program:
1. Names of exportable functions in the DLL begin with the "EN" prefix.
2. All other function names are lowercase.
3. Global variable names begin with an uppercase letter.
4. Local variable names are all lowercase.
5. Declared constants and enumerated values defined in TYPES.H are
   all uppercase.
6. String constants defined in TEXT.H begin with a lower case character
   followed by an underscore and then all uppercase characters (e.g.
   t_HEADLOSS)

--------------------------------------------------------------------------

This is the main module of the EPANET program. It uses a series of
functions, all beginning with the letters EN, to control program behavior.
See the main() and ENepanet() functions below for the simplest example of
these.

This module calls the following functions that reside in other modules:
   RULES.C
     initrules()
     allocrules()
     closerules()
   INPUT1.C
     getdata()
     initreport()
   INPUT2.C
     netsize()
     setreport()
   HYDRAUL.C
     openhyd()
     inithyd()
     runhyd()
     nexthyd()
     closehyd()
     resistance()
     tankvolume()
     getenergy()
     setlinkstatus()
     setlinksetting()
   QUALITY.C
     openqual()
     initqual()
     runqual()
     nextqual()
     stepqual()
     closequal()
   REPORT.C
     writeline()
     writelogo()
     writereport()
   HASH.C
     HTcreate()
     HTfind()
     HTfree()

The macro ERRCODE(x) is defined in TYPES.H. It says if the current
value of the error code variable (errcode) is not fatal (< 100) then
execute function x and set the error code equal to its return value.

*******************************************************************************
*/

//#define LINUX // Un-comment this line to compile in Linux
#define DLL     /* Un-comment this line to compile as a DLL */
#ifdef DLL
#ifndef LINUX
#include <windows.h>

/*** Updated 9/7/00 ***/
#include <float.h>

#endif // #ifndef LINUX
#endif // #ifdef  DLL


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include "toolkit.h"
#include "hash.h"
#include "text.h"
#include "types.h"
#include "enumstxt.h"
#include "funcs.h"
#define  EXTERN
#include "vars.h"

void (* viewprog) (char *);     /* Pointer to progress viewing function */
#if 0
static int   openhydfile(void);                /* Opens hydraulics file      */
static int   openoutfile(void);                /* Opens binary output file   */
#endif
static char*   geterrmsg(int, char *, size_t); /* Gets text of error message */

/*
  ----------------------------------------------------------------
  Entry point used to compile a Windows DLL
  ----------------------------------------------------------------
*/

#ifndef LINUX
#ifdef DLL
#pragma argsused
int WINAPI DllEntryPoint(HINSTANCE hinst, unsigned long reason, void* reserved)
{
    viewprog = NULL;
    return 1;
}
#endif

#endif // #ifndef LINUX

/*
  ----------------------------------------------------------------
  Entry point used to compile a stand-alone executable.
  ----------------------------------------------------------------
*/

#ifndef DLL
int   main(int argc, char *argv[])
/*--------------------------------------------------------------
**  Input:   argc    = number of command line arguments
**           *argv[] = array of command line arguments
**  Output:  none
**  Purpose: main program segment
**
**  Command line for stand-alone operation is:
**    progname f1  f2  f3
**  where progname = name of executable this code was compiled to,
**  f1 = name of input file, f2 = name of report file, and
**  f3 = name of binary output file (optional).
**--------------------------------------------------------------
*/
{
    char *f1,*f2,*f3;
    char blank[] = "";
    int  errcode;

/* Check for proper number of command line arguments */
    if (argc < 3) writecon(FMT03);
    else
    {

        /* Call the main control function */
        f1 = argv[1];
        f2 = argv[2];
        if (argc > 3) f3 = argv[3];
        else          f3 = blank;
        writecon(FMT01);
        errcode = ENepanet(f1,f2,f3,NULL);
        if (errcode > 0) writecon(FMT11);
        /* Negative errorcode replaces Warnflag usage.  */
        else if (errorcode < 0) writecon(FMT10);
        else writecon(FMT09);
    }
    return(0);
}                                       /* End of main */
#endif


/*
  ----------------------------------------------------------------
  Functions for opening & closing the EPANET system
  ----------------------------------------------------------------
*/


/*** updated 3/1/01 ***/
int DLLEXPORT ENepanet(char *f1, char *f2, char *f3, void (*pviewprog) (char *))

/*------------------------------------------------------------------------
**   Input:   f1 = pointer to name of input file
**            f2 = pointer to name of report file
**            f3 = pointer to name of binary output file
**            pviewprog = see note below
**   Output:  none
**  Returns: error code
**  Purpose: runs a complete EPANET simulation
**
**  The pviewprog() argument is a pointer to a callback function
**  that takes a character string (char *) as its only parameter.
**  The function would reside in and be used by the calling
**  program to display the progress messages that EPANET generates
**  as it carries out its computations. If this feature is not
**  needed then the argument should be NULL.
**-------------------------------------------------------------------------
*/
{
    int  errcode = 0;
    viewprog = pviewprog;
    ERRCODE(ENopen(f1,f2,f3));
#if 0
    if (Hydflag != USE) 
#endif
        ERRCODE(ENsolveH());
#if 0
    ERRCODE(ENsolveQ());
    ERRCODE(ENreport());
#endif
    ENclose();
    return(errcode);
}


int DLLEXPORT ENopen(char *f1, char *f2, char *f3)
/*----------------------------------------------------------------
**  Input:   f1 = pointer to name of input file
**           f2 = pointer to name of report file
**           f3 = pointer to name of binary output file
**  Output:  none
**  Returns: error code
**  Purpose: opens EPANET input file & reads in network data
**----------------------------------------------------------------
*/
{
    int  errcode = 0;

/*** Updated 9/7/00 ***/
/* Reset math coprocessor */
#ifndef LINUX
#ifdef DLL
    _fpreset();
#endif
#endif // #ifndef LINUX

/* Set system flags */
    Openflag  = FALSE;
#if 0
    OpenQflag = FALSE;
    SaveQflag = FALSE;

/*** Updated 9/7/00 ***/
    Messageflag = TRUE;
/* If binary output file being used, then   */
/* do not write full results to Report file */
/* (use it only for status reports).        */
    Rptflag = 0;
    if (strlen(f3) == 0) Rptflag = 1;

/*** Updated 9/7/00 ***/
/*** Previous code segment ignored. ***/
/*** Rptflag now always set to 1.   ***/
    Rptflag = 1;
#endif

/* Initialize global pointers to NULL. */
    initpointers();

/* Open input & report files */
    ERRCODE(openfiles(f1,f2,f3));
    if (errcode > 0)
    {
#if 0
        errmsg(errcode);
#endif
        return(errcode);
    }
#if 0
    writelogo();
#endif

/* Find network size & allocate memory for data */
    writecon(FMT02);
#ifndef LINUX
    writewin(FMT100);
#endif
    ERRCODE(netsize());
    ERRCODE(allocdata());

/* Retrieve input data */
    ERRCODE(getdata());

/* Free temporary linked lists used for Patterns & Curves */
    freeTmplist(Patlist);
    freeTmplist(Curvelist);
#if 0
/* If using previously saved hydraulics then open its file */
    if (Hydflag == USE) ERRCODE(openhydfile());
    Statflag  = FALSE;
#endif

/* Write input summary to report file */
    if (!errcode)
    {
#if 0
        if (Summaryflag) writesummary();
        writetime(FMT104);
#endif
        Openflag = TRUE;
    }
#if 0
    else errmsg(errcode);
#endif
    return(errcode);
}

#if 0
int DLLEXPORT ENsaveinpfile(char *filename)
/*----------------------------------------------------------------
**  Input:   filename = name of INP file
**  Output:  none
**  Returns: error code
**  Purpose: saves current data base to file
**----------------------------------------------------------------
*/
{
    if (!Openflag) return(102);
    return(saveinpfile(filename));
}
#endif

int DLLEXPORT ENclose()
/*----------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Returns: error code
**  Purpose: frees all memory & files used by EPANET
**----------------------------------------------------------------
*/
{
#if 0
    if (Openflag) writetime(FMT105);
#endif
    freedata();
    if (InFile  != NULL) fclose(InFile);
#if 0
    if (RptFile != NULL) fclose(RptFile);
    if (HydFile != NULL) fclose(HydFile);
    if (OutFile != NULL) fclose(OutFile);
#endif
    Openflag  = FALSE;
#if 0
    OpenQflag = FALSE;
    SaveQflag = FALSE;
#endif
    return(0);
}


/*
  ----------------------------------------------------------------
  Functions for running a hydraulic analysis
  ----------------------------------------------------------------
*/


int DLLEXPORT ENsolveH()
/*----------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Returns: error code
**  Purpose: solves for network hydraulics in all time periods
**----------------------------------------------------------------
*/
{
    int  errcode;
    long t, tstep;
    ENsimulation_t simulation;

    /* Open hydraulics solver */
    errcode = ENopenH (&simulation);
    if (!errcode)
    {
        /* Initialize hydraulics */
        errcode = ENinitH (simulation, EN_SAVE);
        writecon(FMT14);

        /* Analyze each hydraulic period */
        if (!errcode) do
        {
#if 0
            ENclockstring_t atime;
            /* Display progress message */

/*** Updated 6/24/02 ***/
            sprintf(Msg,"%-10s",clocktime (&atime, simulation->Htime));

            writecon(Msg);
            sprintf(Msg,FMT101, atime);
#ifndef LINUX
            writewin(Msg);
#endif
#endif
            /* Solve for hydraulics & advance to next time period */
            ERRCODE(ENrunH (simulation, &t));
            ERRCODE(ENnextH(simulation, &tstep));

/*** Updated 6/24/02 ***/
            writecon("\b\b\b\b\b\b\b\b\b\b");
        }
        while (tstep > 0);
    }

/* Close hydraulics solver */

/*** Updated 6/24/02 ***/
    writecon("\b\b\b\b\b\b\b\b                     ");

    /* Return negative errcode if warnings but not errors.  */
    errcode = (errcode >= simulation->Warnflag) 
        ? errcode : - simulation->Warnflag;
    ENcloseH (simulation);
    return(errcode);
}

#if 0
int DLLEXPORT ENsaveH()
/*----------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Returns: error code
**  Purpose: saves hydraulic results to binary file.
**
**  Must be called before ENreport() if no WQ simulation made.
**  Should not be called if ENsolveQ() will be used.
**----------------------------------------------------------------
*/
{
    char tmpflag;
    int  errcode;

/* Check if hydraulic results exist */
    if (!simulation->SaveHflag) return(104);

/* Temporarily turn off WQ analysis */
    tmpflag = Qualflag;
    Qualflag = NONE;
/* Call WQ solver to simply transfer results */
/* from Hydraulics file to Output file at    */
/* fixed length reporting time intervals.    */
    errcode = ENsolveQ();
/* Restore WQ analysis option */
    Qualflag = tmpflag;
    if (errcode) errmsg(errcode);
    return(errcode);
}
#endif

int DLLEXPORT ENopenH(ENsimulation_t *sim)
/*----------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Returns: error code
**  Purpose: sets up data structures for hydraulic analysis
**----------------------------------------------------------------
*/
{
    int  errcode = 0;

/* Check that input data exists */
    if (!Openflag) return(102);

#if 0
/* Check that previously saved hydraulics file not in use */
    if (Hydflag == USE) return(107);
#endif

    *sim = NULL;
/* Open hydraulics solver */
    *sim = openhyd();

/* FIXME   if (*sim == NULL) errmsg(errcode); 
        (*sim)->SaveHflag = FALSE;
*/
    return errcode;
}


/*** Updated 3/1/01 ***/
int DLLEXPORT ENinitH(ENsimulation_t sim, int flag)
/*----------------------------------------------------------------
**  Input:   flag = 2-digit flag where 1st (left) digit indicates
**                  if link flows should be re-initialized (1) or
**                  not (0) and 2nd digit indicates if hydraulic
**                  results should be saved to file (1) or not (0)
**  Output:  none
**  Returns: error code
**  Purpose: initializes hydraulic analysis
**----------------------------------------------------------------
*/
{
    int errcode = 0;
    int sflag, fflag;

    /* Check that hydraulics solver was opened */
    if (sim == NULL) return(103);

/* Reset status flags */
    sim->SaveHflag = FALSE;
    sim->Warnflag  = FALSE;

    sim->Nwarnings = 0; /* Reset number of warnings.  */
    sim->TotalSystemDemand = 0.0;
    sim->TotalSystemInflow = 0.0;
    sim->TotalSystemLeakage = 0.0;

/* Get values of save-to-file flag and reinitialize-flows flag */
    fflag = flag/EN_INITFLOW;
    sflag = flag - fflag*EN_INITFLOW;


/* Open hydraulics file */
    sim->Saveflag = FALSE;
    if (sflag > 0)
    {
#if 0
        errcode = openhydfile();
#endif
        if (!errcode) sim->Saveflag = TRUE;
#if 0
        else errmsg(errcode);
#endif
    }

/* Initialize hydraulics */
    inithyd (sim, fflag);
#if 0
    if (Statflag > 0) writeheader(STATHDR,0);
#endif
    return(errcode);
}


int DLLEXPORT ENrunH(ENsimulation_t sim, long *t)
/*----------------------------------------------------------------
**  Input:   none (no need to supply a value for *t)
**  Output:  *t = current simulation time (seconds)
**  Returns: error/warning code
**  Purpose: solves hydraulics for conditions at time t.
**
**  This function is used in a loop with ENnextH() to run
**  an extended period hydraulic simulation.
**  See ENsolveH() for an example.
**----------------------------------------------------------------
*/
{
    int errcode;
    *t = 0;
    if (sim == NULL) return(103);
    errcode = runhyd (sim, t);
#if 0
    if (errcode) errmsg(errcode);
#endif
    return(errcode);
}


int DLLEXPORT ENnextH(ENsimulation_t sim, long *tstep)
/*----------------------------------------------------------------
**  Input:   none (no need to supply a value for *tstep)
**  Output:  *tstep = time (seconds) until next hydraulic event
**                    (0 marks end of simulation period)
**  Returns: error code
**  Purpose: determines time until next hydraulic event.
**
**  This function is used in a loop with ENrunH() to run
**  an extended period hydraulic simulation.
**  See ENsolveH() for an example.
**----------------------------------------------------------------
*/
{
    int errcode;
    *tstep = 0;
    if (sim == NULL) return(103);
    errcode = nexthyd (sim, tstep);
#if 0
    if (errcode) errmsg(errcode);
    else
#endif 
        if (sim->Saveflag && *tstep == 0) sim->SaveHflag = TRUE;
    return(errcode);
}


int DLLEXPORT ENcloseH(ENsimulation_t sim)
/*----------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Returns: error code
**  Purpose: frees data allocated by hydraulics solver
**----------------------------------------------------------------
*/
{
    if (!Openflag) return(102);
    if (sim == NULL) return(103);
    closehyd (sim);

    return(0);
}

#if 0
int DLLEXPORT ENsavehydfile(char *filename)
/*----------------------------------------------------------------
**  Input:   filename = name of file
**  Output:  none
**  Returns: error code
**  Purpose: copies binary hydraulics file to disk
**----------------------------------------------------------------
*/
{
    FILE *f;
    int   c;

/* Check that hydraulics results exist */
    if (HydFile == NULL || !simulation->SaveHflag) return(104);

/* Open file */
    if ( (f = fopen(filename,"w+b")) == NULL) return(305);

/* Copy from HydFile to f */
    fseek(HydFile, 0, SEEK_SET);
    while ( (c = fgetc(HydFile)) != EOF) fputc(c, f);
    fclose(f);
    return(0);
}
#endif

#if 0
int DLLEXPORT ENusehydfile(char *filename)
/*----------------------------------------------------------------
**  Input:   filename = name of file
**  Output:  none
**  Returns: error code
**  Purpose: opens previously saved binary hydraulics file
**----------------------------------------------------------------
*/
{
    int errcode;

/* Check that input data exists & hydraulics system closed */
    if (!Openflag) return(102);
    if (simulation != NULL) return(108);

/* Try to open hydraulics file */
    strncpy(HydFname, filename, MAXFNAME);
    Hydflag = USE;
    simulation->SaveHflag = TRUE;
    errcode = openhydfile();

/* If error, then reset flags */
    if (errcode)
    {
        strcpy(HydFname, "");
        Hydflag = SCRATCH;
        simulation->SaveHflag = FALSE;
    }
    return(errcode);
}
#endif

/*
  ----------------------------------------------------------------
  Functions for running a WQ analysis
  ----------------------------------------------------------------
*/

#if 0
int DLLEXPORT ENsolveQ()
/*----------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Returns: error code
**  Purpose: solves for network water quality in all time periods
**----------------------------------------------------------------
*/
{
    int  errcode;
    long t, tstep;

/* Open WQ solver */
    errcode = ENopenQ();
    if (!errcode)
    {
        /* Initialize WQ */
        errcode = ENinitQ(EN_SAVE);
        if (Qualflag) writecon(FMT15);
        else
        {
            writecon(FMT16);
#ifndef LINUX
            writewin(FMT103);
#endif
        }

        /* Analyze each hydraulic period */
        if (!errcode) do
        {
            ENclockstring_t atime;
            /* Display progress message */

/*** Updated 6/24/02 ***/
            sprintf(Msg,"%-10s",clocktime (&atime, simulation->Htime));

            writecon(Msg);
            if (Qualflag)
            {
                sprintf(Msg,FMT102, atime);
#ifndef LINUX
                writewin(Msg);
#endif
            }

            /* Retrieve current network solution & update WQ to next time period */
            tstep = 0;
            ERRCODE(ENrunQ(&t));
            ERRCODE(ENnextQ(&tstep));

/*** Updated 6/24/02 ***/
            writecon("\b\b\b\b\b\b\b\b\b\b");

        }  while (tstep > 0);

    }

/* Close WQ solver */

/*** Updated 6/24/02 ***/
    writecon("\b\b\b\b\b\b\b\b                     ");
    ENcloseQ();
    return(errcode);
}
#endif
#if 0
int DLLEXPORT ENopenQ()
/*----------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Returns: error code
**  Purpose: sets up data structures for WQ analysis
**----------------------------------------------------------------
*/
{
    int errcode = 0;

/* Check that hydraulics results exist */
    OpenQflag = FALSE;
    SaveQflag = FALSE;
    if (!Openflag) return(102);
    if (!simulation->SaveHflag) return(104);

/* Open WQ solver */
    ERRCODE(openqual());
    if (!errcode) OpenQflag = TRUE;
    else errmsg(errcode);
    return(errcode);
}
#endif
#if 0
int DLLEXPORT ENinitQ(int saveflag)
/*----------------------------------------------------------------
**  Input:   saveflag = EN_SAVE (1) if results saved to file,
**                      EN_NOSAVE (0) if not
**  Output:  none
**  Returns: error code
**  Purpose: initializes WQ analysis
**----------------------------------------------------------------
*/
{
    int errcode = 0;
    if (!OpenQflag) return(105);
    initqual();
    SaveQflag = FALSE;
    simulation->Saveflag = FALSE;
    if (saveflag)
    {
#if 0
        errcode = openoutfile();
#endif
        if (!errcode) simulation->Saveflag = TRUE;
    }
    return(errcode);
}
#endif

#if 0
int DLLEXPORT ENrunQ(long *t)
/*----------------------------------------------------------------
**  Input:   none (no need to supply a value for *t)
**  Output:  *t = current simulation time (seconds)
**  Returns: error code
**  Purpose: retrieves hydraulic & WQ results at time t.
**
**  This function is used in a loop with ENnextQ() to run
**  an extended period WQ simulation. See ENsolveQ() for
**  an example.
**----------------------------------------------------------------
*/
{
    int errcode;
    *t = 0;
    if (!OpenQflag) return(105);
    errcode = runqual(t);
    if (errcode) errmsg(errcode);
    return(errcode);
}
#endif

#if 0
int DLLEXPORT ENnextQ(long *tstep)
/*----------------------------------------------------------------
**  Input:   none (no need to supply a value for *tstep)
**  Output:  *tstep = time (seconds) until next hydraulic event
**                    (0 marks end of simulation period)
**  Returns: error code
**  Purpose: advances WQ simulation to next hydraulic event.
**
**  This function is used in a loop with ENrunQ() to run
**  an extended period WQ simulation. See ENsolveQ() for
**  an example.
**----------------------------------------------------------------
*/
{
    int errcode;
    *tstep = 0;
    if (!OpenQflag) return(105);
    errcode = nextqual(tstep);
    if (!errcode && simulation->Saveflag && *tstep == 0) SaveQflag = TRUE;
    if (errcode) errmsg(errcode);
    return(errcode);
}

#endif
#if 0
int DLLEXPORT ENstepQ(long *tleft)
/*----------------------------------------------------------------
**  Input:   none
**  Output:  *tleft = time left in overall simulation (seconds)
**  Returns: error code
**  Purpose: advances WQ simulation by a single WQ time step
**
**  This function is used in a loop with ENrunQ() to run
**  an extended period WQ simulation.
**----------------------------------------------------------------
*/
{
    int errcode;
    *tleft = 0;
    if (!OpenQflag) return(105);
    errcode = stepqual(tleft);
    if (!errcode && simulation->Saveflag && *tleft == 0) SaveQflag = TRUE;
    if (errcode) errmsg(errcode);
    return(errcode);
}

#endif
#if 0
int DLLEXPORT ENcloseQ()
/*----------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Returns: error code
**  Purpose: frees data allocated by WQ solver
**----------------------------------------------------------------
*/
{
    if (!Openflag) return(102);
    closequal();
    OpenQflag = FALSE;
    return(0);
}
#endif

/*
  ----------------------------------------------------------------
  Functions for generating an output report
  ----------------------------------------------------------------
*/

#if 0
int DLLEXPORT ENwriteline(char *line)
/*----------------------------------------------------------------
**  Input:   line = text string
**  Output:  none
**  Returns: error code
**  Purpose: writes line of text to report file
**----------------------------------------------------------------
*/
{
    if (!Openflag) return(102);
    writeline(line);
    return(0);
}
#endif
#if 0
int DLLEXPORT ENreport()
/*----------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Returns: error code
**  Purpose: writes report to report file
**----------------------------------------------------------------
*/
{
    int  errcode;

/* Check if results saved to binary output file */
    if (!SaveQflag) return(106);
    errcode = writereport();
    if (errcode) errmsg(errcode);
    return(errcode);
}
#endif

#if 0
int  DLLEXPORT ENresetreport()
/*----------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Returns: error code
**  Purpose: resets report options to default values
**----------------------------------------------------------------
*/
{
    int i;
    if (!Openflag) return(102);
    initreport();
    for (i=1; i<=Nnodes; i++) Node[i].Rpt = 0;
    for (i=1; i<=Nlinks; i++) Link[i].Rpt = 0;
    return(0);
}
#endif

#if 0
int  DLLEXPORT ENsetreport(char *s)
/*----------------------------------------------------------------
**  Input:   s = report format command
**  Output:  none
**  Returns: error code
**  Purpose: processes a reporting format command
**----------------------------------------------------------------
*/
{
    char s1[MAXLINE+1];
    if (!Openflag) return(102);
    if (strlen(s) > MAXLINE) return(250);
    strcpy(s1,s);
    if (setreport(s1) > 0) return(250);
    else return(0);
}
#endif

/*
  ----------------------------------------------------------------
  Functions for retrieving network information
  ----------------------------------------------------------------
*/

/*** Updated 10/25/00 ***/
int DLLEXPORT ENgetversion(int *v)
/*----------------------------------------------------------------
**  Input:    none
**  Output:   *v = version number of the source code
**  Returns:  error code (should always be 0)
**  Purpose:  retrieves a number assigned to the most recent
**            update of the source code. This number, set by the
**            constant CODEVERSION found in TYPES.H,  began with
**            20001 and increases by 1 with each new update.
**----------------------------------------------------------------
*/
{
    *v = CODEVERSION;
    return(0);
}


int DLLEXPORT ENgetcontrol(ENsimulation_t sim, 
                           int cindex, int *ctype, int *lindex,
                           float *setting, int *nindex, float *level)
/*----------------------------------------------------------------
**  Input:   cindex   = control index (position of control statement
**                      in the input file, starting from 1)
**  Output:  *ctype   = control type code (see TOOLKIT.H)
**           *lindex  = index of controlled link
**           *setting = control setting on link
**           *nindex  = index of controlling node (0 for TIMER
**                      or TIMEOFDAY control)
**           *level   = control level (tank level, junction
**                      pressure, or time (seconds))
**  Returns: error code
**  Purpose: retrieves parameters that define a simple control
**----------------------------------------------------------------
*/
{
    *ctype = 0;
    *lindex = 0;
    *setting = 0.0;
    *nindex = 0;
    *level = 0.0;
    if (!Openflag) return(102);
    if (cindex < 1 || cindex > sim->Ncontrols) return(241);
    *ctype = sim->Control[cindex].Type;
    *lindex = sim->Control[cindex].Link;
    *setting = sim->Control[cindex].Setting;
    if (sim->Control[cindex].Setting != MISSING) 
        switch (Link[*lindex].Type)
    {
    case PRV:
    case PSV:
    case PBV: *setting *= Ucf[PRESSURE]; break;
    case FCV: *setting *= Ucf[FLOW];
    }
    else if (sim->Control[cindex].Status == OPEN) *setting = 1.0;

/*** Updated 3/1/01 ***/
    else *setting = 0.0;

    *nindex = sim->Control[cindex].Node;
    if (*nindex > Nnodes)
        return 0; /* This node is not valid, surely the control is disabled */
    else if (*nindex > Njuncs)
        *level = (sim->Control[cindex].Grade - Node[*nindex].El)*Ucf[ELEV];
    else if (*nindex > 0)
        *level = (sim->Control[cindex].Grade - Node[*nindex].El)*Ucf[PRESSURE];
    else
        *level = (float) sim->Control[cindex].Time;
    return(0);
}


int DLLEXPORT ENgetcount(int code, int *count)
/*----------------------------------------------------------------
**  Input:   code = component code (see TOOLKIT.H)
**  Output:  *count = number of components in network
**  Returns: error code
**  Purpose: retrieves the number of components of a
**           given type in the network
**----------------------------------------------------------------
*/
{
    *count = 0;
    if (!Openflag) return(102);
    switch (code)
    {
    case EN_NODECOUNT:    *count = Nnodes;    break;
    case EN_TANKCOUNT:    *count = Ntanks;    break;
    case EN_LINKCOUNT:    *count = Nlinks;    break;
    case EN_PATCOUNT:     *count = MaxPats;     break;
    case EN_CURVECOUNT:   *count = Ncurves;   break;
    case EN_CONTROLCOUNT: *count = MaxControls; break;
    case EN_PUMPCOUNT:    *count = Npumps;    break;
    case EN_RESERVCOUNT:  *count = Nreservs;  break;
    case EN_JUNCSCOUNT:   *count = Njuncs;  break;
    default: return(251);
    }
    return(0);
}


int  DLLEXPORT ENgetoption(int code, float *value)
/*----------------------------------------------------------------
**  Input:   code = option code (see TOOLKIT.H)
**  Output:  *value = option value
**  Returns: error code
**  Purpose: gets value for an analysis option
**----------------------------------------------------------------
*/
{
    if (!Openflag) return(102);
    *value = 0.0;
    switch (code)
    {
    case EN_TRIALS:     *value = (float)MaxIter;
        break;
    case EN_ACCURACY:   *value = Hacc;
        break;
#if 0
    case EN_TOLERANCE:  *value = Ctol*Ucf[QUALITY];
        break;
#endif
    case EN_EMITEXPON:  if (Qexp > 0.0) *value = 1.0/Qexp;
        break;
    case EN_DEMANDMULT: *value = Dmult;
        break;
    default:            return(251);
    }
    return(0);
}


int DLLEXPORT ENgettimeparam(int code, long *value)
/*----------------------------------------------------------------
**  Input:   code = time parameter code (see TOOLKIT.H)
**  Output:  *value = value of time parameter
**  Returns: error code
**  Purpose: retrieves value of specific time parameter
**----------------------------------------------------------------
*/
{
    *value = 0;
    if (!Openflag) return(102);
    switch (code)
    {
    case EN_DURATION:     *value = Dur;       break;
    case EN_HYDSTEP:      *value = Hstep;     break;
#if 0
    case EN_QUALSTEP:     *value = Qstep;     break;
#endif
    case EN_PATTERNSTEP:  *value = Pstep;     break;
    case EN_PATTERNSTART: *value = Pstart;    break;
    case EN_REPORTSTEP:   *value = Rstep;     break;
    case EN_REPORTSTART:  *value = Rstart;    break;
    case EN_STATISTIC:    *value = Tstatflag; break;
    case EN_PERIODS:      *value = Nperiods;  break;
    case EN_CLOCKSTART:   *value = Tstart;    break;
    default: return(251);
    }
    return(0);
}


int DLLEXPORT ENgetflowunits(int *code)
/*----------------------------------------------------------------
**  Input:   none
**  Output:  *code = code of flow units in use
**                   (see TOOLKIT.H or TYPES.H)
**  Returns: error code
**  Purpose: retrieves flow units code
**----------------------------------------------------------------
*/
{
    *code = -1;
    if (!Openflag) return(102);
    *code = Flowflag;
    return(0);
}


int  DLLEXPORT  ENgetpatternindex(char *id, int *index)
/*----------------------------------------------------------------
**  Input:   id     = time pattern ID
**  Output:  *index = index of time pattern in list of patterns
**  Returns: error code
**  Purpose: retrieves index of time pattern with specific ID
**----------------------------------------------------------------
*/
{
    int i;
    *index = 0;
    if (!Openflag) return(102);
    for (i=1; i <= MaxPats; i++)
    {
        if (strcmp(id, InpPattern[i].ID) == 0)
        {
            *index = i;
            return(0);
        }
    }
    *index = 0;
    return(205);
}


int DLLEXPORT ENgetpatternid(int index, char *id)
/*----------------------------------------------------------------
**  Input:   index = index of time pattern
**  Output:  id    = pattern ID
**  Returns: error code
**  Purpose: retrieves ID of a time pattern with specific index
**
**  NOTE: 'id' must be able to hold MAXID characters
**----------------------------------------------------------------
*/
{
    strcpy(id,"");
    if (!Openflag) return(102);
    if (index < 1 || index > MaxPats) return(205);
    strcpy (id, InpPattern[index].ID);
    return(0);
}


int DLLEXPORT ENgetpatternlen(int index, int *len)
/*----------------------------------------------------------------
**  Input:   index = index of time pattern
**  Output:  *len  = pattern length (number of multipliers)
**  Returns: error code
**  Purpose: retrieves number of multipliers in a time pattern
**----------------------------------------------------------------
*/
{
    if (!Openflag) return(102);
    if (index < 1 || index > MaxPats) return(205);
    *len = InpPattern[index].Length;
    return(0);
}


int DLLEXPORT ENgetpatternvalue(int index, int period, float *value)
/*----------------------------------------------------------------
**  Input:   index  = index of time pattern
**           period = pattern time period
**  Output:  *value = pattern multiplier
**  Returns: error code
**  Purpose: retrieves multiplier for a specific time period
**           and pattern
**----------------------------------------------------------------
*/
{
    if (!Openflag) return(102);
    if (index < 1 || index > MaxPats) return(205);
    if (period < 1 || period > InpPattern[index].Length) return(251);
    *value = InpPattern[index].F[period-1];
    return(0);
}

#if 0
int  DLLEXPORT ENgetqualtype(int *qualcode, int *tracenode)
/*----------------------------------------------------------------
**  Input:   none
**  Output:  *qualcode  = WQ analysis code number (see TOOLKIT.H)
**           *tracenode = index of node being traced (if
**                        qualocode = WQ tracing)
**  Returns: error code
**  Purpose: retrieves type of quality analysis called for
**----------------------------------------------------------------
*/
{
    *tracenode = 0;
    if (!Openflag) return(102);
    *qualcode = Qualflag;
    if (Qualflag == TRACE) *tracenode = TraceNode;
    return(0);
}
#endif

int  DLLEXPORT ENgeterror(int errcode, char *errmsg, int n)
/*----------------------------------------------------------------
**  Input:   errcode = error/warning code number
**           n       = maximum length of string errmsg
**  Output:  errmsg  = text of error/warning message
**  Returns: error code
**  Purpose: retrieves text of error/warning message
**----------------------------------------------------------------
*/
{
    switch (errcode)
    {
    case 1:  strncpy(errmsg,WARN1,n);   break;
    case 2:  strncpy(errmsg,WARN2,n);   break;
    case 3:  strncpy(errmsg,WARN3,n);   break;
    case 4:  strncpy(errmsg,WARN4,n);   break;
    case 5:  strncpy(errmsg,WARN5,n);   break;
    case 6:  strncpy(errmsg,WARN6,n);   break;
    default: geterrmsg (errcode, errmsg, n);
    }
    if (strlen(errmsg) == 0) return(251);
    else return(0);
}


/*
  ----------------------------------------------------------------
  Functions for retrieving node data
  ----------------------------------------------------------------
*/


int DLLEXPORT ENgetnodeindex(char *id, int *index)
/*----------------------------------------------------------------
**  Input:   id = node ID
**  Output:  *index = index of node in list of nodes
**  Returns: error code
**  Purpose: retrieves index of a node with specific ID
**----------------------------------------------------------------
*/
{
    *index = 0;
    if (!Openflag) return(102);
    *index = findnode(id);
    if (*index == 0) return(203);
    else return(0);
}


int DLLEXPORT ENgetnodeid(int index, char *id)
/*----------------------------------------------------------------
**  Input:   index = index of node in list of nodes
**  Output:  id = node ID
**  Returns: error code
**  Purpose: retrieves ID of a node with specific index
**
**  NOTE: 'id' must be able to hold MAXID characters
**----------------------------------------------------------------
*/
{
    strcpy(id,"");
    if (!Openflag) return(102);
    if (index < 1 || index > Nnodes) return(203);
    strcpy(id,Node[index].ID);
    return(0);
}


int  DLLEXPORT ENgetnodetype(int index, int *code)
/*----------------------------------------------------------------
**  Input:   index = node index
**  Output:  *code = node type code number (see TOOLKIT.H)
**  Returns: error code
**  Purpose: retrieves node type of specific node
**----------------------------------------------------------------
*/
{
    *code = -1;
    if (!Openflag) return(102);
    if (index < 1 || index > Nnodes) return(203);
    if (index <= Njuncs) *code = EN_JUNCTION;
    else
    {
        if (Tank[index-Njuncs].A == 0.0) *code = EN_RESERVOIR;
        else *code = EN_TANK;
    }
    return(0);
}


int DLLEXPORT ENgetnodevalue(ENsimulation_t sim, 
                             int index, int code, float *value)
/*----------------------------------------------------------------
**  Input:   index = node index
**           code  = node parameter code (see TOOLKIT.H)
**  Output:  *value = value of node's parameter
**  Returns: error code
**  Purpose: retrieves parameter value for a node
**----------------------------------------------------------------
*/
{
    Pdemand demand;
#if 0
    Psource source;
#endif

/* Check for valid arguments */
    *value = 0.0;
    if (!Openflag) return(102);
    if (index <= 0 || index > Nnodes) return(203);

/* Retrieve called-for parameter */
    switch (code)
    {
    case EN_ELEVATION:
        *value = Node[index].El*Ucf[ELEV];
        break;

    case EN_BASEDEMAND:
        *value = 0.0;
        /* NOTE: primary demand category is last on demand list */
        if (index <= Njuncs)
            for (demand = Node[index].D; demand != NULL; demand = demand->next)
                *value = (demand->Base);
        *value *= Ucf[FLOW];
        break;

    case EN_PATTERN:
        *value = 0.0;
        /* NOTE: primary demand category is last on demand list */
        if (index <= Njuncs)
        {
            for (demand = Node[index].D; demand != NULL; demand = demand->next)
                *value = (float)(demand->Pat);
        }
        else *value = (float)(Tank[index-Njuncs].Pat);
        break;

    case EN_EMITTER:
        *value = 0.0;
        if (Node[index].Ke > 0.0)
            *value = Ucf[FLOW]/pow((Ucf[PRESSURE]*Node[index].Ke),(1.0/Qexp));
        break;
#if 0
    case EN_INITQUAL:
        *value = Node[index].C0*Ucf[QUALITY];
        break;

    case EN_SOURCEQUAL:
        source = Node[index].S;
        if (source == NULL) return(240);
        *value = source->C0;
        break;

    case EN_SOURCETYPE:
        source = Node[index].S;
        if (source == NULL) return(240);
        *value = source->Type;
        break;

    case EN_SOURCEMASS:
        source = Node[index].S;
        if (source == NULL) return(240);
        *value = source->Smass*60.0;
        break;

    case EN_SOURCEPAT:
        source = Node[index].S;
        if (source == NULL) return(240);
        *value = source->Pat;
        break;
#endif
    case EN_TANKLEVEL:
        if (index <= Njuncs) return(251);
        *value = (Tank[index-Njuncs].H0 - Node[index].El)*Ucf[ELEV];
        break;

    case EN_MAXLEVEL:
        if (index <= Njuncs) return(251);
        *value = (Tank[index-Njuncs].Hmax - Node[index].El)*Ucf[ELEV];
        break;

    case EN_MINLEVEL:
        if (index <= Njuncs) return(251);
        *value = (Tank[index-Njuncs].Hmin - Node[index].El)*Ucf[ELEV];
        break;

    case EN_INITVOL:
        if(index <= Njuncs) return(251);
        *value = (Tank[index-Njuncs].V0)*Ucf[VOLUME];
        break;

    case EN_VOLUME:
        if(index <= Njuncs) return(251);
        *value = tankvolume(index-Njuncs, sim->H[index]) * Ucf[VOLUME];
        break;

    case EN_DEMAND:
        *value = sim->D[index] * Ucf[FLOW];
        break;

    case EN_HEAD:
        *value = sim->H[index] * Ucf[HEAD];
        break;

    case EN_PRESSURE:
        *value = (sim->H[index] - Node[index].El)*Ucf[PRESSURE];
        break;
#if 0
    case EN_QUALITY:
        *value = sim->C[index] * Ucf[QUALITY];
        break;
#endif
    default: return(251);
    }
    return(0);
}


/*
  ----------------------------------------------------------------
  Functions for retrieving link data
  ----------------------------------------------------------------
*/


int DLLEXPORT ENgetlinkindex(char *id, int *index)
/*----------------------------------------------------------------
**  Input:   id = link ID
**  Output:  *index = index of link in list of links
**  Returns: error code
**  Purpose: retrieves index of a link with specific ID
**----------------------------------------------------------------
*/
{
    *index = 0;
    if (!Openflag) return(102);
    *index = findlink(id);
    if (*index == 0) return(204);
    else return(0);
}


int DLLEXPORT ENgetlinkid(int index, char *id)
/*----------------------------------------------------------------
**  Input:   index = index of link in list of links
**  Output:  id = link ID
**  Returns: error code
**  Purpose: retrieves ID of a link with specific index
**
**  NOTE: 'id' must be able to hold MAXID characters
**----------------------------------------------------------------
*/
{
    strcpy(id,"");
    if (!Openflag) return(102);
    if (index < 1 || index > Nlinks) return(204);
    strcpy(id,Link[index].ID);
    return(0);
}


int  DLLEXPORT ENgetlinktype(int index, int *code)
/*------------------------------------------------------------------
**  Input:   index = link index
**  Output:  *code = link type code number (see TOOLKIT.H)
**  Returns: error code
**  Purpose: retrieves link type of specific link
**------------------------------------------------------------------
*/
{
    *code = -1;
    if (!Openflag) return(102);
    if (index < 1 || index > Nlinks) return(204);
    *code = Link[index].Type;
    return(0);
}


int  DLLEXPORT ENgetlinknodes(int index, int *node1, int *node2)
/*----------------------------------------------------------------
**  Input:   index = link index
**  Output:  *node1 = index of link's starting node
**           *node2 = index of link's ending node
**  Returns: error code
**  Purpose: retrieves end nodes of a specific link
**----------------------------------------------------------------
*/
{
    *node1 = 0;
    *node2 = 0;
    if (!Openflag) return(102);
    if (index < 1 || index > Nlinks) return(204);
    *node1 = Link[index].N1;
    *node2 = Link[index].N2;
    return(0);
}


int DLLEXPORT ENgetlinkvalue(ENsimulation_t sim, 
                             int index, int code, float *value)
/*------------------------------------------------------------------
**  Input:   index = link index
**           code  = link parameter code (see TOOLKIT.H)
**  Output:  *value = value of link's parameter
**  Returns: error code
**  Purpose: retrieves parameter value for a link
**------------------------------------------------------------------
*/
{
    float a,h,q;
    int pump_idx, i;

/* Check for valid arguments */
    *value = 0.0;
    if (!Openflag) return(102);
    if (index <= 0 || index > Nlinks) return(204);

/* Retrieve called-for parameter */
    switch (code)
    {
    case EN_DIAMETER:
        if (Link[index].Type == PUMP) *value = 0.0;
        else *value = Link[index].Diam*Ucf[DIAM];
        break;

    case EN_LENGTH:
        *value = Link[index].Len*Ucf[ELEV];
        break;

    case EN_ROUGHNESS:
        if (Link[index].Type <= PIPE)
        {
            if (Formflag == DW)
                *value = Link[index].Kc*(1000.0*Ucf[ELEV]);
            else *value = Link[index].Kc;
        }
        else *value = 0.0;
        break;

    case EN_MINORLOSS:
        if (Link[index].Type != PUMP)
        {
            *value = Link[index].Km;
            *value *= (SQR(Link[index].Diam)*SQR(Link[index].Diam)/0.02517);
        }
        else *value = 0.0;
        break;

    case EN_INITSTATUS:
        if (Link[index].Stat <= CLOSED) *value = 0;
        else *value = 1;
        break;

    case EN_INITSETTING:
        if (Link[index].Type == PIPE || Link[index].Type == CV)
            return ENgetlinkvalue (sim, index, EN_ROUGHNESS, value);
        *value = Link[index].Kc;
        switch (Link[index].Type)
        {
        case PRV:
        case PSV:
        case PBV: *value *= Ucf[PRESSURE]; break;
        case FCV: *value *= Ucf[FLOW];
        }
        break;
#if 0
    case EN_KBULK:
        *value = Link[index].Kb*SECperDAY;
        break;

    case EN_KWALL:
        *value = Link[index].Kw*SECperDAY;
        break;
#endif
    case EN_FLOW:

/*** Updated 10/25/00 ***/
        if (sim->S[index] <= CLOSED) *value = 0.0;

        else *value = sim->Q[index] * Ucf[FLOW];
        break;

    case EN_VELOCITY:
        if (Link[index].Type == PUMP) *value = 0.0;

/*** Updated 11/19/01 ***/
        else if (sim->S[index] <= CLOSED) *value = 0.0;

        else
        {
            q = ABS (sim->Q[index]);
            a = PI*SQR(Link[index].Diam)/4.0;
            *value = q/a*Ucf[VELOCITY];
        }
        break;

    case EN_HEADLOSS:

/*** Updated 11/19/01 ***/
        if (sim->S[index] <= CLOSED) *value = 0.0;

        else
        {
            h = sim->H[Link[index].N1] - sim->H[Link[index].N2];
            if (Link[index].Type != PUMP) h = ABS(h);
            *value = h*Ucf[HEADLOSS];
        }
        break;

    case EN_STATUS:
        if (sim->S[index] <= CLOSED) *value = 0;
        else *value = 1;
        break;

    case EN_SETTING:
        if (Link[index].Type == PIPE || Link[index].Type == CV)
            return ENgetlinkvalue (sim, index, EN_ROUGHNESS, value);
        if (sim->K[index] == MISSING) *value = 0.0;
        else                     *value = sim->K[index];
        switch (Link[index].Type)
        {
        case PRV:
        case PSV:
        case PBV: *value *= Ucf[PRESSURE]; break;
        case FCV: *value *= Ucf[FLOW];
        }
        break;

    case EN_SCHEDULE:
        if(Link[index].Type != PUMP) return(204);

        pump_idx = (int) Link[index].Diam;
        for(i=0; i < 24; i++)
            value[i] = sim->Pump[pump_idx].Schedule[i];
        break;


    case EN_ENERGY:
        getenergy (sim, index, value, &a);
        break;

    default: return(251);
    }
    return(0);
}


/*
  ----------------------------------------------------------------
  Functions for changing network data
  ----------------------------------------------------------------
*/


int DLLEXPORT ENsetcontrol(ENsimulation_t sim, 
                           int cindex, int ctype, int lindex,
                           float setting, int nindex, float level)
/*----------------------------------------------------------------
**  Input:   cindex  = control index (position of control statement
**                     in the input file, starting from 1)
**           ctype   = control type code (see TOOLKIT.H)
**           lindex  = index of controlled link
**           setting = control setting applied to link
**           nindex  = index of controlling node (0 for TIMER
**                     or TIMEOFDAY control)
**           level   = control level (tank level, junction pressure,
**                     or time (seconds))
**  Output:  none
**  Returns: error code
**  Purpose: specifies parameters that define a simple control
**----------------------------------------------------------------
*/
{
    char  status = ACTIVE;
    long  t = 0;

    int errcode = 0;

/* Check that input file opened */
    if (!Openflag) return(102);

    /* Check that control exists */
    if (cindex > sim->Ncontrols) {
        sim->Control = realloc(sim->Control, 
                                      (cindex + 1) * sizeof(Scontrol));
        ERRCODE(MEMCHECK(sim->Control));
        if(errcode) return(errcode);
        sim->Ncontrols = cindex;
    }

    if (cindex < 1 || cindex > sim->Ncontrols) return(241);

/* Check that controlled link exists */
    if (lindex == 0)
    {
        sim->Control[cindex].Link = 0;
        return(0);
    }
    if (lindex < 0 || lindex > Nlinks) return(204);

/* Cannot control check valve. */
    if (Link[lindex].Type == CV) return(207);

/* Check for valid parameters */
    if (ctype < 0 || ctype > EN_TIMEOFDAY) return(251);
    if (ctype == EN_LOWLEVEL || ctype == EN_HILEVEL)
    {
        if (nindex < 1 || nindex > Nnodes) return(203);
    }
    else nindex = 0;
    if (setting < 0.0 || level < 0.0) return(202);

/* Adjust units of control parameters */
    switch (Link[lindex].Type)
    {
    case PRV:
    case PSV:
    case PBV:  setting /= Ucf[PRESSURE];
        break;
    case FCV:  setting /= Ucf[FLOW];
        break;

/*** Updated 9/7/00 ***/
    case GPV:  if (setting == 0.0) status = CLOSED;
    else if (setting == 1.0) status = OPEN;
    else return(202);
        setting = Link[lindex].Kc;
        break;

    case PIPE:
    case PUMP: status = OPEN;
        if (setting == 0.0) status = CLOSED;
    }
    if (ctype == LOWLEVEL || ctype == HILEVEL)
    {
        if (nindex > Njuncs) level = Node[nindex].El + level/Ucf[ELEV];
        else level = Node[nindex].El + level/Ucf[PRESSURE];
    }
    if (ctype == TIMER)     t = (long)ROUND(level);
    if (ctype == TIMEOFDAY) t = (long)ROUND(level) % SECperDAY;

/* Reset control's parameters */
    sim->Control[cindex].Type = (char)ctype;
    sim->Control[cindex].Link = lindex;
    sim->Control[cindex].Node = nindex;
    sim->Control[cindex].Status = status;
    sim->Control[cindex].Setting = setting;
    sim->Control[cindex].Grade = level;
    sim->Control[cindex].Time = t;
    return(0);
}

#if 0
int DLLEXPORT ENsetnodevalue(int index, int code, float value)
/*----------------------------------------------------------------
**  Input:   index = node index
**           code  = node parameter code (see TOOLKIT.H)
**           value = parameter value
**  Output:  none
**  Returns: error code
**  Purpose: sets input parameter value for a node
**----------------------------------------------------------------
*/
{
    int  j;
    Pdemand demand;
#if 0
    Psource source;
#endif

    if (!Openflag) return(102);
    if (index <= 0 || index > Nnodes) return(203);
    switch (code)
    {
    case EN_ELEVATION:
        if (index <= Njuncs) Node[index].El = value/Ucf[ELEV];
        else
        {
            value = (value/Ucf[ELEV]) - Node[index].El;
            j = index - Njuncs;
            Tank[j].H0 += value;
            Tank[j].Hmin += value;
            Tank[j].Hmax += value;
            Node[index].El += value;
            simulation->H[index] += value;
        }
        break;

    case EN_BASEDEMAND:
        /* NOTE: primary demand category is last on demand list */
        if (index <= Njuncs)
        {
            for (demand = Node[index].D; demand != NULL; demand = demand ->next)
            {
                if (demand->next == NULL) demand->Base = value/Ucf[FLOW];
            }
        }
        break;

    case EN_PATTERN:
        /* NOTE: primary demand category is last on demand list */
        j = ROUND(value);
        if (j < 0 || j > MaxPats) return(205);
        if (index <= Njuncs)
        {
            for (demand = Node[index].D; demand != NULL; demand = demand ->next)
            {
                if (demand->next == NULL) demand->Pat = j;
            }
        }
        else Tank[index-Njuncs].Pat = j;
        break;

    case EN_EMITTER:
        if (index > Njuncs) return(203);
        if (value < 0.0) return(202);
        if (value > 0.0)
            value = pow((Ucf[FLOW]/value),Qexp)/Ucf[PRESSURE];
        Node[index].Ke = value;
        break;
#if 0
    case EN_INITQUAL:
        if (value < 0.0) return(202);
        Node[index].C0 = value/Ucf[QUALITY];
        if (index > Njuncs) Tank[index-Njuncs].C = Node[index].C0;
        break;

    case EN_SOURCEQUAL:
    case EN_SOURCETYPE:
    case EN_SOURCEPAT:
        if (value < 0.0) return(202);
        source = Node[index].S;
        if (source == NULL)
        {
            source = (struct Ssource *) malloc(sizeof(struct Ssource));
            if (source == NULL) return(101);
            source->Type = CONCEN;
            source->C0 = 0.0;
            source->Pat = 0;
            Node[index].S = source;
        }
        if (code == EN_SOURCEQUAL) source->C0 = value;
        else if (code == EN_SOURCEPAT)
        {
            j = ROUND(value);
            if (j < 0 || j > Npats) return(205);
            source->Pat = j;
        }
        else
        {
            j = ROUND(value);
            if ( j < CONCEN || j > FLOWPACED) return(251);
            else source->Type = (char)j;
        }
        return(0);
#endif
    case EN_TANKLEVEL:
        if (index <= Njuncs) return(251);
        j = index - Njuncs;
        if (Tank[j].A == 0.0)  /* Tank is a reservoir */
        {
            Tank[j].H0 = value/Ucf[ELEV];
            Tank[j].Hmin = Tank[j].H0;
            Tank[j].Hmax = Tank[j].H0;
            Node[index].El = Tank[j].H0;
            simulation->H[index] = Tank[j].H0;
        }
        else
        {
            value = Node[index].El + value/Ucf[ELEV];
            if (value > Tank[j].Hmax
                ||  value < Tank[j].Hmin) return(202);
            Tank[j].H0 = value;
            Tank[j].V0 = tankvolume(j, Tank[j].H0);
            simulation->H[index] = Tank[j].H0;
        }
        break;

    default: return(251);
    }
    return(0);
}
#endif

int DLLEXPORT ENsetlinkvalue(ENsimulation_t sim, 
                             int index, int code, float value)
/*----------------------------------------------------------------
**  Input:   index = link index
**           code  = link parameter code (see TOOLKIT.H)
**           value = parameter value
**  Output:  none
**  Returns: error code
**  Purpose: sets input parameter value for a link
**----------------------------------------------------------------
*/
{
    char  s;
    float r;

    if (!Openflag) return(102);
    if (index <= 0 || index > Nlinks) return(204);
    switch (code)
    {
    case EN_DIAMETER:
        if (Link[index].Type != PUMP)
        {
            if (value <= 0.0) return(202);
            value /= Ucf[DIAM];              /* Convert to feet */
            r = Link[index].Diam/value;      /* Ratio of old to new diam */
            Link[index].Km *= SQR(r)*SQR(r); /* Adjust minor loss factor */
            Link[index].Diam = value;        /* Update diameter */
            resistance(index);               /* Update resistance factor */
        }
        break;

    case EN_LENGTH:
        if (Link[index].Type <= PIPE)
        {
            if (value <= 0.0) return(202);
            Link[index].Len = value/Ucf[ELEV];
            resistance(index);
        }
        break;

    case EN_ROUGHNESS:
        if (Link[index].Type <= PIPE)
        {
            if (value <= 0.0) return(202);
            Link[index].Kc = value;
            if (Formflag  == DW) Link[index].Kc /= (1000.0*Ucf[ELEV]);
            resistance(index);
        }
        break;

    case EN_MINORLOSS:
        if (Link[index].Type != PUMP)
        {
            if (value <= 0.0) return(202);
            Link[index].Km = 0.02517*value/SQR(Link[index].Diam)/SQR(Link[index].Diam);
        }
        break;

    case EN_INITSTATUS:
    case EN_STATUS:
        /* Cannot set status for a check valve */
        if (Link[index].Type == CV) return(207);
        s = (char)ROUND(value);
        if (s < 0 || s > 1) return(251);
        if (code == EN_INITSTATUS)
            setlinkstatus (sim, index, s, 
                           &Link[index].Stat, &Link[index].Kc);
        else
            setlinkstatus (sim, index, s, 
                           &(sim->S[index]), &(sim->K[index]));
        break;

    case EN_INITSETTING:
    case EN_SETTING:
        if (value < 0.0) return(202);
        if (Link[index].Type == PIPE || Link[index].Type == CV)
            return ENsetlinkvalue(sim, index, EN_ROUGHNESS, value);
        else
        {
            switch (Link[index].Type)
            {
            case PUMP: break;
            case PRV:
            case PSV:
            case PBV: value /= Ucf[PRESSURE]; break;
            case FCV: value /= Ucf[FLOW]; break;
            case TCV: break;

/***  Updated 9/7/00  ***/
            case GPV: return(202);  /* Cannot modify setting for GPV */

            default:  return(251);
            }
            if (code == EN_INITSETTING)
                setlinksetting (sim, index, value, &Link[index].Stat, &Link[index].Kc);
            else
                setlinksetting (sim, index, value, &(sim->S[index]), &(sim->K[index]));
        }
        break;
#if 0
    case EN_KBULK:
        if (Link[index].Type <= PIPE) Link[index].Kb = value/SECperDAY;
        break;

    case EN_KWALL:
        if (Link[index].Type <= PIPE) Link[index].Kw = value/SECperDAY;
        break;
#endif
    case EN_UPATTERN:
        if(Link[index].Type != PUMP) return(204);
        Pump[(int) Link[index].Diam].Upat = (int) value;
        break;

    default: return(251);
    }
    return(0);
}

/*
**--------------------------------------------------------------------
**  Input:   id    = pattern ID label
**  Output:  returns error code
**  Purpose: adds a new pattern to the database
**--------------------------------------------------------------------
*/
int  DLLEXPORT  ENaddpattern(char *id)
{
    int errcode;
    int index = 0;

    errcode = ENgetpatternindex(id, &index);

    if(index == 0) {
        errcode = 0;

        InpPattern = realloc(InpPattern, (MaxPats+2)*sizeof(Spattern));
        ERRCODE(MEMCHECK(InpPattern));

        if(errcode) return(errcode);

        (MaxPats)++;
        strncpy (InpPattern[MaxPats].ID, id, MAXID);
        InpPattern[MaxPats].Length = 0;
        InpPattern[MaxPats].F = NULL;
    }
    return(errcode);
}

int  DLLEXPORT  ENsetpattern(ENsimulation_t sim, int index, float *f, int n)
/*----------------------------------------------------------------
**   Input:   index = time pattern index
**            *f    = array of pattern multipliers
**            n     = number of time periods in pattern
**   Output:  none
**   Returns: error code
**   Purpose: sets multipliers for a specific time pattern
**----------------------------------------------------------------
*/
{
    int j;

/* Check for valid arguments */
    if (!Openflag) return(102);
    if (index <= 0 || index > sim->Npats) return(205);
    if (n <= 0) return(202);

/* Re-set number of time periods & reallocate memory for multipliers */
    sim->Pattern[index].Length = n;
    sim->Pattern[index].F = realloc (sim->Pattern[index].F, n * sizeof(float));
    if (sim->Pattern[index].F == NULL) return(101);

/* Load multipliers into pattern */
    for (j=0; j<n; j++) sim->Pattern[index].F[j] = f[j];
    return(0);
}


int  DLLEXPORT  ENsetpatternvalue(ENsimulation_t sim, 
                                  int index, int period, float value)
/*----------------------------------------------------------------
**  Input:   index  = time pattern index
**           period = time pattern period
**           value  = pattern multiplier
**  Output:  none
**  Returns: error code
**  Purpose: sets multiplier for a specific time period and pattern
**----------------------------------------------------------------
*/
{
    if (!Openflag) return(102);
    if (index  <= 0 || index  > sim->Npats) return(205);
    if (period <= 0 || period > sim->Pattern[index].Length) return(251);
    sim->Pattern[index].F[period-1] = value;
    return(0);
}


int  DLLEXPORT  ENsettimeparam(int code, long value)
/*----------------------------------------------------------------
**  Input:   code  = time parameter code (see TOOLKIT.H)
**           value = time parameter value
**  Output:  none
**  Returns: error code
**  Purpose: sets value for time parameter
**----------------------------------------------------------------
*/
{
    if (!Openflag) return(102);
    /* FIXME: not sure about the proper behaviour here. If you call
       ENsettimeparam in the middle of a simulation, anything may
       happen.  */
    /*    if (OpenHflag || OpenQflag) return(109);*/
    if (value < 0) return(202);
    switch(code)
    {
    case EN_DURATION:      Dur = value;
        if (Rstart > Dur) Rstart = 0;
        break;
    case EN_HYDSTEP:       if (value == 0) return(202);
        Hstep = value;
        Hstep = MIN(Pstep, Hstep);
        Hstep = MIN(Rstep, Hstep);
#if 0
        Qstep = MIN(Qstep, Hstep);
#endif
        break;
#if 0
    case EN_QUALSTEP:      if (value == 0) return(202);
        Qstep = value;
        Qstep = MIN(Qstep, Hstep);
        break;
#endif
    case EN_PATTERNSTEP:   if (value == 0) return(202);
        Pstep = value;
        if (Hstep > Pstep) Hstep = Pstep;
        break;
    case EN_PATTERNSTART:  Pstart = value;
        break;
    case EN_REPORTSTEP:    if (value == 0) return(202);
        Rstep = value;
        if (Hstep > Rstep) Hstep = Rstep;
        break;
    case EN_REPORTSTART:   if (Rstart > Dur) return(202);
        Rstart = value;
        break;
    case EN_RULESTEP:      if (value == 0) return(202);
        Rulestep = value;
        Rulestep = MIN(Rulestep, Hstep);
        break;
    case EN_STATISTIC:     if (value > RANGE) return(202);
        Tstatflag = (char)value;
        break;
    default:               return(251);
    }
    return(0);
}


int  DLLEXPORT ENsetoption(int code, float value)
/*----------------------------------------------------------------
**  Input:   code  = option code (see TOOLKIT.H)
**           value = option value
**  Output:  none
**  Returns: error code
**  Purpose: sets value for an analysis option
**----------------------------------------------------------------
*/
{
    if (!Openflag) return(102);
    switch (code)
    {
    case EN_TRIALS:     if (value < 1.0) return(202);
        MaxIter = (int)value;
        break;
    case EN_ACCURACY:   if (value < 1.e-5 || value > 1.e-1) return(202);
        Hacc = value;
        break;
#if 0
    case EN_TOLERANCE:  if (value < 0.0) return(202);
        Ctol = value/Ucf[QUALITY];
        break;
#endif
#if 0
    case EN_EMITEXPON:
    {   
        int   i,j;
        float Ke,n,ucf;
        if (value <= 0.0) return(202);
        n = 1.0/value;
        ucf = pow(Ucf[FLOW],n)/Ucf[PRESSURE];
        for (i=1; i<=Njuncs; i++)
        {
            j = ENgetnodevalue(i,EN_EMITTER,&Ke);
            if (j == 0 && Ke > 0.0) Node[i].Ke = ucf/pow(Ke,n);
        }
        Qexp = n;
        break;
    }
#endif
    case EN_DEMANDMULT: if (value <= 0.0) return(202);
        Dmult = value;
        break;
    default:            return(251);
    }
    return(0);
}

#if 0
int  DLLEXPORT ENsetstatusreport(int code)
/*----------------------------------------------------------------
**  Input:   code = status reporting code (0, 1, or 2)
**  Output:  none
**  Returns: error code
**  Purpose: sets level of hydraulic status reporting
**----------------------------------------------------------------
*/
{
    int errcode = 0;
    if (code >= 0 && code <= 2) Statflag = (char)code;
    else errcode = 202;
    return(errcode);
}
#endif
#if 0
int  DLLEXPORT ENsetqualtype(int qualcode, char *chemname,
                             char *chemunits, char *tracenode)
/*----------------------------------------------------------------
**  Input:   qualcode  = WQ parameter code (see TOOLKIT.H)
**           chemname  = name of WQ constituent
**           chemunits = concentration units of WQ constituent
**           tracenode = ID of node being traced
**  Output:  none
**  Returns: error code
**  Purpose: sets type of quality analysis called for
**
**  NOTE: chemname and chemunits only apply when WQ analysis
**        is for chemical. tracenode only applies when WQ
**        analysis is source tracing.
**----------------------------------------------------------------
*/
{
/*** Updated 3/1/01 ***/
    float ccf = 1.0;

    if (!Openflag) return(102);
    if (qualcode < EN_NONE || qualcode > EN_TRACE) return(251);
    Qualflag = (char)qualcode;
    if (Qualflag == CHEM)                   /* Chemical constituent */
    {
        strncpy(ChemName,chemname,MAXID);
        strncpy(ChemUnits,chemunits,MAXID);

/*** Updated 3/1/01 ***/
        strncpy(Field[QUALITY].Units,ChemUnits,MAXID);
        strncpy(Field[REACTRATE].Units,ChemUnits,MAXID);
        strcat(Field[REACTRATE].Units,t_PERDAY);
        ccf =  1.0/LperFT3;

    }
    if (Qualflag == TRACE)                  /* Source tracing option */
    {
        TraceNode = findnode(tracenode);
        if (TraceNode == 0) return(203);
        strncpy(ChemName,u_PERCENT,MAXID);
        strncpy(ChemUnits,tracenode,MAXID);

/*** Updated 3/1/01 ***/
        strcpy(Field[QUALITY].Units,u_PERCENT);
    }
    if (Qualflag == AGE)                    /* Water age analysis */
    {
        strncpy(ChemName,w_AGE,MAXID);
        strncpy(ChemUnits,u_HOURS,MAXID);

/*** Updated 3/1/01 ***/
        strcpy(Field[QUALITY].Units,u_HOURS);
    }

/*** Updated 3/1/01 ***/
    Ucf[QUALITY]   = ccf;
    Ucf[LINKQUAL]  = ccf;
    Ucf[REACTRATE] = ccf;

    return(0);
}
#endif

/*
  ----------------------------------------------------------------
  Functions for opening files
  ----------------------------------------------------------------
*/


int   openfiles(char *f1, char *f2, char *f3)
/*----------------------------------------------------------------
**  Input:   f1 = pointer to name of input file
**           f2 = pointer to name of report file
**           f3 = pointer to name of binary output file
**  Output:  none
**  Returns: error code
**  Purpose: opens input & report files
**----------------------------------------------------------------
*/
{
/* Initialize file pointers to NULL */
    InFile = NULL;
#if 0
    RptFile = NULL;
    OutFile = NULL;
    HydFile = NULL;
#endif

/* Save file names */
    strncpy(InpFname,f1,MAXFNAME);
    strncpy(Rpt1Fname,f2,MAXFNAME);
    strncpy(OutFname,f3,MAXFNAME);

    /* Check that file names are not identical, unless both report and
       binary output files are the empty string "" */
    if (strcomp(f1,f2) || strcomp(f1,f3) ||
        (strlen(f2) != 0 && strlen(f3) != 0 && strcomp(f2,f3)) ) {
        writecon(FMT04);
        return(301);
    }

/* Attempt to open input and report files */
    if ((InFile = fopen(f1,"rt")) == NULL)
    {
        writecon(FMT05);
        writecon(f1);
        return(302);
    }
#if 0
    if (strlen(f2) == 0) RptFile = stdout;
    else if ((RptFile = fopen(f2,"wt")) == NULL)
    {
        writecon(FMT06);
        return(303);
    }
#endif
    return(0);
}                                       /* End of openfiles */

#if 0
static int  openhydfile()
/*----------------------------------------------------------------
** Input:   none
** Output:  none
** Returns: error code
** Purpose: opens file that saves hydraulics solution
**----------------------------------------------------------------
*/
{
    INT4 nsize[6];                       /* Temporary array */
    INT4 magic;
    INT4 version;
    int errcode = 0;

/* If HydFile currently open, then close it if its not a scratch file */
    if (HydFile != NULL)
    {
        if (Hydflag == SCRATCH) return(0);
        fclose(HydFile);
    }

/* Use Hydflag to determine the type of hydraulics file to use. */
/* Write error message if the file cannot be opened.            */
    HydFile = NULL;
    switch(Hydflag)
    {
    case SCRATCH:  HydFile = tmpfile();
        break;
    case SAVE:     HydFile = fopen(HydFname,"w+b");
        break;
    case USE:      HydFile = fopen(HydFname,"rb");
        break;
    }
    if (HydFile == NULL) return(305);

/* If a previous hydraulics solution is not being used, then */
/* save the current network size parameters to the file.     */
    if (Hydflag != USE)
    {
        magic = MAGICNUMBER;
        version = VERSION;
        nsize[0] = Nnodes;
        nsize[1] = Nlinks;
        nsize[2] = Ntanks;
        nsize[3] = Npumps;
        nsize[4] = Nvalves;
        nsize[5] = Dur;
        fwrite(&magic,sizeof(INT4),1,HydFile);
        fwrite(&version,sizeof(INT4),1,HydFile);
        fwrite(nsize,sizeof(INT4),6,HydFile);
    }

/* If a previous hydraulics solution is being used, then */
/* make sure its network size parameters match those of  */
/* the current network.                                  */
    if (Hydflag == USE)
    {
        fread(&magic,sizeof(INT4),1,HydFile);
        if (magic != MAGICNUMBER) return(306);
        fread(&version,sizeof(INT4),1,HydFile);
        if (version != VERSION) return(306);
        if (fread(nsize,sizeof(INT4),6,HydFile) < 6) return(306);
        if (nsize[0] != Nnodes  || nsize[1] != Nlinks ||
            nsize[2] != Ntanks  || nsize[3] != Npumps ||
            nsize[4] != Nvalves || nsize[5] != Dur) return(306);
        simulation->SaveHflag = TRUE;
    }

/* Save current position in hydraulics file  */
/* where storage of hydraulic results begins */
    HydOffset = ftell(HydFile);
    return(errcode);
}
#endif

#if 0
static int  openoutfile()
/*----------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Returns: error code
**  Purpose: opens binary output file.
**----------------------------------------------------------------
*/
{
    int errcode = 0;

/* Close output file if already opened */
    if (OutFile != NULL) fclose(OutFile);
    OutFile = NULL;
    if (TmpOutFile != NULL) fclose(TmpOutFile);
    TmpOutFile = NULL;

/* If output file name was supplied, then attempt to */
/* open it. Otherwise open a temporary output file.  */
    if (strlen(OutFname) != 0)
    {
        if ( (OutFile = fopen(OutFname,"w+b")) == NULL)
        {
            writecon(FMT07);
            errcode = 304;
        }
    }
    else if ( (OutFile = tmpfile()) == NULL)
    {
        writecon(FMT08);
        errcode = 304;
    }

/* Save basic network data & energy usage results */
    ERRCODE(savenetdata());
    OutOffset1 = ftell(OutFile);
    ERRCODE(saveenergy());
    OutOffset2 = ftell(OutFile);

/* Open temporary file if computing time series statistic */
    if (!errcode)
    {
        if (Tstatflag != SERIES)
        {
            if ( (TmpOutFile = tmpfile()) == NULL) errcode = 304;
        }
        else TmpOutFile = OutFile;
    }
    return(errcode);
}
#endif

/*
  ----------------------------------------------------------------
  Global memory management functions
  ----------------------------------------------------------------
*/


void initpointers()
/*----------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Purpose: initializes global pointers to NULL
**----------------------------------------------------------------
*/
{
    Node     = NULL;
    Link     = NULL;
    Tank     = NULL;
    Pump     = NULL;
    Valve    = NULL;
    Curve    = NULL;

    Patlist  = NULL;
    Curvelist = NULL;
    Nht      = NULL;
    Lht      = NULL;
    initrules();
}


int  allocdata()
/*----------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Returns: error code
**  Purpose: allocates memory for network data structures
**----------------------------------------------------------------
*/
{
    int n;
    int errcode = 0;

/* Allocate node & link ID hash tables */
    Nht = HTcreate();
    Lht = HTcreate();
    ERRCODE(MEMCHECK(Nht));
    ERRCODE(MEMCHECK(Lht));

/* Allocate memory for network nodes */
/*************************************************************
 NOTE: Because network components of a given type are indexed
       starting from 1, their arrays must be sized 1
       element larger than the number of components.
*************************************************************/
    if (!errcode)
    {
        n = MaxNodes + 1;
        Node = (Snode *) calloc(n, sizeof(Snode));
        ERRCODE(MEMCHECK(Node));
    }

/* Allocate memory for network links */
    if (!errcode)
    {
        n = MaxLinks + 1;
        Link = (Slink *) calloc(n, sizeof(Slink));
        ERRCODE(MEMCHECK(Link));
    }

/* Allocate memory for tanks, sources, pumps, valves,   */
/* controls, demands, time patterns, & operating curves */
    if (!errcode)
    {
        Tank    = (Stank *)    calloc(MaxTanks+1,   sizeof(Stank));
        Pump    = (Spump *)    calloc(MaxPumps+1,   sizeof(Spump));
        Valve   = (Svalve *)   calloc(MaxValves+1,  sizeof(Svalve));
        InpControl = calloc (MaxControls + 1, sizeof(Scontrol));
        InpPattern = calloc (MaxPats+1,       sizeof(Spattern));
        Curve   = (Scurve *)   calloc(MaxCurves+1,  sizeof(Scurve));
        ERRCODE(MEMCHECK(Tank));
        ERRCODE(MEMCHECK(Pump));
        ERRCODE(MEMCHECK(Valve));
        ERRCODE(MEMCHECK(InpControl));
        ERRCODE(MEMCHECK(InpPattern));
        ERRCODE(MEMCHECK(Curve));
    }

/* Initialize pointers used in patterns, curves, and demand category lists */
    if (!errcode)
    {
        for (n=0; n<=MaxPats; n++)
        {
            InpPattern[n].Length = 0;
            InpPattern[n].F = NULL;
        }
        for (n=0; n<=MaxCurves; n++)
        {
            Curve[n].Npts = 0;
            Curve[n].Type = -1;
            Curve[n].X = NULL;
            Curve[n].Y = NULL;
        }
        for (n=0; n<=MaxNodes; n++) Node[n].D = NULL;
    }

/* Allocate memory for rule base (see RULES.C) */
    if (!errcode) errcode = allocrules();
    return(errcode);
}                                       /* End of allocdata */


void  freeTmplist(STmplist *t)
/*----------------------------------------------------------------
**  Input:   t = pointer to start of a temporary list
**  Output:  none
**  Purpose: frees memory used for temporary storage
**           of pattern & curve data
**----------------------------------------------------------------
*/
{
    STmplist   *tnext;
    while (t != NULL)
    {
        tnext = t->next;
        freeFloatlist(t->x);
        freeFloatlist(t->y);
        free(t);
        t = tnext;
    }
}


void  freeFloatlist(SFloatlist *f)
/*----------------------------------------------------------------
**  Input:   f = pointer to start of list of floats
**  Output:  none
**  Purpose: frees memory used for storing list of floats
**----------------------------------------------------------------
*/
{
    SFloatlist *fnext;
    while (f != NULL)
    {
        fnext = f->next;
        free(f);
        f = fnext;
    }
}


void  freedata()
/*----------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Purpose: frees memory allocated for network data structures.
**----------------------------------------------------------------
*/
{
    int j;
    Pdemand demand, nextdemand;
#if 0
    Psource source;
#endif

/* Free memory for node data */
    if (Node != NULL)
    {
        for (j=0; j<=MaxNodes; j++)
        {
            /* Free memory used for demand category list */
            demand = Node[j].D;
            while (demand != NULL)
            {
                nextdemand = demand->next;
                free(demand);
                demand = nextdemand;
            }
#if 0
            /* Free memory used for WQ source data */
            source = Node[j].S;
            if (source != NULL) free(source);
#endif
        }
        free(Node);
    }

/* Free memory for other network objects */
    free(Link);
    free(Tank);
    free(Pump);
    free(Valve);
    free(InpControl);

/* Free memory for time patterns */
    if (InpPattern != NULL)
    {
        for (j=0; j<=MaxPats; j++) free(InpPattern[j].F);
        free(InpPattern);
    }

/* Free memory for curves */
    if (Curve != NULL)
    {
        for (j=0; j<=MaxCurves; j++)
        {
            free(Curve[j].X);
            free(Curve[j].Y);
        }
        free(Curve);
    }

/* Free memory for rule base (see RULES.C) */
    freerules();

/* Free hash table memory */
    if (Nht != NULL) HTfree(Nht);
    if (Lht != NULL) HTfree(Lht);
}


/*
  ----------------------------------------------------------------
  General purpose functions
  ----------------------------------------------------------------
*/


int  strcomp(char *s1, char *s2)
/*---------------------------------------------------------------
**  Input:   s1 = character string
**           s2 = character string
**  Output:  none
**  Returns: 1 if s1 is same as s2, 0 otherwise
**  Purpose: case insensitive comparison of strings s1 & s2
**---------------------------------------------------------------
*/
{
    int i;
    for (i=0; UCHAR(s1[i]) == UCHAR(s2[i]); i++)
        if (!s1[i+1] && !s2[i+1]) return(1);
    return(0);
}                                       /*  End of strcomp  */


float  interp(int n, float x[], float y[], float xx)
/*----------------------------------------------------------------
**  Input:   n  = number of data pairs defining a curve
**           x  = x-data values of curve
**           y  = y-data values of curve
**           xx = specified x-value
**  Output:  none
**  Returns: y-value on curve at x = xx
**  Purpose: uses linear interpolation to find y-value on a
**           data curve corresponding to specified x-value.
**  NOTE:    does not extrapolate beyond endpoints of curve.
**----------------------------------------------------------------
*/
{
    int    k,m;
    float  dx,dy;

    m = n - 1;                          /* Highest data index      */
    if (xx <= x[0]) return(y[0]);       /* xx off low end of curve */
    for (k=1; k<=m; k++)                /* Bracket xx on curve     */
    {
        if (x[k] >= xx)                 /* Interp. over interval   */
        {
            dx = x[k]-x[k-1];
            dy = y[k]-y[k-1];
            if (ABS(dx) < TINY) return(y[k]);
            else return(y[k] - (x[k]-xx)*dy/dx);
        }
    }
    return(y[m]);                       /* xx off high end of curve */
}                       /* End of interp */


int   findnode(char *id)
/*----------------------------------------------------------------
**  Input:   id = node ID
**  Output:  none
**  Returns: index of node with given ID, or 0 if ID not found
**  Purpose: uses hash table to find index of node with given ID
**----------------------------------------------------------------
*/
{
    return(HTfind(Nht,id));
}


int  findlink(char *id)
/*----------------------------------------------------------------
**  Input:   id = link ID
**  Output:  none
**  Returns: index of link with given ID, or 0 if ID not found
**  Purpose: uses hash table to find index of link with given ID
**----------------------------------------------------------------
*/
{
    return(HTfind(Lht,id));
}


static char *geterrmsg(int errcode, char *buf, size_t buflen)
/*----------------------------------------------------------------
**  Input:   errcode = error code
**  Output:  none
**  Returns: pointer to string with error message
**  Purpose: retrieves text of error message
**----------------------------------------------------------------
*/
{
    switch (errcode)
    {                                   /* Warnings */
/*
  case 1:     strcpy(buf, WARN1, buflen);   break;
  case 2:     strcpy(buf, WARN2, buflen);   break;
  case 3:     strcpy(buf, WARN3, buflen);   break;
  case 4:     strcpy(buf, WARN4, buflen);   break;
  case 5:     strcpy(buf, WARN5, buflen);   break;
  case 6:     strcpy(buf, WARN6, buflen);   break;
*/
        /* System Errors */
    case 101:   strncpy (buf, ERR101, buflen);  break;
    case 102:   strncpy (buf, ERR102, buflen);  break;
    case 103:   strncpy (buf, ERR103, buflen);  break;
    case 104:   strncpy (buf, ERR104, buflen);  break;
    case 105:   strncpy (buf, ERR105, buflen);  break;
    case 106:   strncpy (buf, ERR106, buflen);  break;
    case 107:   strncpy (buf, ERR107, buflen);  break;
    case 108:   strncpy (buf, ERR108, buflen);  break;
    case 109:   strncpy (buf, ERR109, buflen);  break;
    case 110:   strncpy (buf, ERR110, buflen);  break;
    case 120:   strncpy (buf, ERR120, buflen);  break;

        /* Input Errors */
    case 200:  strncpy (buf, ERR200, buflen);   break;
    case 223:  strncpy (buf, ERR223, buflen);   break;
    case 224:  strncpy (buf, ERR224, buflen);   break;

        /* Toolkit function errors */
    case 202:  snprintf(buf, buflen, ERR202,t_FUNCCALL,""); break;
    case 203:  snprintf(buf, buflen, ERR203,t_FUNCCALL,""); break;
    case 204:  snprintf(buf, buflen, ERR204,t_FUNCCALL,""); break;
    case 205:  snprintf(buf, buflen, ERR205,t_FUNCCALL,""); break;
    case 207:  snprintf(buf, buflen, ERR207,t_FUNCCALL,""); break;
    case 240:  snprintf(buf, buflen, ERR240,t_FUNCCALL,""); break;
    case 241:  snprintf(buf, buflen, ERR241,t_FUNCCALL,""); break;
    case 250:  snprintf(buf, buflen, ERR250);  break;
    case 251:  snprintf(buf, buflen, ERR251);  break;

        /* File Errors */
    case 301:  strncpy (buf, ERR301, buflen);   break;
    case 302:  strncpy (buf, ERR302, buflen);   break;
    case 303:  strncpy (buf, ERR303, buflen);   break;
    case 304:  strncpy (buf, ERR304, buflen);   break;
    case 305:  strncpy (buf, ERR305, buflen);   break;
    case 306:  strncpy (buf, ERR306, buflen);   break;
    case 307:  strncpy (buf, ERR307, buflen);   break;
    case 308:  strncpy (buf, ERR308, buflen);   break;
    case 309:  strncpy (buf, ERR309, buflen);   break;
    default:   strncpy (buf, "", buflen);
    }
    return buf;
}

#if 0
void  errmsg(int errcode)
/*----------------------------------------------------------------
**  Input:   errcode = error code
**  Output:  none
**  Purpose: writes error message to report file
**----------------------------------------------------------------
*/
{
    if (errcode == 309)    /* Report file write error -  */
    {                      /* Do not write msg to file.  */
        writecon("\n  ");
        writecon(geterrmsg(errcode));
    }
    else if (RptFile != NULL && Messageflag)
    {
        writeline(geterrmsg(errcode));
    }
}
#endif

void  writecon(char *s)
/*----------------------------------------------------------------
**  Input:   text string
**  Output:  none
**  Purpose: writes string of characters to console
**----------------------------------------------------------------
*/
{
#ifndef DLL
    fprintf(stdout,s);
    fflush(stdout);
#endif
}


#ifndef LINUX
void writewin(char *s)
/*----------------------------------------------------------------
**  Input:   text string
**  Output:  none
**  Purpose: passes character string to viewprog() in
**           application which calls the EPANET DLL
**----------------------------------------------------------------
*/
{
    char progmsg[MAXMSG+1];
#ifdef DLL
    if (viewprog != NULL)
    {
        strncpy(progmsg,s,MAXMSG);
        viewprog(progmsg);
    }
#endif
}
#endif
/*----------------------------------------------------------------
**  Input:   none
**  Output:  total energy cost per pump plus demand cost
**----------------------------------------------------------------
*/
int  DLLEXPORT ENgettotalenergycost(ENsimulation_t sim, float *cost)
{
    int   j;

    if (!Openflag) return(102);
    if (sim == NULL) return(103);

    if (Npumps == 0) {
        *cost = 0.0;
        return(0);
    }

    *cost = 0.0;
    for (j=1; j<=Npumps; j++)
        *cost += sim->Pump[j].Energy[5];

    *cost +=  (sim->Emax * Dcost);

    return(0);
}
/*----------------------------------------------------------------
**  Input:   link index of a pump
**  Output:  number of switches (OFF->ON) of a pump
**----------------------------------------------------------------
*/
int DLLEXPORT ENgetpumpswitches(ENsimulation_t sim, int index, int *value)
{
    if (!Openflag) return(102);
    if (sim == NULL) return(103);
    if (index <= 0 || index > Nlinks) return(204);
    *value = sim->Pump[PUMPINDEX(index)].Nswitches;
    return(0);
}

/*----------------------------------------------------------------
**  Input:   pump index of a pump [1, Npumps]
**  Output:  link index of a pump
**----------------------------------------------------------------
*/
int DLLEXPORT ENgetpumpindex(int pump_index, int *link_index)
{
    if (!Openflag) return(102);
    if((pump_index <= 0)||(pump_index > Npumps)) return(401);
    *link_index = Pump[pump_index].Link;
    return(0);
}

/*----------------------------------------------------------------
**  Input:   tank index of a tank (not reservoirs) [1, Ntanks]
**  Output:  node index of a tank (not reservoirs)
**----------------------------------------------------------------
*/
int DLLEXPORT ENgettankindex(int tank_index, int *node_index)
{
    int typecode;
    int temp_index;
    int k;

    if (!Openflag) return(102);
    if((tank_index <= 0)||(tank_index > Ntanks)) return(402);

    temp_index = 0;
    k = 1;
    while(temp_index < tank_index) {
        *node_index = Tank[k++].Node;
        ENgetnodetype(*node_index, &typecode);
        if(typecode == EN_TANK) {
            temp_index++;
        }
    }
    return(0);
}
#if 0
/*----------------------------------------------------------------
** Adds rule with format:

 IF SYSTEM CLOCKTIME >= start_time (in seconds)
 AND SYSTEM CLOCKTIME <  stop_time  (in seconds)
 AND TANK id(tank_index) LEVEL [BELOW|ABOVE] level
 THEN PUMP id(pump_index) STATUS IS status

 In case that start_time > stop_time, then adds one rule for the
 interval [start_time, SECperDAY - 1] and another for [0, stop_time].
 In case that start_time == stop_time, don't do anything.
 **----------------------------------------------------------------
*/
int DLLEXPORT ENaddleveltrigger(int pump_index, int tank_index,
                                int start_time, int stop_time,
                                float level, int status)
{
    int error;

    start_time = start_time % SECperDAY;
    stop_time  = stop_time  % SECperDAY;

    if (start_time > stop_time) {
        if (start_time != SECperDAY - 1) {
            error = addlevelbasedtrigger (pump_index, tank_index,
                                          start_time, SECperDAY - 1,
                                          level, status);
            if (error) return (error);
        }
        if (stop_time != 0)
            return addlevelbasedtrigger (pump_index, tank_index,
                                         0, stop_time,
                                         level, status);
    } else if (start_time < stop_time)
        return addlevelbasedtrigger (pump_index, tank_index,
                                     start_time, stop_time,
                                     level, status);
    return 0;
}
#endif
/*----------------------------------------------------------------
**  Input:   node index of a pump
**  Output:  shortest time interval (seconds) that pump was idle
**----------------------------------------------------------------
*/
int DLLEXPORT ENgetminstoptime(ENsimulation_t sim, int index, int *value)
{
    int s,  pump_idx, schedule_idx, errorcode;
    int *schedule;

    float initial_status;

    int temp_value;

    pump_idx = PUMPINDEX(index);
    schedule_idx = sim->Pump[pump_idx].Schedule_idx;
    schedule = sim->Pump[pump_idx].Schedule;

    errorcode = ENgetlinkvalue (sim, index, EN_INITSTATUS, &initial_status);
    if(errorcode > 0) return (errorcode);

    if (schedule_idx == 0) { // No pump status changes
        *value = (initial_status == 1) ? (0) : (SECperDAY);

    } else if (schedule_idx == 1) {
        *value = (initial_status == 1)
            ? ( SECperDAY - schedule[0] )
            : ( schedule[0] );

    } else {
        // Initial idle interval
        *value = (initial_status == 1)
            ? (schedule[1] - schedule[0])
            : (schedule[0]);

        if (schedule_idx % 2 == initial_status) {
            // Last status change turns off, temp_value is idle time
            temp_value = SECperDAY - schedule[schedule_idx-1];

            if (initial_status == 1 && schedule[0] > 0 ) {
                // No cycle: choose shortest
                if (temp_value < *value) *value = temp_value;
            } else { // Cycle: sum final and initial interval
                *value = *value + temp_value;
            }
        }

        // initial_status == 0 && Pump[pump_idx].Schedule[0] == 0
        if (*value == 0) *value = SECperDAY; // set to maximum

        for ( s = 2 + initial_status; s < schedule_idx; s += 2) {
            temp_value = schedule[s] - schedule[s-1];
            if (temp_value < *value)  *value = temp_value;
        }
    }
    return 0;
}

int DLLEXPORT ENrulesclear(void)
{
        clearrules();
        Nrules = 0;
        MaxRules = 0;
        return(0);
}

int DLLEXPORT ENgetnumwarnings(ENsimulation_t sim)
{
    return sim->Nwarnings;
}
int  DLLEXPORT ENgettotaldemand(ENsimulation_t sim, float *demand)
{
    *demand = sim->TotalSystemDemand * Ucf[DEMAND];
    return 0;
}

int  DLLEXPORT ENgettotalleakage(ENsimulation_t sim, float *leakage)
{
    *leakage = sim->TotalSystemLeakage * Ucf[VOLUME];
    return 0;
}

int  DLLEXPORT ENgettotalinflow(ENsimulation_t sim, float *inflow)
{
    *inflow = sim->TotalSystemInflow * Ucf[FLOW];
    return 0;
}

int DLLEXPORT ENgetnode_xcoord(int index, float *x)
{
    if (!Openflag) return(102);
    if (index <= 0 || index > Nnodes) return(203);

    if(x) *x = Node[index].x;
    return 0;
}

int DLLEXPORT ENgetnode_ycoord(int index, float *y)
{
    if (!Openflag) return(102);
    if (index <= 0 || index > Nnodes) return(203);

    if (y) *y = Node[index].y;
    return 0;
}

/*************************** END OF EPANET.C ***************************/
