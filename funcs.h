/*
**************************************************************************

FUNCS.H -- Function Prototypes for EPANET Program

AUTHOR:     M. Lopez-Ibanez (M. L\'opez-Ib\'a\~nez)
            Napier University, Edinburgh, UK.

$LastChangedBy: manu $ $Revision: 148 $
$LastChangedDate: 2007-11-20 17:29:03 +0100 (Tue, 20 Nov 2007) $

---------------------------------------------------------------------

 Copyright (c) 2005, 2006 Manuel Lopez-Ibanez
 TeX: \copyright 2005, 2006 Manuel L\'opez-Ib\'a\~nez

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
DATE:       5/8/00
            9/25/00
            10/25/00
            12/29/00
            3/1/01
AUTHOR:     L. Rossman
            US EPA - NRMRL

**************************************************************************
*/

/* ------- EPANET.C --------------------*/
/*
**  NOTE: The exportable functions that can be called
**        via the DLL are prototyped in TOOLKIT.H.
*/
void    initpointers(void);               /* Initializes pointers       */
int     allocdata(void);                  /* Allocates memory           */
void    freeTmplist(STmplist *);          /* Frees items in linked list */
void    freeFloatlist(SFloatlist *);      /* Frees list of floats       */
void    freedata(void);                   /* Frees allocated memory     */
int     openfiles(char *,char *,char *);  /* Opens input & report files */
int     strcomp(char *, char *);          /* Compares two strings       */
float   interp(int, float *,              /* Interpolates a data curve  */
               float *, float);
int     findnode(char *);                 /* Finds node's index from ID */
int     findlink(char *);                 /* Finds link's index from ID */
#if 0
void    errmsg(int);                      /* Reports program error      */
#endif
void    writecon(char *);                 /* Writes text to console     */
void    writewin(char *);                 /* Passes text to calling app */

/* ------- INPUT1.C --------------------*/
int     getdata(void);                    /* Gets network data          */
void    setdefaults(void);                /* Sets default values        */
#if 0
void    initreport(void);                 /* Initializes report options */
#endif
void    adjustdata(void);                 /* Adjusts input data         */
int     inittanks(void);                  /* Initializes tank levels    */
void    initunits(void);                  /* Determines reporting units */

/* -------- INPUT2.C -------------------*/
int     netsize(void);                    /* Determines network size    */
int     readdata(void);                   /* Reads in network data      */
int     newline(int, char *);             /* Processes new line of data */
int     addnodeID(int, char *);           /* Adds node ID to data base  */
int     addlinkID(int, char *);           /* Adds link ID to data base  */
int     addpattern(char *);               /* Adds pattern to data base  */
int     addcurve(char *);                 /* Adds curve to data base    */
STmplist *findID(char *, STmplist *);     /* Locates ID on linked list  */
int     unlinked(void);                   /* Checks for unlinked nodes  */
int     getpumpparams(void);              /* Computes pump curve coeffs.*/
int     getpatterns(void);                /* Gets pattern data from list*/
int     getcurves(void);                  /* Gets curve data from list  */
int     findmatch(char *,char *[]);       /* Finds keyword in line      */
int     match(char *, char *);            /* Checks for word match      */
int     gettokens(char *);                /* Tokenizes input line       */
int     getfloat(char *, float *);        /* Converts string to float   */
float   hour(char *, char *);             /* Converts time to hours     */
#if 0
int     setreport(char *);                /* Processes reporting command*/
#endif

/* ---------- INPUT3.C -----------------*/
int     juncdata(void);                   /* Processes junction data    */
int     tankdata(void);                   /* Processes tank data        */
int     pipedata(void);                   /* Processes pipe data        */
int     pumpdata(void);                   /* Processes pump data        */
int     valvedata(void);                  /* Processes valve data       */
int     patterndata(void);                /* Processes pattern data     */
int     curvedata(void);                  /* Processes curve data       */
int     demanddata(void);                 /* Processes demand data      */
int     controldata(void);                /* Processes simple controls  */
int     energydata(void);                 /* Processes energy data      */
int     sourcedata(void);                 /* Processes source data      */
int     emitterdata(void);                /* Processes emitter data     */
#if 0
int     qualdata(void);                   /* Processes quality data     */
int     reactdata(void);                  /* Processes reaction data    */
int     mixingdata(void);                 /* Processes tank mixing data */
#endif
int     statusdata(void);                 /* Processes link status data */
#if 0
int     reportdata(void);                 /* Processes report options   */
#endif
int     timedata(void);                   /* Processes time options     */
int     optiondata(void);                 /* Processes analysis options */
int     optionchoice(int);                /* Processes option choices   */
int     optionvalue(int);                 /* Processes option values    */
int     powercurve(float, float, float,   /* Coeffs. of power pump curve*/
                   float, float, float *,
                   float *, float *);
int     valvecheck(int, int, int);        /* Checks valve placement     */
void    changestatus(int, char, float);   /* Changes status of a link   */
int     coordinates(void);                // Process coordinates of nodes

/* -------------- RULES.C --------------*/
void    initrules(void);                  /* Initializes rule base      */
void    addrule(char *);                  /* Adds rule to rule base     */
int     allocrules(void);                 /* Allocates memory for rule  */
int     ruledata(void);                   /* Processes rule input data  */
int     checkrules(ENsimulation_t, long); /* Checks all rules           */
void    freerules(void);                  /* Frees rule base memory     */
int     addlevelbasedtrigger(int pump_index, int tank_index,
                             int start_time, int stop_time,
                             float level, int status);
void    clearrules(void);


/* ------------- REPORT.C --------------*/
#if 0
int     writereport(void);                /* Writes formatted report    */
void    writelogo(void);                  /* Writes program logo        */
void    writesummary(void);               /* Writes network summary     */
#endif
void    writehydstat(ENsimulation_t,      /* Writes hydraulic status    */
                     int,float);    
#if 0
void    writeenergy(void);                /* Writes energy usage        */
int     writeresults(void);               /* Writes node/link results   */
void    writeheader(int,int);             /* Writes heading on report   */
void    writeline(char *);                /* Writes line to report file */
void    writerelerr(int, float);          /* Writes convergence error   */
void    writestatchange(int,char,char);   /* Writes link status change  */
void    writecontrolaction(int, int);     /* Writes control action taken*/
void    writeruleaction(int, char *);     /* Writes rule action taken   */
#endif
int     writehydwarn(ENsimulation_t,      /* Writes hydraulic warnings  */
                     int,float);    
void    writehyderr(ENsimulation_t, int); /* Writes hydraulic error msg.*/
#if 0
void    writelimits(int,int);             /* Writes reporting limits    */
#endif
int     checklimits(float *,int,int);     /* Checks variable limits     */
#if 0
void    writetime(char *);                /* Writes current clock time  */
#endif
char    *clocktime(ENclockstring_t *, long);/* Converts time to hrs:min   */
char    *fillstr(char *, char, int);      /* Fills string with character*/
int     getnodetype(int);                 /* Determines node type       */

/* --------- HYDRAUL.C -----------------*/
ENsimulation_t openhyd(void);             /* Opens hydraulics solver    */

/*** Updated 3/1/01 ***/
void    inithyd(ENsimulation_t, int);     /* Re-sets initial conditions */

int     runhyd(ENsimulation_t, long *);   /* Solves 1-period hydraulics */
int     nexthyd(ENsimulation_t, long *);  /* Moves to next time period  */
void    closehyd(ENsimulation_t);         /* Closes hydraulics solver   */
void    getenergy(ENsimulation_t, int,    /* Computes link energy use   */
                  float *, float *); 
void    setlinkstatus(ENsimulation_t, int,/* Sets link status           */ 
                      char, char *, float *);
void    setlinksetting(ENsimulation_t, int,/* Sets pump/valve setting    */
                       float, char *, float *);
void    resistance(int);                  /* Computes resistance coeff. */
float   tankvolume(int,float);            /* Finds tank vol. from grade */


#if 0
/* ----------- QUALITY.C ---------------*/
int     openqual(void);                   /* Opens WQ solver system     */
void    initqual(void);                   /* Initializes WQ solver      */
int     runqual(long *);                  /* Gets current WQ results    */
int     nextqual(long *);                 /* Updates WQ by hyd.timestep */
int     stepqual(long *);                 /* Updates WQ by WQ time step */
int     closequal(void);                  /* Closes WQ solver system    */
int     gethyd(long *, long *);           /* Gets next hyd. results     */
char    setReactflag(void);               /* Checks for reactive chem.  */
void    transport(long);                  /* Transports mass in network */
void    initsegs(void);                   /* Initializes WQ segments    */
void    reorientsegs(void);               /* Re-orients WQ segments     */
void    updatesegs(long);                 /* Updates quality in segments*/
void    removesegs(int);                  /* Removes a WQ segment       */
void    addseg(int,float,float);          /* Adds a WQ segment to pipe  */
void    accumulate(long);                 /* Sums mass flow into node   */
void    updatenodes(long);                /* Updates WQ at nodes        */
void    sourceinput(long);                /* Computes source inputs     */
void    release(long);                    /* Releases mass from nodes   */
void    updatetanks(long);                /* Updates WQ in tanks        */
void    updatesourcenodes(long);          /* Updates WQ at source nodes */
void    tankmix1(int, long);              /* Complete mix tank model    */
void    tankmix2(int, long);              /* 2-compartment tank model   */
void    tankmix3(int, long);              /* FIFO tank model            */
void    tankmix4(int, long);              /* LIFO tank model            */
float   sourcequal(Psource);              /* Finds WQ input from source */
float   avgqual(int);                     /* Finds avg. quality in pipe */
void    ratecoeffs(void);                 /* Finds wall react. coeffs.  */
float   piperate(int);                    /* Finds wall react. coeff.   */
float   pipereact(int,float,float,long);  /* Reacts water in a pipe     */
float   tankreact(float,float,float,
                  long);                  /* Reacts water in a tank     */
float   bulkrate(float,float,float);      /* Finds bulk reaction rate   */
float   wallrate(float,float,float,float);/* Finds wall reaction rate   */
#endif

/* ------------ OUTPUT.C ---------------*/
#if 0
int     savenetdata(void);                /* Saves basic data to file   */
int     savehyd(long *);                  /* Saves hydraulic solution   */
int     savehydstep(long *);              /* Saves hydraulic timestep   */
int     saveenergy(void);                 /* Saves energy usage         */
int     readhyd(long *);                  /* Reads hydraulics from file */
int     readhydstep(long *);              /* Reads time step from file  */
int     saveoutput(void);                 /* Saves results to file      */
int     nodeoutput(int,float *,float);    /* Saves node results to file */
int     linkoutput(int,float *,float);    /* Saves link results to file */
int     savefinaloutput(void);            /* Finishes saving output     */
#endif
int     savetimestat(float *, char);      /* Saves time stats to file   */
int     savenetreacts(REAL, REAL, REAL,
                      REAL);              /* Saves react. rates to file */
#if 0
int     saveepilog(void);                 /* Saves output file epilog   */
#endif

/* ------------ INPFILE.C ---------------*/
#if 0
int     saveinpfile(char *);             /* Saves network to text file  */
#endif
