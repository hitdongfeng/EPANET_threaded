/*
*******************************************************************************

QUALITY.C -- Water Quality Simulator for EPANET Program

AUTHOR:     M. Lopez-Ibanez (M. L\'opez-Ib\'a\~nez)
            Napier University, Edinburgh, UK.

$LastChangedBy: manu $ $Revision: 136 $
$LastChangedDate: 2007-11-16 14:20:36 +0100 (Fri, 16 Nov 2007) $

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
DATE:       5/29/00
            9/7/00
            10/25/00
AUTHOR:     L. Rossman
            US EPA - NRMRL

  This module contains the network water quality simulator.

  For each time period, hydraulic results are read in from the
  binary file HydFile, hydraulic and water quality results are
  written to the binary output file OutFile (if the current period
  is a reporting period), and the water quality is transported
  and reacted over the duration of the time period.

  The entry points for this module are:
    openqual()   -- called from ENopenQ() in EPANET.C
    initqual()   -- called from ENinitQ() in EPANET.C
    runqual()    -- called from ENrunQ() in EPANET.C
    nextqual()   -- called from ENnextQ() in EPANET.C
    stepqual()   -- called from ENstepQ() in EPANET.C
    closequal()  -- called from ENcloseQ() in EPANET.C

  Calls are made to:
    AllocInit()
    Alloc()
    AllocFree()
  in MEMPOOL.C to utilize a memory pool to prevent excessive malloc'ing
  when constantly creating and destroying pipe sub-segments during
  the water quality transport calculations.

  Calls are also made to:
    readhyd()
    readhydstep()
    savenetdata()
    saveoutput()
    savefinaloutput()
  in OUTPUT.C to retrieve hydraulic results and save all results.

*******************************************************************************
*/
#if 0
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include "toolkit.h"
#include "hash.h"
#include "text.h"
#include "types.h"
#include "funcs.h"
#define  EXTERN  extern
#include "vars.h"
#include "mempool.h"

/*
** Macros to identify upstream & downstream nodes of a link
** under the current flow and to compute link volume
*/
#define   UP_NODE(x)   ( (FlowDir[(x)]=='+') ? Link[(x)].N1 : Link[(x)].N2 )
#define   DOWN_NODE(x) ( (FlowDir[(x)]=='+') ? Link[(x)].N2 : Link[(x)].N1 )
#define   LINKVOL(k)   ( 0.785398*Link[(k)].Len*SQR(Link[(k)].Diam) )

Pseg      FreeSeg;              /* Pointer to unused segment               */
Pseg      *FirstSeg,            /* First (downstream) segment in each pipe */
          *LastSeg;             /* Last (upstream) segment in each pipe    */
char      *FlowDir;             /* Flow direction for each pipe            */
float     *VolIn;               /* Total volume inflow to node             */
float     *MassIn;              /* Total mass inflow to node               */
float     Sc;                   /* Schmidt Number                          */
float     Bucf;                 /* Bulk reaction units conversion factor   */
float     Tucf;                 /* Tank reaction units conversion factor   */
char      Reactflag;            /* Reaction indicator                      */
char      OutOfMemory;          /* Out of memory indicator                 */


int  openqual()
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  returns error code
**   Purpose: opens WQ solver system
**--------------------------------------------------------------
*/
{
    int errcode = 0;
    int n;

    /* Allocate memory pool for WQ segments */
    OutOfMemory = FALSE;
    if (AllocInit() == NULL) errcode = 101;

    simulation->C    = calloc(MaxNodes + 1, sizeof(float));
    ERRCODE(MEMCHECK(simulation->C));

    /* Allocate reaction rate array */
    simulation->R  = calloc ((Nlinks+1), sizeof(float));
    ERRCODE(MEMCHECK(simulation->R));

    /* Allocate memory for WQ solver */
    n        = Nlinks+Ntanks+1;
    FirstSeg = (Pseg *) calloc(n, sizeof(Pseg));
    LastSeg  = (Pseg *) calloc(n, sizeof(Pseg));
    FlowDir  = (char *) calloc(n, sizeof(char));
    n        = Nnodes+1;
    VolIn    = (float *) calloc(n, sizeof(float));
    MassIn   = (float *) calloc(n, sizeof(float));
    ERRCODE(MEMCHECK(FirstSeg));
    ERRCODE(MEMCHECK(LastSeg));
    ERRCODE(MEMCHECK(FlowDir));
    ERRCODE(MEMCHECK(VolIn));
    ERRCODE(MEMCHECK(MassIn));
    return(errcode);
}


/* Local function to compute unit conversion factor for bulk reaction rates */
float getucf(float order)
{
    if (order < 0.0) order = 0.0;
    if (order == 1.0) return(1.0);
    else return(1./pow(LperFT3,(order-1.0)));
}


void  initqual()
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  none
**   Purpose: re-initializes WQ solver system
**--------------------------------------------------------------
*/
{
    int i;

    /* Initialize quality, tank volumes, & source mass flows */
    for (i=1; i<=Nnodes; i++) simulation->C[i] = Node[i].C0;
    for (i=1; i<=Ntanks; i++) Tank[i].C = Node[Tank[i].Node].C0;
    for (i=1; i<=Ntanks; i++) simulation->V[i] = Tank[i].V0;
    for (i=1; i<=Nnodes; i++)
        if (Node[i].S != NULL) Node[i].S->Smass = 0.0;

    /* Set WQ parameters */
    Bucf = 1.0;
    Tucf = 1.0;
    Reactflag = 0;
    if (Qualflag != NONE)
    {
        /* Initialize WQ at trace node (if applicable) */
        if (Qualflag == TRACE) simulation->C[TraceNode] = 100.0;

        /* Compute Schmidt number */
        if (Diffus > 0.0)
            Sc = Viscos/Diffus;
        else
            Sc = 0.0;

        /* Compute unit conversion factor for bulk react. coeff. */
        Bucf = getucf(BulkOrder);
        Tucf = getucf(TankOrder);

        /* Check if modeling a reactive substance */
        Reactflag = setReactflag();

        /* Reset memory pool */
        FreeSeg = NULL;
        AllocReset();
    }

    /* Initialize avg. reaction rates */
    Wbulk = 0.0;
    Wwall = 0.0;
    Wtank = 0.0;
    Wsource = 0.0;
#if 0
    /* Re-position hydraulics file */
    fseek(HydFile,HydOffset,SEEK_SET);
#endif
    /* Set elapsed times to zero */
    simulation->Htime = 0;
    Qtime = 0;
    simulation->Rtime = Rstart;
    Nperiods = 0;
}


int runqual(long *t)
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  t = pointer to current simulation time (sec)
**   Returns: error code
**   Purpose: retrieves hydraulics for next hydraulic time step
**            (at time *t) and saves current results to file
**--------------------------------------------------------------
*/
{
    long    hydtime;       /* Hydraulic solution time */
    long    hydstep;       /* Hydraulic time step     */
    int     errcode = 0;

    /* Update reported simulation time */
    *t = Qtime;

    /* Read hydraulic solution from hydraulics file */
    if (Qtime == simulation->Htime)
    {
        errcode = gethyd(&hydtime, &hydstep);
        simulation->Htime = hydtime + hydstep;
    }
    return(errcode);
}


int nextqual(long *tstep)
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  tstep = pointer to time step (sec)
**   Returns: error code
**   Purpose: updates WQ conditions until next hydraulic
**            solution occurs (after *tstep secs.)
**--------------------------------------------------------------
*/
{
    long    hydstep;       /* Hydraulic solution time step */
    int     errcode = 0;

    /* Determine time step */
    *tstep = 0;
    hydstep = simulation->Htime - Qtime;

    /* Perform water quality routing over this time step */
    if (Qualflag != NONE && hydstep > 0) transport(hydstep);

    /* Update current time */
    if (OutOfMemory) errcode = 101;
    if (!errcode) *tstep = hydstep;
    Qtime += hydstep;
#if 0
    /* Save final output if no more time steps */
    if (!errcode && simulation->Saveflag && *tstep == 0) 
        errcode = savefinaloutput();
#endif
    return(errcode);
}


int stepqual(long *tleft)
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  tleft = pointer to time left in simulation
**   Returns: error code
**   Purpose: updates WQ conditions over a single WQ time step
**--------------------------------------------------------------
*/
{  long dt, hstep, t, tstep;
int  errcode = 0;
tstep = Qstep;
do
{
    dt = tstep;
    hstep = simulation->Htime - Qtime;
    if (hstep < dt)
    {
        dt = hstep;
        if (Qualflag != NONE) transport(dt);
        Qtime += dt;
        errcode = runqual(&t);
        Qtime = t;
    }
    else
    {
        if (Qualflag != NONE) transport(dt);
        Qtime += dt;
    }
    tstep -= dt;
    if (OutOfMemory) errcode = 101;
}  while (!errcode && tstep > 0);
*tleft = Dur - Qtime;
#if 0
if (!errcode && simulation->Saveflag && *tleft == 0) 
    errcode = savefinaloutput();
#endif
return(errcode);
}


int closequal()
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  returns error code
**   Purpose: closes WQ solver system
**--------------------------------------------------------------
*/
{
    int errcode = 0;
    AllocFreePool();   /* Free memory pool */
    free (simulation->C);
    free(FirstSeg);
    free(LastSeg);
    free(FlowDir);
    free(VolIn);
    free(MassIn);
    free (simulation->R);
    return(errcode);
}


int  gethyd(long *hydtime, long *hydstep)
/*
**-----------------------------------------------------------
**   Input:   none
**   Output:  hydtime = pointer to hydraulic solution time
**            hydstep = pointer to hydraulic time step
**   Returns: error code
**   Purpose: retrieves hydraulic solution and hydraulic
**            time step for next hydraulic event
**
**   NOTE: when this function is called, WQ results have
**         already been updated to the point in time when
**         the next hydraulic event occurs.
**-----------------------------------------------------------
*/
{
    int errcode = 0;
#if 0
    /* Read hydraulic results from file */
    if (!readhyd(hydtime)) return(307);
    if (!readhydstep(hydstep)) return(307);
#endif
    simulation->Htime = *hydtime;

    /* Save current results to output file */
    if (simulation->Htime >= simulation->Rtime)
    {
        if (simulation->Saveflag)
        {
            errcode = saveoutput();
            Nperiods++;
        }
        simulation->Rtime += Rstep;
    }

    /* If simulating WQ: */
    if (Qualflag != NONE && Qtime < Dur)
    {

        /* Compute reaction rate coeffs. */
        if (Reactflag && Qualflag != AGE) ratecoeffs();

        /* Initialize pipe segments (at time 0) or  */
        /* else re-orient segments if flow reverses.*/
        if (Qtime == 0) initsegs();
        else            reorientsegs();
    }
    return(errcode);
}


char  setReactflag()
/*
**-----------------------------------------------------------
**   Input:   none
**   Output:  returns 1 for reactive WQ constituent, 0 otherwise
**   Purpose: checks if reactive chemical being simulated
**-----------------------------------------------------------
*/
{
    int  i;
    if      (Qualflag == TRACE) return(0);
    else if (Qualflag == AGE)   return(1);
    else
    {
        for (i=1; i<=Nlinks; i++)
        {
            if (Link[i].Type <= PIPE)
            {
                if (Link[i].Kb != 0.0 || Link[i].Kw != 0.0) return(1);
            }
        }
        for (i=1; i<=Ntanks; i++)
            if (Tank[i].Kb != 0.0) return(1);
    }
    return(0);
}


void  transport(long tstep)
/*
**--------------------------------------------------------------
**   Input:   tstep = length of current time step
**   Output:  none
**   Purpose: transports constituent mass through pipe network
**            under a period of constant hydraulic conditions.
**--------------------------------------------------------------
*/
{
    long   qtime, dt;

    /* Repeat until elapsed time equals hydraulic time step */
    qtime = 0;
    while (!OutOfMemory && qtime < tstep)
    {                                  /* Qstep is quality time step */
        dt = MIN(Qstep,tstep-qtime);    /* Current time step */
        qtime += dt;                    /* Update elapsed time */
        if (Reactflag) updatesegs(dt);  /* Update quality in inner link segs */
        accumulate(dt);                 /* Accumulate flow at nodes */
        updatenodes(dt);                /* Update nodal quality */
        sourceinput(dt);                /* Compute inputs from sources */
        release(dt);                    /* Release new nodal flows */
    }
    updatesourcenodes(tstep);          /* Update quality at source nodes */
}


void  initsegs()
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  none
**   Purpose: initializes water quality segments
**--------------------------------------------------------------
*/
{
    int     j,k;
    float   c,v;

    /* Examine each link */
    for (k=1; k<=Nlinks; k++)
    {

        /* Establish flow direction */
        FlowDir[k] = '+';
        if (simulation->Q[k] < 0.) FlowDir[k] = '-';

        /* Set segs to zero */
        LastSeg[k] = NULL;
        FirstSeg[k] = NULL;

        /* Find quality of downstream node */
        j = DOWN_NODE(k);
        if (j <= Njuncs) c = simulation->C[j];
        else             c = Tank[j-Njuncs].C;

        /* Fill link with single segment with this quality */
        addseg(k,LINKVOL(k),c);
    }

    /* Initialize segments in tanks that use them */
    for (j=1; j<=Ntanks; j++)
    {

        /* Skip reservoirs & complete mix tanks */
        if (Tank[j].A == 0.0
            ||  Tank[j].MixModel == MIX1) continue;

        /* Tank segment pointers are stored after those for links */
        k = Nlinks + j;
        c = Tank[j].C;
        LastSeg[k] = NULL;
        FirstSeg[k] = NULL;

        /* Add 2 segments for 2-compartment model */
        if (Tank[j].MixModel == MIX2)
        {
            v = MAX (0, simulation->V[j] - Tank[j].V1max);
            addseg(k,v,c);
            v = simulation->V[j] - v;
            addseg(k,v,c);
        }

        /* Add one segment for FIFO & LIFO models */
        else
        {
            v = simulation->V[j];
            addseg(k,v,c);
        }
    }
}


void  reorientsegs()
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  none
**   Purpose: re-orients segments (if flow reverses)
**--------------------------------------------------------------
*/
{
    Pseg   seg, nseg, pseg;
    int    k;
    char   newdir;

    /* Examine each link */
    for (k=1; k<=Nlinks; k++)
    {

        /* Find new flow direction */
        newdir = '+';
        if (simulation->Q[k] == 0.0)     newdir = FlowDir[k];
        else if (simulation->Q[k] < 0.0) newdir = '-';

        /* If direction changes, then reverse order of segments */
        /* (first to last) and save new direction */
        if (newdir != FlowDir[k])
        {
            seg = FirstSeg[k];
            FirstSeg[k] = LastSeg[k];
            LastSeg[k] = seg;
            pseg = NULL;
            while (seg != NULL)
            {
                nseg = seg->prev;
                seg->prev = pseg;
                pseg = seg;
                seg = nseg;
            }
            FlowDir[k] = newdir;
        }
    }
}


void  updatesegs(long dt)
/*
**-------------------------------------------------------------
**   Input:   t = time from last WQ segment update
**   Output:  none
**   Purpose: reacts material in pipe segments up to time t
**-------------------------------------------------------------
*/
{
    int    k;
    Pseg   seg;
    float  cseg, rsum, vsum;

    /* Examine each link in network */
    for (k=1; k<=Nlinks; k++)
    {

        /* Skip zero-length links (pumps & valves) */
        rsum = 0.0;
        vsum = 0.0;
        if (Link[k].Len == 0.0) continue;

        /* Examine each segment of the link */
        seg = FirstSeg[k];
        while (seg != NULL)
        {

            /* React segment over time dt */
            cseg = seg->c;
            seg->c = pipereact(k,seg->c,seg->v,dt);

            /* Accumulate volume-weighted reaction rate */
            if (Qualflag == CHEM)
            {
                rsum += ABS((seg->c - cseg))*seg->v;
                vsum += seg->v;
            }
            seg = seg->prev;
        }

        /* Normalize volume-weighted reaction rate */
        if (vsum > 0.0) simulation->R[k] = rsum / vsum / dt * SECperDAY;
        else simulation->R[k] = 0.0;
    }
}


void  removesegs(int k)
/*
**-------------------------------------------------------------
**   Input:   k = link index
**   Output:  none
**   Purpose: removes all segments in link k
**-------------------------------------------------------------
*/
{
    Pseg seg;
    seg = FirstSeg[k];
    while (seg != NULL)
    {
        FirstSeg[k] = seg->prev;
        seg->prev = FreeSeg;
        FreeSeg = seg;
        seg = FirstSeg[k];
    }
    LastSeg[k] = NULL;
}


void  addseg(int k, float v, float c)
/*
**-------------------------------------------------------------
**   Input:   k = link segment
**            v = segment volume
**            c = segment quality
**   Output:  none
**   Purpose: adds a segment to start of link k (i.e., upstream
**            of current last segment).
**-------------------------------------------------------------
*/
{
    Pseg seg;

    if (FreeSeg != NULL)
    {
        seg = FreeSeg;
        FreeSeg = seg->prev;
    }
    else
    {
        seg = (struct Sseg *) Alloc(sizeof(struct Sseg));
        if (seg == NULL)
        {
            OutOfMemory = TRUE;
            return;
        }
    }
    seg->v = v;
    seg->c = c;
    seg->prev = NULL;
    if (FirstSeg[k] == NULL) FirstSeg[k] = seg;
    if (LastSeg[k] != NULL) LastSeg[k]->prev = seg;
    LastSeg[k] = seg;
}


void accumulate(long dt)
/*
**-------------------------------------------------------------
**   Input:   dt = current WQ time step
**   Output:  none
**   Purpose: accumulates mass flow at nodes and updates nodal
**            quality
**-------------------------------------------------------------
*/
{
    int    i,j,k;
    float  cseg,v,vseg;
    Pseg   seg;

    /* Re-set memory used to accumulate mass & volume */
    memset(VolIn,0,(Nnodes+1)*sizeof(float));
    memset(MassIn,0,(Nnodes+1)*sizeof(float));
    memset(simulation->X, 0, (Nnodes+1) * sizeof(float));

    /* Compute average conc. of segments adjacent to each node */
    /* (For use if there is no transport through the node) */
    for (k=1; k<=Nlinks; k++)
    {
        j = DOWN_NODE(k);             /* Downstream node */
        if (FirstSeg[k] != NULL)      /* Accumulate concentrations */
        {
            MassIn[j] += FirstSeg[k]->c;
            VolIn[j]++;
        }
        j = UP_NODE(k);              /* Upstream node */
        if (LastSeg[k] != NULL)      /* Accumulate concentrations */
        {
            MassIn[j] += LastSeg[k]->c;
            VolIn[j]++;
        }
    }
    for (k=1; k<=Nnodes; k++)
        if (VolIn[k] > 0.0) simulation->X[k] = MassIn[k]/VolIn[k];

    /* Move mass from first segment of each pipe into downstream node */
    memset(VolIn,0,(Nnodes+1)*sizeof(float));
    memset(MassIn,0,(Nnodes+1)*sizeof(float));
    for (k=1; k<=Nlinks; k++)
    {
        i = UP_NODE(k);               /* Upstream node */
        j = DOWN_NODE(k);             /* Downstream node */
        v = ABS (simulation->Q[k])*dt;/* Flow volume */

        /* If link volume < flow volume, then transport upstream    */
        /* quality to downstream node and remove all link segments. */
        if (LINKVOL(k) < v)
        {
            VolIn[j] += v;
            seg = FirstSeg[k];
            cseg = simulation->C[i];
            if (seg != NULL) cseg = seg->c;
            MassIn[j] += v*cseg;
            removesegs(k);
        }

        /* Otherwise remove flow volume from leading segments */
        /* and accumulate flow mass at downstream node        */
        else while (v > 0.0)
        {
            /* Identify leading segment in pipe */
            seg = FirstSeg[k];
            if (seg == NULL) break;

            /* Volume transported from this segment is */
            /* minimum of flow volume & segment volume */
            /* (unless leading segment is also last segment) */
            vseg = seg->v;
            vseg = MIN(vseg,v);
            if (seg == LastSeg[k]) vseg = v;

            /* Update volume & mass entering downstream node  */
            cseg = seg->c;
            VolIn[j] += vseg;
            MassIn[j] += vseg*cseg;

            /* Reduce flow volume by amount transported */
            v -= vseg;

            /* If all of segment's volume was transferred, then */
            /* replace leading segment with the one behind it   */
            /* (Note that the current seg is recycled for later use.) */
            if (v >= 0.0 && vseg >= seg->v)
            {
                FirstSeg[k] = seg->prev;
                if (FirstSeg[k] == NULL) LastSeg[k] = NULL;
                seg->prev = FreeSeg;
                FreeSeg = seg;
            }

            /* Otherwise reduce segment's volume */
            else
            {
                seg->v -= vseg;
            }
        }     /* End while */
    }        /* Next link */
}


void updatenodes(long dt)
/*
**---------------------------------------------------------------------------
**   Input:   dt = current WQ time step
**   Output:  none
**   Purpose: updates concentration at all nodes to mixture of accumulated
**            inflow from connecting pipes.
**
**  Note:     Does not account for source flow effects. X[i] contains
**            average concen. of segments adjacent to node i, used in case
**            there was no inflow into i.
**---------------------------------------------------------------------------
*/
{
    int i;

    /* Update junction quality */
    for (i=1; i<=Njuncs; i++)
    {
        if (simulation->D[i] < 0.0) VolIn[i] -= simulation->D[i] * dt;
        if (VolIn[i] > 0.0) simulation->C[i] = MassIn[i] / VolIn[i];
        else                simulation->C[i] = simulation->X[i];
    }

    /* Update tank quality */
    updatetanks(dt);

    /* For flow tracing, set source node concen. to 100. */
    if (Qualflag == TRACE) simulation->C[TraceNode] = 100.0;
}


void sourceinput(long dt)
/*
**---------------------------------------------------------------------
**   Input:   dt = current WQ time step
**   Output:  none
**   Purpose: computes contribution (if any) of mass additions from WQ
**            sources at each node.
**---------------------------------------------------------------------
*/
{
    int   j,n;
    float massadded, s, volout;
    float qout, qcutoff;
    Psource source;

    /* Establish a flow cutoff which indicates no outflow from a node */
    qcutoff = 10.0*TINY;

    /* Zero-out the work array X */
    memset (simulation->X, 0, (Nnodes + 1) * sizeof(float));
    if (Qualflag != CHEM) return;

    /* Consider each node */
    for (n=1; n<=Nnodes; n++)
    {

        /* Skip node if no WQ source */
        source = Node[n].S;
        if (source == NULL) continue;
        if (source->C0 == 0.0) continue;

        /* Find total flow volume leaving node */
        if (n <= Njuncs) volout = VolIn[n];  /* Junctions */
        else volout = VolIn[n] - simulation->D[n] * dt;    /* Tanks */
        qout = volout / (float) dt;

        /* Evaluate source input only if node outflow > cutoff flow */
        if (qout > qcutoff)
        {

            /* Mass added depends on type of source */
            s = sourcequal(source);
            switch(source->Type)
            {
                /* Concen. Source: */
                /* Mass added = source concen. * -(demand) */
            case CONCEN:

                /* Only add source mass if demand is negative */
                if (simulation->D[n] < 0.0)
                {
                    massadded = -s * simulation->D[n] * dt;

                    /* If node is a tank then set concen. to 0. */
                    /* (It will be re-set to true value in updatesourcenodes()) */
                    if (n > Njuncs) simulation->C[n] = 0.0;
                }
                else massadded = 0.0;
                break;

                /* Mass Inflow Booster Source: */
            case MASS:
                massadded = s*dt;
                break;

                /* Setpoint Booster Source: */
                /* Mass added is difference between source */
                /* & node concen. times outflow volume  */
            case SETPOINT:
                if (s > simulation->C[n]) 
                    massadded = (s - simulation->C[n]) * volout;
                else massadded = 0.0;
                break;

                /* Flow-Paced Booster Source: */
                /* Mass added = source concen. times outflow volume */
            case FLOWPACED:
                massadded = s*volout;
                break;

            default:
                fprintf(stderr,
                        "quality.c: sourceinput(): error in source->Type");
                exit(1);
            }

            /* Source concen. contribution = (mass added / outflow volume) */
            simulation->X[n] = massadded / volout;

            /* Update total mass added for time period & simulation */
            source->Smass += massadded;
            if (simulation->Htime >= Rstart) Wsource += massadded;
        }
    }

    /* Add mass inflows from reservoirs to Wsource*/
    if (simulation->Htime >= Rstart)
    {
        for (j=1; j<=Ntanks; j++)
        {
            if (Tank[j].A == 0.0)
            {
                n = Njuncs + j;
                volout = VolIn[n] - simulation->D[n] * dt;
                if (volout > 0.0) Wsource += volout * simulation->C[n];
            }
        }
    }
}


void release(long dt)
/*
**---------------------------------------------------------
**   Input:   dt = current WQ time step
**   Output:  none
**   Purpose: creates new segments in outflow links from nodes.
**---------------------------------------------------------
*/
{
    int    k,n;
    float  c,q,v;
    Pseg   seg;

    /* Examine each link */
    for (k=1; k<=Nlinks; k++)
    {

        /* Ignore links with no flow */
        if (simulation->Q[k] == 0.0) continue;

        /* Find flow volume released to link from upstream node */
        /* (NOTE: Flow volume is allowed to be > link volume.) */
        n = UP_NODE(k);
        q = ABS (simulation->Q[k]);
        v = q*dt;

        /* Include source contribution in quality released from node. */
        c = simulation->C[n] + simulation->X[n];

        /* If link has a last seg, check if its quality     */
        /* differs from that of the flow released from node.*/
        if ( (seg = LastSeg[k]) != NULL)
        {

            /* Quality is the same, so just add flow volume to last seg */
            if (ABS(seg->c - c) < Ctol) seg->v += v;

            /* Otherwise add a new seg to end of link */
            else addseg(k,v,c);
        }

        /* If link has no segs then add a new one. */
        else addseg(k,LINKVOL(k),c);
    }
}


void  updatesourcenodes(long dt)
/*
**---------------------------------------------------
**   Input:   dt = current WQ time step
**   Output:  none
**   Purpose: updates quality at source nodes.
**            (X[n] = concen. added by source at node n)
**---------------------------------------------------
*/
{
    int i,n;
    Psource source;

    if (Qualflag != CHEM) return;

    /* Examine each WQ source node */
    for (n=1; n<=Nnodes; n++)
    {
        source = Node[n].S;
        if (source == NULL) continue;

        /* Add source to current node concen. */
        simulation->C[n] += simulation->X[n];

        /* For tanks, node concen. = internal concen. */
        if (n > Njuncs)
        {
            i = n - Njuncs;
            if (Tank[i].A > 0.0) simulation->C[n] = Tank[i].C;
        }

        /* Normalize mass added at source to time step */
        source->Smass /= (float)dt;
    }
}


void  updatetanks(long dt)
/*
**---------------------------------------------------
**   Input:   dt = current WQ time step
**   Output:  none
**   Purpose: updates tank volumes & concentrations
**---------------------------------------------------
*/
{
    int   i,n;

    /* Examine each reservoir & tank */
    for (i=1; i<=Ntanks; i++)
    {

        /* Use initial quality for reservoirs */
        if (Tank[i].A == 0.0)
        {
            n = Tank[i].Node;
            simulation->C[n] = Node[n].C0;
        }

        /* Update tank WQ based on mixing model */
        else switch(Tank[i].MixModel)
        {
        case MIX2: tankmix2(i,dt); break;
        case FIFO: tankmix3(i,dt); break;
        case LIFO: tankmix4(i,dt); break;
        default:   tankmix1(i,dt); break;
        }
    }
}


void  tankmix1(int i, long dt)
/*
**---------------------------------------------
**   Input:   i = tank index
**            dt = current WQ time step
**   Output:  none
**   Purpose: complete mix tank model
**---------------------------------------------
*/
{
    int   n;
    float cin;

    /* Blend inflow with contents */
    n = Tank[i].Node;
    if (VolIn[n] > 0.0) cin = MassIn[n]/VolIn[n];
    else                 cin = 0.0;
    if (simulation->V[i] > 0.0)
        Tank[i].C = tankreact (Tank[i].C, simulation->V[i], Tank[i].Kb,dt) 
            + (cin - Tank[i].C) * VolIn[n] / simulation->V[i];
    else Tank[i].C = cin;
    Tank[i].C = MAX(0.0, Tank[i].C);

    /* Update tank volume & nodal quality */
    simulation->V[i] += simulation->D[n] * dt;
    simulation->C[n] = Tank[i].C;
}

/*** Updated 10/25/00 ***/
void  tankmix2(int i, long dt)
/*
**------------------------------------------------
**   Input:   i = tank index
**            dt = current WQ time step
**   Output:  none
**   Purpose: 2-compartment tank model
**            (seg1 = mixing zone,
**             seg2 = ambient zone)
**------------------------------------------------
*/
{
    int     k,n;
    long    tstep, tstar;
    float   cin,        /* Inflow quality */
        qin,        /* Inflow rate */
        qout,       /* Outflow rate */
        qnet;       /* Net flow rate */
    Pseg     seg1,seg2;  /* Compartment segments */

    /* Identify segments for each compartment */
    k = Nlinks + i;
    seg1 = LastSeg[k];
    seg2 = FirstSeg[k];
    if (seg1 == NULL || seg2 == NULL) return;

    /* Find inflows & outflows */
    n = Tank[i].Node;
    qnet = simulation->D[n];
    qin = VolIn[n]/(float)dt;
    qout = qin - qnet;
    if (qin > 0.0) cin = MassIn[n]/VolIn[n];
    else           cin = 0.0;
    simulation->V[i] += qnet*dt;

    /* Case of no net volume change */
    if (ABS(qnet) < TINY)
    {
        seg1->c = tankreact(seg1->c,seg1->v,Tank[i].Kb,dt);
        seg2->c = tankreact(seg2->c,seg2->v,Tank[i].Kb,dt);
    }

    /* Case of net filling (qnet > 0) */
    else if (qnet > 0.0)
    {

        /* Case where ambient zone empty & mixing zone filling */
        if (seg2->v <= 0.0)
        {
            tstar = (long) ((Tank[i].V1max-(seg1->v))/qnet);
            tstep = MIN(dt, tstar);
            if (seg1->v > 0.0)
                seg1->c = tankreact(seg1->c,seg1->v,Tank[i].Kb,tstep) +
                    qin*tstep*(cin-(seg1->c))/(seg1->v);
            else seg1->c = cin;
            seg1->c = MAX(0.0, seg1->c);
            seg1->v += qnet*tstep;
            seg2->c = 0.0;
            dt -= tstep;
        }

        /* Case where mixing zone full & ambient zone filling */
        if (dt > 1)
        {
            seg1->c = tankreact(seg1->c,seg1->v,Tank[i].Kb,dt) +
                qin*dt*(cin-(seg1->c))/(seg1->v);
            seg1->c = MAX(0.0, seg1->c);
            if (seg2->v <= 0.0) seg2->c = seg1->c;
            else
                seg2->c = tankreact(seg2->c,seg2->v,Tank[i].Kb,dt) +
                    qnet*dt*((seg1->c)-(seg2->c))/(seg2->v);
            seg2->c = MAX(0.0, seg2->c);
            seg2->v += qnet*dt;
        }
    }

    /* Case of net emptying (qnet < 0) */
    else if (qnet < 0.0 && seg1->v > 0.0)
    {

        /* Case where mixing zone full & ambient zone draining */
        if ((seg2->v) > 0.0)
        {
            tstar = (long)(seg2->v/-qnet);
            tstep = MIN(dt, tstar);
            seg1->c = tankreact(seg1->c,seg1->v,Tank[i].Kb,tstep) +
                (qin*cin - qnet*(seg2->c) - qout*(seg1->c))*
                tstep/(seg1->v);
            seg2->c = tankreact(seg2->c,seg2->v,Tank[i].Kb,tstep);
            seg2->v += qnet*tstep;
            seg1->c = MAX(0.0, seg1->c);
            seg2->c = MAX(0.0, seg2->c);
            seg2->v = MAX(0.0, seg2->v);
            dt -= tstep;
        }

        /* Case where ambient zone empty & mixing zone draining */
        if (dt > 1)
        {
            seg1->c = tankreact(seg1->c,seg1->v,Tank[i].Kb,dt) +
                qin*dt*(cin-(seg1->c))/(seg1->v);
            seg1->v += qnet*dt;
            seg1->c = MAX(0.0, seg1->c);
            seg1->v = MAX(0.0, seg1->v);
            seg2->c = 0.0;
        }
    }

    /* Use quality of mixed compartment (seg1) to */
    /* represent quality of tank since this is where */
    /* outflow begins to flow from */
    Tank[i].C = seg1->c;
    simulation->C[n] = Tank[i].C;
}


void  tankmix3(int i, long dt)
/*
**----------------------------------------------------------
**   Input:   i = tank index
**            dt = current WQ time step
**   Output:  none
**   Purpose: First-In-First-Out (FIFO) tank model
**----------------------------------------------------------
*/
{
    int   k,n;
    float vin,vnet,vout,vseg;
    float cin,vsum,csum;
    Pseg  seg;

    k = Nlinks + i;
    if (LastSeg[k] == NULL || FirstSeg[k] == NULL) return;

    /* React contents of each compartment */
    if (Reactflag)
    {
        seg = FirstSeg[k];
        while (seg != NULL)
        {
            seg->c = tankreact(seg->c,seg->v,Tank[i].Kb,dt);
            seg = seg->prev;
        }
    }

    /* Find inflows & outflows */
    n = Tank[i].Node;
    vnet = simulation->D[n] * dt;
    vin = VolIn[n];
    vout = vin - vnet;
    if (vin > 0.0) cin = MassIn[n]/VolIn[n];
    else           cin = 0.0;
    simulation->V[i] += vnet;

    /* Withdraw flow from first segment */
    vsum = 0.0;
    csum = 0.0;
    while (vout > 0.0)
    {
        seg = FirstSeg[k];
        if (seg == NULL) break;
        vseg = seg->v;           /* Flow volume from leading seg */
        vseg = MIN(vseg,vout);
        if (seg == LastSeg[k]) vseg = vout;
        vsum += vseg;
        csum += (seg->c)*vseg;
        vout -= vseg;            /* Remaining flow volume */
        if (vout >= 0.0 && vseg >= seg->v)  /* Seg used up */
        {
            FirstSeg[k] = seg->prev;
            if (FirstSeg[k] == NULL) LastSeg[k] = NULL;
            seg->prev = FreeSeg;
            FreeSeg = seg;
        }
        else                /* Remaining volume in segment */
        {
            seg->v -= vseg;
        }
    }

    /* Use quality withdrawn from 1st segment */
    /* to represent overall quality of tank */
    if (vsum > 0.0) Tank[i].C = csum/vsum;
    else            Tank[i].C = FirstSeg[k]->c;
    simulation->C[n] = Tank[i].C;

    /* Add new last segment for new flow entering tank */
    if (vin > 0.0)
    {
        if ( (seg = LastSeg[k]) != NULL)
        {
            /* Quality is the same, so just add flow volume to last seg */
            if (ABS(seg->c - cin) < Ctol) seg->v += vin;

            /* Otherwise add a new seg to tank */
            else addseg(k,vin,cin);
        }

        /* If no segs left then add a new one. */
        else addseg(k,vin,cin);
    }
}


void  tankmix4(int i, long dt)
/*
**----------------------------------------------------------
**   Input:   i = tank index
**            dt = current WQ time step
**   Output:  none
**   Purpose: Last In-First Out (LIFO) tank model
**----------------------------------------------------------
*/
{
    int   k, n;
    float vin, vnet, cin, vsum, csum, vseg;
    Pseg  seg, tmpseg;

    k = Nlinks + i;
    if (LastSeg[k] == NULL || FirstSeg[k] == NULL) return;

    /* React contents of each compartment */
    if (Reactflag)
    {
        seg = LastSeg[k];
        while (seg != NULL)
        {
            seg->c = tankreact(seg->c,seg->v,Tank[i].Kb,dt);
            seg = seg->prev;
        }
    }

    /* Find inflows & outflows */
    n = Tank[i].Node;
    vnet = simulation->D[n] * dt;
    vin = VolIn[n];
    if (vin > 0.0) cin = MassIn[n]/VolIn[n];
    else           cin = 0.0;
    simulation->V[i] += vnet;
    Tank[i].C = LastSeg[k]->c;

    /* If tank filling, then create new last seg */
    if (vnet > 0.0)
    {
        if ( (seg = LastSeg[k]) != NULL)
        {
            /* Quality is the same, so just add flow volume to last seg */
            if (ABS(seg->c - cin) < Ctol) seg->v += vnet;

            /* Otherwise add a new last seg to tank */
            /* which points to old last seg */
            else
            {
                tmpseg = seg;
                LastSeg[k] = NULL;
                addseg(k,vnet,cin);
                LastSeg[k]->prev = tmpseg;
            }
        }

        /* If no segs left then add a new one. */
        else addseg(k,vnet,cin);

        /* Update reported tank quality */
        Tank[i].C = LastSeg[k]->c;
    }

    /* If net emptying then remove last segments until vnet consumed */
    else if (vnet < 0.0)
    {
        vsum = 0.0;
        csum = 0.0;
        vnet = -vnet;
        while (vnet > 0.0)
        {
            seg = LastSeg[k];
            if (seg == NULL) break;
            vseg = seg->v;
            vseg = MIN(vseg,vnet);
            if (seg == FirstSeg[k]) vseg = vnet;
            vsum += vseg;
            csum += (seg->c)*vseg;
            vnet -= vseg;
            if (vnet >= 0.0 && vseg >= seg->v)  /* Seg used up */
            {
                LastSeg[k] = seg->prev;
                if (LastSeg[k] == NULL) FirstSeg[k] = NULL;
                seg->prev = FreeSeg;
                FreeSeg = seg;
            }
            else                /* Remaining volume in segment */
            {
                seg->v -= vseg;
            }
        }
        /* Reported tank quality is mixture of flow released and any inflow */
        Tank[i].C = (csum + MassIn[n])/(vsum + vin);
    }
    simulation->C[n] = Tank[i].C;
}


float  sourcequal(Psource source)
/*
**--------------------------------------------------------------
**   Input:   j = source index
**   Output:  returns source WQ value
**   Purpose: determines source concentration in current time period
**--------------------------------------------------------------
*/
{
    int   i;
    long  k;
    float c;

    /* Get source concentration (or mass flow) in original units */
    c = source->C0;

    /* Convert mass flow rate from min. to sec. */
    /* and convert concen. from liters to cubic feet */
    if (source->Type == MASS) c /= 60.0;
    else c /= Ucf[QUALITY];

    /* Apply time pattern if assigned */
    i = source->Pat;
    if (i == 0) return(c);
    k = ((Qtime+Pstart)/Pstep) % (long)Pattern[i].Length;
    return(c*Pattern[i].F[k]);
}


float  avgqual(int k)
/*
**--------------------------------------------------------------
**   Input:   k = link index
**   Output:  returns WQ value
**   Purpose: computes average quality in link k
**--------------------------------------------------------------
*/
{
    float  vsum = 0.0,
        msum = 0.0;
    Pseg   seg;

    if (Qualflag == NONE) return(0.);
    seg = FirstSeg[k];
    while (seg != NULL)
    {
        vsum += seg->v;
        msum += (seg->c)*(seg->v);
        seg = seg->prev;
    }
    if (vsum > 0.0) return(msum/vsum);
    else return (simulation->C[Link[k].N1] + simulation->C[Link[k].N2]) / 2.;
}


void  ratecoeffs()
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  none
**   Purpose: determines wall reaction coeff. for each pipe
**--------------------------------------------------------------
*/
{
    int   k;
    float kw;

    for (k=1; k<=Nlinks; k++)
    {
        kw = Link[k].Kw;
        if (kw != 0.0) kw = piperate(k);
        Link[k].R = kw;
        simulation->R[k] = 0.0;
    }
}                         /* End of ratecoeffs */


float piperate(int k)
/*
**--------------------------------------------------------------
**   Input:   k = link index
**   Output:  returns reaction rate coeff. for 1st-order wall
**            reactions or mass transfer rate coeff. for 0-order
**            reactions
**   Purpose: finds wall reaction rate coeffs.
**--------------------------------------------------------------
*/
{
    float a,d,u,kf,kw,y,Re,Sh;

    d = Link[k].Diam;                    /* Pipe diameter, ft */

/* Ignore mass transfer if Schmidt No. is 0 */
    if (Sc == 0.0)
    {
        if (WallOrder == 0.0) return(BIG);
        else return(Link[k].Kw*(4.0/d)/Ucf[ELEV]);
    }

/* Compute Reynolds No. */
    a = PI*d*d/4.0;
    u = ABS (simulation->Q[k]) / a;
    Re = u*d/Viscos;

/* Compute Sherwood No. for stagnant flow  */
/* (mass transfer coeff. = Diffus./radius) */
    if (Re < 1.0) Sh = 2.0;

/* Compute Sherwood No. for turbulent flow */
/* using the Notter-Sleicher formula.      */
    else if (Re >= 2300.0)
        Sh = 0.0149*pow(Re,0.88)*pow(Sc,0.333);

/* Compute Sherwood No. for laminar flow */
/* using Graetz solution formula.        */
    else
    {
        y = d/Link[k].Len*Re*Sc;
        Sh = 3.65+0.0668*y/(1.0+0.04*pow(y,0.667));
    }

/* Compute mass transfer coeff. (in ft/sec) */
    kf = Sh*Diffus/d;

/* For zero-order reaction, return mass transfer coeff. */
    if (WallOrder == 0.0) return(kf);

/* For first-order reaction, return apparent wall coeff. */
    kw = Link[k].Kw/Ucf[ELEV];       /* Wall coeff, ft/sec */
    kw = (4.0/d)*kw*kf/(kf+ABS(kw)); /* Wall coeff, 1/sec  */
    return(kw);
}                         /* End of piperate */


float  pipereact(int k, float c, float v, long dt)
/*
**------------------------------------------------------------
**   Input:   k = link index
**            c = current WQ in segment
**            v = segment volume
**            dt = time step
**   Output:  returns new WQ value
**   Purpose: computes new quality in a pipe segment after
**            reaction occurs
**------------------------------------------------------------
*/
{
    float cnew, dc, dcbulk, dcwall, rbulk, rwall;

    /* For water age (hrs), update concentration by timestep */
    if (Qualflag == AGE) return(c+(float)dt/3600.0);

    /* Otherwise find bulk & wall reaction rates */
    rbulk = bulkrate(c,Link[k].Kb,BulkOrder)*Bucf;
    rwall = wallrate(c,Link[k].Diam,Link[k].Kw,Link[k].R);

    /* Find change in concentration over timestep */
    dcbulk = rbulk*(float)dt;
    dcwall = rwall*(float)dt;

    /* Update cumulative mass reacted */
    if (simulation->Htime >= Rstart)
    {
        Wbulk += ABS(dcbulk)*v;
        Wwall += ABS(dcwall)*v;
    }

    /* Update concentration */
    dc = dcbulk + dcwall;
    cnew = c + dc;
    cnew = MAX(0.0,cnew);
    return(cnew);
}


float  tankreact(float c, float v, float kb, long dt)
/*
**-------------------------------------------------------
**   Input:   c = current WQ in tank
**            v = tank volume
**            kb = reaction coeff.
**            dt = time step
**   Output:  returns new WQ value
**   Purpose: computes new quality in a tank after
**            reaction occurs
**-------------------------------------------------------
*/
{
    float cnew, dc, rbulk;

/*** Updated 9/7/00 ***/
    /* If no reaction then return current WQ */
    if (!Reactflag) return(c);

    /* For water age, update concentration by timestep */
    if (Qualflag == AGE) return(c + (float)dt/3600.0);

    /* Find bulk reaction rate */
    rbulk = bulkrate(c,kb,TankOrder)*Tucf;

    /* Find concentration change & update quality */
    dc = rbulk*(float)dt;
    if (simulation->Htime >= Rstart) Wtank += ABS(dc)*v;
    cnew = c + dc;
    cnew = MAX(0.0,cnew);
    return(cnew);
}


float  bulkrate(float c, float kb, float order)
/*
**-----------------------------------------------------------
**   Input:   c = current WQ concentration
**            kb = bulk reaction coeff.
**            order = bulk reaction order
**   Output:  returns bulk reaction rate
**   Purpose: computes bulk reaction rate (mass/volume/time)
**-----------------------------------------------------------
*/
{
    float c1;

    /* Find bulk reaction potential taking into account */
    /* limiting potential & reaction order. */

    /* Zero-order kinetics: */
    if (order == 0.0) c = 1.0;

    /* Michaelis-Menton kinetics: */
    else if (order < 0.0)
    {
        c1 = Climit + SGN(kb)*c;
        if (ABS(c1) < TINY) c1 = SGN(c1)*TINY;
        c = c/c1;
    }

    /* N-th order kinetics: */
    else
    {
        /* Account for limiting potential */
        if (Climit == 0.0) c1 = c;
        else c1 = MAX(0.0, SGN(kb)*(Climit-c));

        /* Compute concentration potential */
        if (order == 1.0) c = c1;
        else if (order == 2.0) c = c1*c;
        else c = c1*pow(MAX(0.0,c),order-1.0);
    }

    /* Reaction rate = bulk coeff. * potential) */
    if (c < 0) c = 0;
    return(kb*c);
}


float  wallrate(float c, float d, float kw, float kf)
/*
**------------------------------------------------------------
**   Input:   c = current WQ concentration
**            d = pipe diameter
**            kw = intrinsic wall reaction coeff.
**            kf = mass transfer coeff. for 0-order reaction
**                 (ft/sec) or apparent wall reaction coeff.
**                 for 1-st order reaction (1/sec)
**   Output:  returns wall reaction rate in mass/ft3/sec
**   Purpose: computes wall reaction rate
**------------------------------------------------------------
*/
{
    if (kw == 0.0 || d == 0.0) return(0.0);
    if (WallOrder == 0.0)       /* 0-order reaction */
    {
        kf = SGN(kw)*c*kf;       /* Mass transfer rate (mass/ft2/sec)*/
        kw = kw*SQR(Ucf[ELEV]);  /* Reaction rate (mass/ft2/sec) */
        if (ABS(kf) < ABS(kw))   /* Reaction mass transfer limited */
            kw = kf;
        return(kw*4.0/d);        /* Reaction rate (mass/ft3/sec) */
    }
    else return(c*kf);          /* 1st-order reaction */
}

/************************* End of QUALITY.C ***************************/
#endif
