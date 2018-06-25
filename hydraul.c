/*
*********************************************************************

HYDRAUL.C --  Hydraulic Simulator for EPANET Program

AUTHOR:     M. Lopez-Ibanez (M. L\'opez-Ib\'a\~nez)
            Napier University, Edinburgh, UK.

$LastChangedBy: manu $ $Revision: 172 $
$LastChangedDate: 2008-09-10 15:23:48 +0200 (Wed, 10 Sep 2008) $

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
DATE:       6/5/00
            9/7/00
            10/25/00
            12/29/00
            3/1/01
            11/19/01
            6/24/02
AUTHOR:     L. Rossman
            US EPA - NRMRL

  This module contains the network hydraulic simulator.
  It simulates the network's hydraulic behavior over an
  an extended period of time and writes its results to the
  binary file HydFile.

  The entry points for this module are:
     openhyd()    -- called from ENopenH() in EPANET.C
     inithyd()    -- called from ENinitH() in EPANET.C
     runhyd()     -- called from ENrunH() in EPANET.C
     nexthyd()    -- called from ENnextH() in EPANET.C
     closehyd()   -- called from ENcloseH() in EPANET.C
     tankvolume() -- called from ENsetnodevalue() in EPANET.C
     setlinkstatus(),
     setlinksetting(),
     resistance()-- all called from ENsetlinkvalue() in EPANET.C

  External functions called by this module are:
     createsparse() -- see SMATRIX.C
     freesparse()   -- see SMATRIX.C
     linsolve()     -- see SMATRIX.C
     checkrules()   -- see RULES.C
     interp()       -- see EPANET.C
     savehyd()      -- see OUTPUT.C
     savehydstep()  -- see OUTPUT.C
     writehydstat() -- see REPORT.C
     writehyderr()  -- see REPORT.C
     writehydwarn() -- see REPORT.C
*******************************************************************
*/

#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include "toolkit.h"
#include "hash.h"
#include "text.h"
#include "types.h"
#include "funcs.h"
#define  EXTERN  extern
#include "vars.h"

#define   QZERO  1.e-6  /* Equivalent to zero flow */
#define   CBIG   1.e8   /* Big coefficient         */
#define   INVCBIG   1.e-8   /* 1.0 / Big coefficient         */
#define   CSMALL 1.e-6  /* Small coefficient       */

/* Constants used for computing Darcy-Weisbach friction factor */
#define A1  0.314159265359e04  /* 1000*PI */
#define A2  0.157079632679e04  /* 500*PI  */
#define A3  0.502654824574e02  /* 16*PI   */
#define A4  6.283185307        /* 2*PI    */
#define A8  4.61841319859      /* 5.74*(PI/4)^.9 */
#define A9  -8.685889638e-01   /* -2/ln(10)      */
#define AA  -1.5634601348      /* -2*.9*2/ln(10) */
#define AB  3.28895476345e-03  /* 5.74/(4000^.9) */
#define AC  -5.14214965799e-03 /* AA*AB */

#if __GNUC__ >= 3
# define inline		inline __attribute__ ((always_inline))
#else
# define inline		/* no inline */
#endif

/* Functions defined and called only in this file.  */
static int  allocmatrix(ENsimulation_t); /* Allocates matrix coeffs.  */
static void freematrix(ENsimulation_t);  /* Frees matrix coeffs.       */
static void pumpswitch(ENsimulation_t, int index, char status);
static inline void    
pipecoeff(ENsimulation_t, int, REAL);    /* Computes pipe coeff.       */
static void newcoeffs(ENsimulation_t sim);  /* Computes matrix coeffs.    */

static REAL DWcoeff(char type, float diam, float flow, 
                    float roughness, REAL *dfdq); /* Computes D-W coeff.   */
static float
initlinkflow(int i, char s, float k);  /* Initializes link flow      */

static void    addenergy(ENsimulation_t, long);/* Accumulates energy usage   */
static void    setlinkflow(ENsimulation_t, int,/* Sets link flow via headloss*/
                           float);
static void    demands(ENsimulation_t);     /* Computes current demands   */
static int     controls(ENsimulation_t);    /* Controls link settings     */
static long    timestep(ENsimulation_t);    /* Computes new time step     */
static void
tanktimestep(ENsimulation_t, long *);      /* Time till tanks fill/drain */
static void  
controltimestep(ENsimulation_t, long *);  /* Time till control action   */
static void  
ruletimestep(ENsimulation_t, long *);     /* Time till rule action      */
static void
tanklevels(ENsimulation_t, long);         /* Computes new tank levels   */
static float   tankgrade(int,float);      /* Finds tank grade from vol. */
static int     pswitch(ENsimulation_t);   /* Pressure switch controls   */
static int     netsolve(ENsimulation_t,   /* Solves network equations   */
                        int *,float *);
static int     badvalve(ENsimulation_t,   /* Checks for bad valve       */
                        int);           
static int    valvestatus(ENsimulation_t);/* Updates valve status       */
static int    linkstatus(ENsimulation_t); /* Updates link status        */
static char   cvstatus(char,float,float); /* Updates CV status          */
static char   pumpstatus(ENsimulation_t,  /* Updates pump status        */
                         int,float);
static char    prvstatus(ENsimulation_t,  /* Updates PRV status         */
                         int,char,float, float,float);
static char    psvstatus(ENsimulation_t,  /* Updates PSV status         */
                         int,char,float,float,float);
static char    fcvstatus(ENsimulation_t,  /* Updates FCV status         */
                         int,char,float,float);
static void    tankstatus(ENsimulation_t, /* Checks if tank full/empty  */
                          int,int,int);
static REAL    newflows(ENsimulation_t);  /* Updates link flows         */
static void    linkcoeffs(ENsimulation_t);/* Computes link coeffs.      */
static void    nodecoeffs(ENsimulation_t);/* Computes node coeffs.      */
static void    valvecoeffs(ENsimulation_t);/* Computes valve coeffs.     */
static void    pumpcoeff(ENsimulation_t, int); /* Computes pump coeff.       */
/*** Updated 10/25/00 ***/
/*** Updated 12/29/00 ***/
static void    curvecoeff(                    /* Computes curve coeffs.     */
                          int,float,float *, float *);
static void    gpvcoeff(ENsimulation_t, int); /* Computes GPV coeff.        */
static void    pbvcoeff(ENsimulation_t, int); /* Computes PBV coeff.        */
static void    tcvcoeff(ENsimulation_t, int); /* Computes TCV coeff.        */
static void    prvcoeff(ENsimulation_t,       /* Computes PRV coeff.        */
                        int,int,int);         
static void    psvcoeff(ENsimulation_t,       /* Computes PSV coeff.        */
                        int,int,int);         
static void    fcvcoeff(ENsimulation_t,       /* Computes FCV coeff.        */
                        int,int,int);         
static void    emittercoeffs(ENsimulation_t); /* Computes emitter coeffs.   */
static REAL    emitflowchange(ENsimulation_t, /* Computes new emitter flow  */
                              int); 


ENsimulation_t openhyd()
/*
 *--------------------------------------------------------------
 *  Input:   none
 *  Output:  returns error code
 *  Purpose: opens hydraulics solver system
 *--------------------------------------------------------------
 */
{
    int  i;
    int  n;
    int  errcode = 0;
    ENsimulation_t sim;

    sim = malloc (sizeof(struct _ENsimulation_t));

    sim->smatrix = createsparse();     /* See SMATRIX.C  */

    sim->S    = calloc (MaxLinks + 1, sizeof(char));
    ERRCODE(MEMCHECK(sim->S));

    sim->OldStat = calloc (Nlinks + Ntanks + 1, sizeof(char));
    ERRCODE(MEMCHECK(sim->OldStat));

/* Allocate memory for network nodes */
    n = MaxNodes + 1;
    sim->H    = calloc (n, sizeof(REAL));
    sim->D    = calloc (n, sizeof(float));
    
    ERRCODE(MEMCHECK(sim->H));
    ERRCODE(MEMCHECK(sim->D));

    sim->E   = calloc (Nnodes + 1, sizeof(float));
    ERRCODE(MEMCHECK(sim->E));

    sim->V   = malloc ((MaxTanks + 1) * sizeof(float));
    ERRCODE(MEMCHECK(sim->V));

    sim->Pump = malloc ((MaxPumps + 1) * sizeof(SimPump_t));

/* Allocate memory for network links */
    if (!errcode)
    {
        n = MaxLinks + 1;
        sim->K    = calloc (n, sizeof(float));
        sim->Q    = calloc (n, sizeof(float));
        ERRCODE(MEMCHECK(sim->K));
        ERRCODE(MEMCHECK(sim->Q));
    }

    sim->Ncontrols = MaxControls;
    sim->Control = malloc ((MaxControls + 1) * sizeof(Scontrol));
    memcpy (sim->Control, InpControl, 
            (MaxControls + 1) * sizeof(Scontrol));
    
    sim->Npats     = MaxPats;
    sim->Pattern = malloc ((MaxPats + 1) * sizeof(Spattern));
    memcpy (sim->Pattern, InpPattern, (MaxPats + 1) * sizeof(Scontrol));
    for (i = 0; i <= sim->Npats; i++)
    {
        sim->Pattern[i].F = malloc (InpPattern[i].Length 
                                           * sizeof(float));
        memcpy (sim->Pattern[i].F, InpPattern[i].F, 
                InpPattern[i].Length * sizeof(float));        
    }
    
    ERRCODE(allocmatrix (sim));      /* Allocate solution matrices */
    for (i=1; i<=Nlinks; i++)    /* Initialize flows */
        sim->Q[i] = initlinkflow (i, Link[i].Stat, Link[i].Kc);

    return sim;
}


/*** Updated 3/1/01 ***/
void inithyd(ENsimulation_t sim, int initflag)
/*
**--------------------------------------------------------------
**  Input:   initflag > 0 if link flows should be re-initialized
**                    = 0 if not
**  Output:  none
**  Purpose: initializes hydraulics solver system
**--------------------------------------------------------------
*/
{
    int i,j,s;

    /* Initialize tanks */
    for (i=1; i<=Ntanks; i++)
    {
        sim->V[i] = Tank[i].V0;
        sim->H[Tank[i].Node] = Tank[i].H0;
        sim->D[Tank[i].Node] = 0.0;
        sim->OldStat[Nlinks+i] = TEMPCLOSED;
    }


    /* Initialize emitter flows */
    memset (sim->E, 0, (Nnodes+1) * sizeof(float));
    for (i=1; i<=Njuncs; i++)
        if (Node[i].Ke > 0.0) sim->E[i] = 1.0;

    /* Initialize links */
    for (i=1; i<=Nlinks; i++)
    {
        /* Initialize status and setting */
        sim->S[i] = Link[i].Stat;
        sim->K[i] = Link[i].Kc;

        /* Start active PRVs & PSVs in OPEN position */
        if ((Link[i].Type == PRV || Link[i].Type == PSV)
            && (Link[i].Kc != MISSING)) 
            sim->S[i] = OPEN;

/*** Updated 3/1/01 ***/
        /* Initialize flows if necessary */
        if (sim->S[i] <= CLOSED) sim->Q[i] = QZERO;
        else if (ABS (sim->Q[i]) <= QZERO || initflag > 0)
            sim->Q[i] = initlinkflow (i, sim->S[i], sim->K[i]);

        /* Save initial status */
        sim->OldStat[i] = sim->S[i];
    }

    /* Reset pump energy usage
       and pump switches       */
    for (i=1; i<=Npumps; i++)
    {
        for (j=0; j<6; j++) sim->Pump[i].Energy[j] = 0.0;
        sim->Pump[i].Nswitches = 0;
        sim->Pump[i].Schedule_idx = 0;
        for(s=0; s < 24; s++) sim->Pump[i].Schedule[s] = -1;
    }
    sim->Emax = 0.0;             /* Zero peak energy usage         */

#if 0
    /* Re-position hydraulics file */
    if (sim->Saveflag) fseek(HydFile,HydOffset,SEEK_SET);
#endif
/*** Updated 3/1/01 ***/
    /* Initialize current time */
    sim->Haltflag = 0;
    sim->Htime = 0;
    sim->Rtime = Rstep;
    checkrules (sim, 0);
}


int   runhyd(ENsimulation_t sim, long *t)
/*
**--------------------------------------------------------------
**  Input:   none
**  Output:  t = pointer to current time (in seconds)
**  Returns: error code
**  Purpose: solves network hydraulics in a single time period
**--------------------------------------------------------------
*/
{
    int   iter;                          /* Iteration count   */
    int   errcode;                       /* Error code        */
    float relerr;                        /* Solution accuracy */

    /* Find new demands & control actions */
    *t = sim->Htime;
    demands (sim);
    controls (sim);

    /* Solve network hydraulic equations */
    errcode = netsolve(sim, &iter,&relerr);
    if (!errcode)
    {
#if 0
        /* Report new status & save results */
        if (Statflag) writehydstat(iter,relerr);
#endif
/*** Updated 3/1/01 ***/
        /* If system unbalanced and no extra trials */
        /* allowed, then activate the Haltflag.     */
        if (relerr > Hacc && ExtraIter == -1) sim->Haltflag = 1;

        /* Report any warning conditions */
        if (!errcode) errcode = writehydwarn (sim, iter, relerr);
    }
    return(errcode);
}                               /* end of runhyd */


int  nexthyd(ENsimulation_t sim, long *tstep)
/*
**--------------------------------------------------------------
**  Input:   none
**  Output:  tstep = pointer to time step (in seconds)
**  Returns: error code
**  Purpose: finds length of next time step & updates tank
**           levels and rule-based contol actions
**--------------------------------------------------------------
*/
{
    long  hydstep;         /* Actual time step  */
    int   errcode = 0;     /* Error code        */
    int   i;
    double tmp_float;

/*** Updated 3/1/01 ***/
    /* Save current results to hydraulics file and */
    /* force end of simulation if Haltflag is active */
#if 0
    if (sim->Saveflag) errcode = savehyd (&sim->Htime);
#endif
    if (sim->Haltflag) sim->Htime = Dur;

    /* Compute next time step & update tank levels */
    *tstep = 0;
    hydstep = 0;
    if (sim->Htime < Dur) hydstep = timestep (sim);
#if 0
    if (sim->Saveflag) errcode = savehydstep(&hydstep);
#endif

    /* Compute pumping energy */
    if (Dur == 0) addenergy (sim, 0);
    else if (sim->Htime < Dur) addenergy (sim, hydstep);

    if (sim->Htime < Dur) 
    {
        sim->TotalSystemDemand += sim->Dsystem * hydstep;

        for (i=1, tmp_float=0.0; i<=Njuncs; i++) 
            tmp_float += sim->E[i];
        sim->TotalSystemLeakage += tmp_float * hydstep;

        for (i=1, tmp_float=0.0; i<=Ntanks; i++)
            if (Tank[i].A == 0.0)
                tmp_float += sim->D[Tank[i].Node];

        sim->TotalSystemInflow += tmp_float * hydstep;

    /* Update current time. */
        sim->Htime += hydstep;
        if (sim->Htime >= sim->Rtime) 
            sim->Rtime += Rstep;
    }
    else
    {
        sim->Htime++;          /* Force completion of analysis */
    }
    *tstep = hydstep;
    return(errcode);
}


void  closehyd(ENsimulation_t sim)
/*
**--------------------------------------------------------------
**  Input:   none
**  Output:  returns error code
**  Purpose: closes hydraulics solver system
**--------------------------------------------------------------
*/
{
    int j;

    freesparse (sim->smatrix);           /* see SMATRIX.C */
    freematrix (sim);

    free(sim->Control);

    /* Free memory for time patterns */
    if (sim->Pattern != NULL)
    {
        for (j=0; j<= sim->Npats; j++) free (sim->Pattern[j].F);
        free(sim->Pattern);
    }

    /* Free memory for computed results */
    free (sim->Pump);
    free (sim->S);
    free (sim->OldStat);
    free (sim->H);
    free (sim->D);
    free (sim->E);
    free (sim->K);
    free (sim->Q);
    free (sim->V);

    free (sim);
}


static int  allocmatrix(ENsimulation_t sim)
/*
**--------------------------------------------------------------
**  Input:   none
**  Output:  returns error code
**  Purpose: allocates memory used for solution matrix coeffs.
**--------------------------------------------------------------
*/
{
    int errcode = 0;
    /*  Purpose: allocates memory used for solution matrix coeffs.  */
    sim->Aii = calloc (Nnodes + 1, sizeof(REAL));
    sim->Aij = calloc (smatrix_non_zero_coeffs (sim->smatrix) + 1, 
                              sizeof(REAL));
    sim->F   = calloc (Nnodes + 1, sizeof(REAL));
    sim->P   = calloc (Nlinks + 1, sizeof(float));
    sim->Y   = calloc (Nlinks + 1, sizeof(float));
    sim->X   = calloc (MAX ((Nnodes + 1), (Nlinks + 1)), sizeof(float));
    ERRCODE(MEMCHECK(sim->Aii));
    ERRCODE(MEMCHECK(sim->Aij));
    ERRCODE(MEMCHECK(sim->F));
    ERRCODE(MEMCHECK(sim->P));
    ERRCODE(MEMCHECK(sim->Y));
    ERRCODE(MEMCHECK(sim->X));
    return(errcode);
}                               /* end of allocmatrix */


static void  freematrix(ENsimulation_t sim)
/*
**--------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Purpose: frees memory used for solution matrix coeffs.
**--------------------------------------------------------------
*/
{
    free (sim->Aii);
    free (sim->Aij);
    free (sim->F);
    free (sim->P);
    free (sim->Y);
    free (sim->X);
}                               /* end of freematrix */


static float
initlinkflow(int i, char s, float k)
/*
**--------------------------------------------------------------------
**  Input:   i = link index
**           s = link status
**           k = link setting (i.e., pump speed)
**  Output:  none
**  Purpose: sets initial flow in link to QZERO if link is closed,
**           to design flow for a pump, or to flow at velocity of
**           1 fps for other links.
**--------------------------------------------------------------------
*/
{
    if (s == CLOSED) return QZERO;
    else if (Link[i].Type == PUMP) return k * Pump[PUMPINDEX(i)].Q0;
    else return PI * SQR(Link[i].Diam) / 4.0;
}


/*** Updated 9/7/00 ***/
static void  setlinkflow(ENsimulation_t sim, int k, float dh)
/*
**--------------------------------------------------------------
**  Input:   k = link index
**           dh = headloss across link
**  Output:  none
**  Purpose: sets flow in link based on current headloss
**--------------------------------------------------------------
*/
{
    int   i,p;
    float h0;
    REAL  x,y;

    switch (Link[k].Type)
    {
    case CV:
    case PIPE:

        /* For Darcy-Weisbach formula: */
        /* use approx. inverse of formula. */
        if (Formflag == DW)
        {
            x = -log(sim->K[k] / 3.7 / Link[k].Diam);
            y = sqrt(ABS(dh)/Link[k].R/1.32547);
            sim->Q[k] = x*y;
        }

        /* For Hazen-Williams or Manning formulas: */
        /* use inverse of formula. */
        else
        {
            /* Precision warning: Link[k].R is usually very small and
               both dh and Link[k].R are float while x is REAL.
               Therefore, there is some problem with the precision and
               both fabs() and fabsf() break compatibility.
            */
            x = ABS(dh)/Link[k].R;
            y = 1.0/Hexp;
            sim->Q[k] = pow(x,y);
        }

        /* Change sign of flow to match sign of headloss */
        if (dh < 0.0) sim->Q[k] = -sim->Q[k];
        break;

    case PUMP:

        /* Convert headloss to pump head gain */
        dh = -dh;
        p = PUMPINDEX(k);

        /* For custom pump curve, interpolate from curve */
        if (Pump[p].Ptype == CUSTOM)
        {
            dh = -dh * Ucf[HEAD] / SQR (sim->K[k]);
            i = Pump[p].Hcurve;
            sim->Q[k] = interp(Curve[i].Npts,Curve[i].Y,Curve[i].X, dh)
                * sim->K[k] / Ucf[FLOW];
        }

        /* Otherwise use inverse of power curve */
        else
        {
            h0 = -SQR (sim->K[k]) * Pump[p].H0;
            x = pow (sim->K[k], 2.0 - Pump[p].N);
            x = ABS(h0-dh)/(Pump[p].R*x);
            y = 1.0/Pump[p].N;
            sim->Q[k] = pow(x,y);
        }
        break;
    }
}


void  setlinkstatus(ENsimulation_t sim, int index, char value, 
                    char *s, float *k)
/*----------------------------------------------------------------
**  Input:   index  = link index
**           value  = 0 (CLOSED) or 1 (OPEN)
**           s      = pointer to link status
**           k      = pointer to link setting
**  Output:  none
**  Purpose: sets link status to OPEN or CLOSED
**----------------------------------------------------------------
*/
{
    /* Status set to open */
    if (value == 1)
    {
        /* Adjust link setting for pumps & valves */
        if (Link[index].Type == PUMP) *k = 1.0;

/*** Updated 9/7/00 ***/
        if (Link[index].Type >  PUMP
            &&  Link[index].Type != GPV) *k = MISSING;

        /* Reset link flow if it was originally closed */
        if (*s <= CLOSED) {
            if (Link[index].Type == PUMP) pumpswitch (sim, index, OPEN);
            sim->Q[index] = initlinkflow (index, OPEN, *k);
        }
        *s = OPEN;
    }

    /* Status set to closed */
    else if (value == 0)
    {
        /* Adjust link setting for pumps & valves */
        if (Link[index].Type == PUMP) *k = 0.0;

/*** Updated 9/7/00 ***/
        if (Link[index].Type >  PUMP
            &&  Link[index].Type != GPV) *k = MISSING;

        /* Reset link flow if it was originally open */
        if (*s > CLOSED) {
            if(Link[index].Type == PUMP) pumpswitch (sim, index, CLOSED);
            sim->Q[index] = initlinkflow(index, CLOSED, *k);
        }
        *s = CLOSED;
    }
}


void  setlinksetting(ENsimulation_t sim, int index, float value, 
                     char *status, float *k)
/*----------------------------------------------------------------
**  Input:   index  = link index
**           value  = pump speed or valve setting
**           s      = pointer to link status
**           k      = pointer to link setting
**  Output:  none
**  Purpose: sets pump speed or valve setting, adjusting link
**           status and flow when necessary
**----------------------------------------------------------------
*/
{
    /* For a pump, status is OPEN if speed > 0, CLOSED otherwise */
    if (Link[index].Type == PUMP)
    {
        *k = value;
        if (value > 0 && *status <= CLOSED)
        {
            *status = OPEN;
            pumpswitch (sim, index, OPEN);
            sim->Q[index] = initlinkflow (index, OPEN, value);
        }
        else if (value == 0 && *status > CLOSED)
        {
            *status = CLOSED;
            pumpswitch (sim, index, CLOSED);
            sim->Q[index] = initlinkflow (index, CLOSED, value);
        }
    }

/***  Updated 9/7/00  ***/
    /* For FCV, activate it */
    else if (Link[index].Type == FCV)
    {
        if (*status <= CLOSED) 
            sim->Q[index] = initlinkflow (index, OPEN, value);
        *k = value;
        *status = ACTIVE;
    }

    /* Open closed control valve with fixed status (setting = MISSING) */
    else
    {
        if (*k == MISSING && *status <= CLOSED)
        {
            sim->Q[index] = initlinkflow (index, OPEN, value);
            *status = OPEN;
        }
        *k = value;
    }
}


void  resistance(int k)
/*
**--------------------------------------------------------------------
**  Input:   k = link index
**  Output:  none
**  Purpose: computes link flow resistance
**--------------------------------------------------------------------
*/
{
    float e,d,L;
    Link[k].R = CSMALL;
    if (Link[k].Type == PIPE || Link[k].Type == CV)
        switch (Link[k].Type)
        {

            /* Link is a pipe. Compute resistance based on headloss formula. */
            /* Friction factor for D-W formula gets included during solution */
            /* process in pipecoeff() function.                              */
        case CV:
        case PIPE:
            e = Link[k].Kc;                 /* Roughness coeff. */
            d = Link[k].Diam;               /* Diameter */
            L = Link[k].Len;                /* Length */
            switch(Formflag)
            {
            case HW: Link[k].R = 4.727*L/pow(e,Hexp)/pow(d,4.871);
                break;
            case DW: Link[k].R = L/2.0/32.2/d/SQR(PI*SQR(d)/4.0);
                break;
            case CM: Link[k].R = SQR(4.0*e/(1.49*PI*d*d))*
                         pow((d/4.0),-1.333)*L;
            }
            break;

            /* Link is a pump. Use negligible resistance. */
        case PUMP:
            Link[k].R = CSMALL;
            break;

            /* Link is a valve. Compute resistance for open valve assuming  */
            /* length is 2*diameter and friction factor is 0.02. Use with   */
            /* other formulas as well since resistance should be negligible.*/
        default:
            d = Link[k].Diam;
            L = 2.0*d;
            Link[k].R = 0.02*L/2.0/32.2/d/SQR(PI*SQR(d)/4.0);
            break;
        }
}


static void  demands(ENsimulation_t sim)
/*
**--------------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Purpose: computes demands at nodes during current time period
**--------------------------------------------------------------------
*/
{
    int i,j,n;
    long k,p;
    float djunc, sum;
    Pdemand demand;

    /* Determine total elapsed number of pattern periods */
    p = (sim->Htime + Pstart) / Pstep;

    /* Update demand at each node according to its assigned pattern */
    sim->Dsystem = 0.0;          /* System-wide demand */
    for (i=1; i<=Njuncs; i++)
    {
        sum = 0.0;
        for (demand = Node[i].D; demand != NULL; demand = demand->next)
        {
            /*
              pattern period (k) = (elapsed periods) modulus
              (periods per pattern)
            */
            j = demand->Pat;
            k = p % (long) sim->Pattern[j].Length;
            djunc = (demand->Base) * sim->Pattern[j].F[k]*Dmult;
            if (djunc > 0.0) sim->Dsystem += djunc;
            sum += djunc;
        }
        sim->D[i] = sum;
    }

    /* Update head at fixed grade nodes with time patterns. */
    for (n=1; n<=Ntanks; n++)
    {
        if (Tank[n].A == 0.0)
        {
            j = Tank[n].Pat;
            if (j > 0)
            {
                k = p % (long) sim->Pattern[j].Length;
                i = Tank[n].Node;
                sim->H[i] = Node[i].El * sim->Pattern[j].F[k];
            }
        }
    }

    /* Update status of pumps with utilization patterns */
    for (n=1; n<=Npumps; n++)
    {
        j = Pump[n].Upat;
        if (j > 0)
        {
            i = Pump[n].Link;
            k = p % (long) sim->Pattern[j].Length;
            setlinksetting (sim, i, sim->Pattern[j].F[k], 
                            &(sim->S[i]), &(sim->K[i]));
        }
    }
}                        /* End of demands */


static int  controls(ENsimulation_t sim)
/*
**---------------------------------------------------------------------
**  Input:   none
**  Output:  number of links whose setting changes
**  Purpose: implements simple controls based on time or tank levels
**---------------------------------------------------------------------
*/
{
    int   i, k, n, reset, setsum;
    float h, vplus;
    float v1, v2;
    float k1, k2;
    char  s1, s2;

    /* Examine each control statement */
    setsum = 0;
    for (i = 1; i <= sim->Ncontrols; i++)
    {
        /* Make sure that link is defined */
        reset = 0;
        if ( (k = sim->Control[i].Link) <= 0) continue;

        /* Link is controlled by tank level */
        if ((n = sim->Control[i].Node) > 0 && n > Njuncs)
        {
            h = sim->H[n];
            vplus = ABS (sim->D[n]);
            v1 = tankvolume (n - Njuncs, h);
            v2 = tankvolume (n - Njuncs, sim->Control[i].Grade);
            if (sim->Control[i].Type == LOWLEVEL && v1 <= v2 + vplus)
                reset = 1;
            if (sim->Control[i].Type == HILEVEL && v1 >= v2 - vplus)
                reset = 1;
        }

        /* Link is time-controlled */
        if (sim->Control[i].Type == TIMER)
        {
            if (sim->Control[i].Time == sim->Htime) reset = 1;
        }

        /* Link is time-of-day controlled */
        if (sim->Control[i].Type == TIMEOFDAY)
        {
            if ((sim->Htime + Tstart) % SECperDAY == sim->Control[i].Time)
                reset = 1;
        }

        /* Update link status & pump speed or valve setting */
        if (reset == 1)
        {
            if (sim->S[k] <= CLOSED) s1 = CLOSED;
            else                     s1 = OPEN;
            s2 = sim->Control[i].Status;
            k1 = sim->K[k];
            k2 = k1;
            if (Link[k].Type > PIPE) k2 = sim->Control[i].Setting;
            if (s1 != s2 || k1 != k2)
            {
                sim->S[k] = s2;
                sim->K[k] = k2;
#if 0
                if (Statflag) writecontrolaction(k,i);
#endif
                if (s1 != s2) {
                    if (Link[k].Type == PUMP) pumpswitch (sim, k, s2);
                    sim->Q[k]  = initlinkflow (k, sim->S[k], sim->K[k]);
                }
                setsum++;
            }
        }
    }
    return(setsum);
}                        /* End of controls */


static long  timestep(ENsimulation_t sim)
/*
**----------------------------------------------------------------
**  Input:   none
**  Output:  returns time step until next change in hydraulics
**  Purpose: computes time step to advance hydraulic simulation
**----------------------------------------------------------------
*/
{
    long   n,t,tstep;

    /* Normal time step is hydraulic time step */
    tstep = Hstep;

    /* Revise time step based on time until next demand period */
    n = ((sim->Htime + Pstart) / Pstep) + 1;  /* Next pattern period  */
    t = n * Pstep - sim->Htime;              /* Time till next period */
    if (t > 0 && t < tstep) tstep = t;

    /* Revise time step based on time until next reporting period */
    t = sim->Rtime - sim->Htime;
    if (t > 0 && t < tstep) tstep = t;

    /* Revise time step based on smallest time to fill or drain a tank */
    tanktimestep (sim, &tstep);

    /* Revise time step based on smallest time to activate a control */
    controltimestep (sim, &tstep);

    /* Evaluate rule-based controls (which will also update tank levels) */
    if (Nrules > 0) ruletimestep (sim, &tstep);
    else tanklevels(sim, tstep);
    return(tstep);
}


static void  tanktimestep(ENsimulation_t sim, long *tstep)
/*
**-----------------------------------------------------------------
**  Input:   *tstep = current time step
**  Output:  *tstep = modified current time step
**  Purpose: revises time step based on shortest time to fill or
**           drain a tank
**-----------------------------------------------------------------
*/
{
    int    i,n;
    float  h,q,v;
    long   t;

    /* (D[n] is net flow rate into (+) or out of (-) tank at node n) */
    for (i=1; i<=Ntanks; i++)
    {
        if (Tank[i].A == 0.0) continue;           /* Skip reservoirs     */
        n = Tank[i].Node;
        h = sim->H[n];                     /* Current tank grade  */
        q = sim->D[n];                     /* Flow into tank      */
        if (ABS(q) <= QZERO) continue;
        if (q > 0.0 && h < Tank[i].Hmax)
        {
            v = Tank[i].Vmax - sim->V[i];  /* Volume to fill      */
        }
        else if (q < 0.0 && h > Tank[i].Hmin)
        {
            v = Tank[i].Vmin - sim->V[i];  /* Volume to drain (-) */
        }
        else continue;
        t = (long)ROUND(v/q);                     /* Time to fill/drain  */
        if (t > 0 && t < *tstep) *tstep = t;
    }
}


static void  controltimestep(ENsimulation_t sim, long *tstep)
/*
**------------------------------------------------------------------
**  Input:   *tstep = current time step
**  Output:  *tstep = modified current time step
**  Purpose: revises time step based on shortest time to activate
**           a simple control
**------------------------------------------------------------------
*/
{
    int   i,j,k,n;
    float h,q,v;
    long  t,t1,t2;

    for (i = 1; i <= sim->Ncontrols; i++)
    {
        if (sim->Control[i].Link <= 0) 
            continue; /* This control is removed so I do not look at it */

        t = 0;
        if ( (n = sim->Control[i].Node) > 0) /* Node control:       */
        {
            if ((j = n-Njuncs) <= 0) continue;     /* Node is a tank      */
            h = sim->H[n];                  /* Current tank grade  */
            q = sim->D[n];                  /* Flow into tank      */
            if (ABS(q) <= QZERO) continue;

            if ((h < sim->Control[i].Grade 
                 && sim->Control[i].Type == HILEVEL /* Tank below hi level */
                 && q > 0.0)                          /* & is filling        */
                || (h > sim->Control[i].Grade 
                    && sim->Control[i].Type == LOWLEVEL /* Tank above low level */ 
                    && q < 0.0))                     /* & is emptying        */
            {                                      /* Time to reach level  */
                v = tankvolume (j, sim->Control[i].Grade) - sim->V[j];
                t = (long)ROUND(v/q);
            }
        }

        if (sim->Control[i].Type == TIMER)    /* Time control:        */
        {
            if (sim->Control[i].Time > sim->Htime)
                t = sim->Control[i].Time - sim->Htime;
        }

        if (sim->Control[i].Type == TIMEOFDAY)/* Time-of-day control: */
        {
            t1 = (sim->Htime + Tstart) % SECperDAY;
            t2 = sim->Control[i].Time;
            if (t2 >= t1) t = t2 - t1;
            else t = SECperDAY - t1 + t2;
        }

        if (t > 0 && t < *tstep)               /* Revise time step     */
        {
            /* Check if rule actually changes link status or setting */
            k = sim->Control[i].Link;
            if ((Link[k].Type > PIPE && sim->K[k] != sim->Control[i].Setting) 
                || sim->S[k] != sim->Control[i].Status)
                *tstep = t;
        }
    }
}                        /* End of timestep */


static void  ruletimestep(ENsimulation_t sim, long *tstep)
/*
**--------------------------------------------------------------
**  Input:   *tstep = current time step (sec)
**  Output:  *tstep = modified time step
**  Purpose: updates next time step by checking if any rules
**           will fire before then; also updates tank levels.
**--------------------------------------------------------------
*/
{
    long tnow,      /* Start of time interval for rule evaluation */
        tmax,      /* End of time interval for rule evaluation   */
        dt,        /* Normal time increment for rule evaluation  */
        dt1;       /* Actual time increment for rule evaluation  */

    /* Find interval of time for rule evaluation */
    tnow = sim->Htime;
    tmax = tnow + *tstep;

    /* If no rules, then time increment equals current time step */
    if (Nrules == 0)
    {
        dt = *tstep;
        dt1 = dt;
    }

    /* Otherwise, time increment equals rule evaluation time step and */
    /* first actual increment equals time until next even multiple of */
    /* Rulestep occurs. */
    else
    {
        dt = Rulestep;
        dt1 = Rulestep - (tnow % Rulestep);
    }

    /* Make sure time increment is no larger than current time step */
    dt = MIN(dt, *tstep);
    dt1 = MIN(dt1, *tstep);
    if (dt1 == 0) dt1 = dt;

    /* Step through time, updating tank levels, until either  */
    /* a rule fires or we reach the end of evaluation period. */
    /*
    ** Note: we are updating the global simulation time (Htime)
    **       here because it is used by functions in RULES.C
    **       to evaluate rules when checkrules() is called.
    **       It is restored to its original value after the
    **       rule evaluation process is completed (see below).
    **       Also note that dt1 will equal dt after the first
    **       time increment is taken.
    */
    do
    {
        /* FIXME: try passing Htime to checkrules to avoid the reset
           of simulation->H time */
        sim->Htime += dt1;   /* Update simulation clock */
        tanklevels (sim, dt1);            /* Find new tank levels    */
        if (checkrules (sim, dt1)) break; /* Stop if rules fire      */
        dt = MIN(dt, tmax - sim->Htime); /* Update time increment   */
        dt1 = dt;                   /* Update actual increment */
    }  while (dt > 0);             /* Stop if no time left    */

    /* Compute an updated simulation time step (*tstep) */
    /* and return simulation time to its original value */
    *tstep = sim->Htime - tnow;
    sim->Htime = tnow;
}


static void  addenergy(ENsimulation_t sim, long hstep)
/*
**-------------------------------------------------------------
**  Input:   hstep = time step (sec)
**  Output:  none
**  Purpose: accumulates pump energy usage
**-------------------------------------------------------------
*/
{
    int   i,j,k;
    long  m,n;
    float c0,c,             /* Energy cost (cost/kwh) */
        f0,               /* Energy cost factor */
        dt,               /* Time interval (hr) */
        e,                /* Pump efficiency (fraction) */
        q,                /* Pump flow (cfs) */
        p,                /* Pump energy (kw) */
        psum = 0.0;       /* Total energy (kw) */

    /* Determine current time interval in hours */
    if      (Dur == 0)    dt = 1.0;
    else if (sim->Htime < Dur) dt = (float) hstep / 3600.0;
    else                  dt = 0.0;
    if      (dt == 0.0)   return;
    n = (sim->Htime + Pstart) / Pstep;

    /* Compute default energy cost at current time */
    c0 = Ecost;
    f0 = 1.0;
    if (Epat > 0)
    {
        m = n % (long) sim->Pattern[Epat].Length;
        f0 = sim->Pattern[Epat].F[m];
    }

    /* Examine each pump */
    for (j=1; j<=Npumps; j++)
    {
        /* Skip closed pumps */
        k = Pump[j].Link;
        if (sim->S[k] <= CLOSED) continue;
        q = MAX(QZERO, fabs (sim->Q[k]));

        /* Find pump-specific energy cost */
        if (Pump[j].Ecost > 0.0) c = Pump[j].Ecost;
        else c = c0;

        if ( (i = Pump[j].Epat) > 0)
        {
            m = n % (long) sim->Pattern[i].Length;
            c *= sim->Pattern[i].F[m];
        }
        else c *= f0;

        /* Find pump energy & efficiency */
        getenergy (sim, k, &p, &e);
        psum += p;

        /* Update pump's cumulative statistics */
        sim->Pump[j].Energy[0] += dt;            /* Time on-line */
        sim->Pump[j].Energy[1] += e*dt;          /* Effic.-hrs   */
        sim->Pump[j].Energy[2] += p/q*dt;        /* kw/cfs-hrs   */
        sim->Pump[j].Energy[3] += p*dt;          /* kw-hrs       */
        sim->Pump[j].Energy[4] = MAX (sim->Pump[j].Energy[4],p);
        sim->Pump[j].Energy[5] += c*p*dt;        /* cost-hrs.    */
    }

    /* Update maximum kw value */
    sim->Emax = MAX (sim->Emax, psum);
}                       /* End of pumpenergy */


void  getenergy(ENsimulation_t sim, int k, float *kw, float *eff)
/*
**----------------------------------------------------------------
**  Input:   k    = link index
**  Output:  *kw  = kwatt energy used
**           *eff = efficiency (pumps only)
**  Purpose: computes flow energy associated with link k
**----------------------------------------------------------------
*/
{
    int   i,j;
    float dh, q, e;

/*** Updated 6/24/02 ***/
    /* No energy if link is closed */
    if (sim->S[k] <= CLOSED)
    {
        *kw = 0.0;
        *eff = 0.0;
        return;
    }
/*** End of update ***/

    /* Determine flow and head difference */
    q = fabs (sim->Q[k]);
    dh = ABS (sim->H[Link[k].N1] - sim->H[Link[k].N2]);

    /* For pumps, find effic. at current flow */
    if (Link[k].Type == PUMP)
    {
        j = PUMPINDEX(k);
        e = Epump;
        if ( (i = Pump[j].Ecurve) > 0)
            e = interp(Curve[i].Npts,Curve[i].X,Curve[i].Y,q*Ucf[FLOW]);
        e = MIN(e, 100.0);
        e = MAX(e, 1.0);
        e /= 100.0;
    }
    else e = 1.0;

    /* Compute energy */
    *kw = dh*q*SpGrav/8.814/e*KWperHP;
    *eff = e;
}


static void  tanklevels(ENsimulation_t sim, long tstep)
/*
**----------------------------------------------------------------
**  Input:   tstep = current time step
**  Output:  none
**  Purpose: computes new water levels in tanks after current
**           time step
**----------------------------------------------------------------
*/
{
    int   i,n;
    float dv;

    for (i=1; i<=Ntanks; i++)
    {

        /* Skip reservoirs */
        if (Tank[i].A == 0.0) continue;

        /* Update the tank's volume & water elevation */
        n = Tank[i].Node;
        dv = sim->D[n] * tstep;
        sim->V[i] += dv;

        /*** Updated 6/24/02 ***/
        /* Check if tank full/empty within next second */
        if (sim->V[i] + sim->D[n] >= Tank[i].Vmax) 
            sim->V[i] = Tank[i].Vmax;
        if (sim->V[i] - sim->D[n] <= Tank[i].Vmin) 
            sim->V[i] = Tank[i].Vmin;

        sim->H[n] = tankgrade (i, sim->V[i]);
    }
}                       /* End of tanklevels */


float  tankvolume(int i, float h)
/*
**--------------------------------------------------------------------
**  Input:   i = tank index
**           h = water elevation in tank
**  Output:  returns water volume in tank
**  Purpose: finds water volume in tank i corresponding to elev. h.
**--------------------------------------------------------------------
*/
{
    int j;

    /* Use level*area if no volume curve */
    j = Tank[i].Vcurve;
    if (j == 0) return(Tank[i].Vmin + (h - Tank[i].Hmin)*Tank[i].A);

    /* If curve exists, interpolate on h to find volume v */
    /* remembering that volume curve is in original units.*/
    else return(interp(Curve[j].Npts, Curve[j].X, Curve[j].Y,
                       (h-Node[Tank[i].Node].El)*Ucf[HEAD])/Ucf[VOLUME]);

}                       /* End of tankvolume */


static float  tankgrade(int i, float v)
/*
**-------------------------------------------------------------------
**  Input:   i = tank index
**           v = volume in tank
**  Output:  returns water level in tank
**  Purpose: finds water level in tank i corresponding to volume v.
**-------------------------------------------------------------------
*/
{
    int j;

    /* Use area if no volume curve */
    j = Tank[i].Vcurve;
    if (j == 0) return(Tank[i].Hmin + (v - Tank[i].Vmin)/Tank[i].A);

    /* If curve exists, interpolate on volume (originally the Y-variable */
    /* but used here as the X-variable) to find new level above bottom.  */
    /* Remember that volume curve is stored in original units.           */
    else return(Node[Tank[i].Node].El +
                interp(Curve[j].Npts, Curve[j].Y, Curve[j].X, v*Ucf[VOLUME])/Ucf[HEAD]);

}                        /* End of tankgrade */


static int  netsolve(ENsimulation_t sim, int *iter, float *relerr)
/*
**-------------------------------------------------------------------
**  Input:   none
**  Output:  *iter   = # of iterations to reach solution
**           *relerr = convergence error in solution
**           returns error code
**  Purpose: solves network nodal equations for heads and flows
**           using Todini's Gradient algorithm
**
*** Updated 9/7/00 ***
**  Notes:   Status checks on PRVs and PSVs made every iteration.
**           Checks on other valves, pumps, and pipes to tanks
**           are made every CheckFreq iteration, up until
**           MaxCheck iterations are reached. After this, a
**           status check is made when convergence is achieved.
**           If convergence acheived but status changes, then
**           at least two more iterations are required before
**           solution is reached, even if system converges.
**           If convergence not achieved in MaxIter and ExtraIter
**           > 0 then another ExtraIter trials are made with no
**           status changes made to any links.
**
**   This procedure calls linsolve() which appears in SMATRIX.C.
**-------------------------------------------------------------------
*/
{
    int   i;                     /* Node index */
    int   errcode;               /* Node causing solution error */
    int   nextcheck;             /* Next status check trial */
    int   concheck = 1;          /* Convergence checking flag */
    int   maxtrials;             /* Max. trials for convergence */
    REAL  newerr;                /* New convergence error */

/*** Updated 9/7/00 ***/
    int   valvechange;           /* Valve status change flag */

/*** Updated 9/7/00 ***/
    nextcheck = CheckFreq;

    /* Repeat iterations until convergence or trial limit is exceeded. */
    /* (ExtraIter used to increase trials in case of status cycling.)  */
#if 0
    if (Statflag == FULL) writerelerr(0,0);
#endif
    maxtrials = MaxIter;
    if (ExtraIter > 0) maxtrials += ExtraIter;
    *iter = 1;
    do {
        /*
        ** Compute coefficient matrices A & F and solve A*H = F
        ** where H = heads, A = Jacobian coeffs. derived from
        ** head loss gradients, & F = flow correction terms.
        ** Solution for H is returned in F from call to linsolve().
        */
        newcoeffs (sim);
        errcode = linsolve (sim->smatrix, Njuncs, sim->Aii, sim->Aij, sim->F);

        /* Take action depending on error code */
        if (errcode < 0) break;    /* Memory allocation problem */
        if (errcode > 0)           /* Ill-conditioning problem */
        {
            /* If control valve causing problem, fix its status & continue, */
            /* otherwise end the iterations with no solution.               */
            if (badvalve(sim, errcode)) continue;
            else break;
        }

        /* Update current solution. */
        /* (Row[i] = row of solution matrix corresponding to node i). */
        for (i=1; i<=Njuncs; i++) 
            sim->H[i] = sim->F[sim->smatrix->Row[i]];   /* Update heads */
        
        newerr = newflows (sim);                        /* Update flows */
        *relerr = newerr;
#if 0
        /* Write convergence error to status report if called for */
        if (Statflag == FULL) writerelerr(*iter,*relerr);
#endif
/*** Updated 9/7/00 ***/
        /* Check if a PRV or PSV changes status */
        valvechange = valvestatus (sim);

        /* Check for convergence (concheck equals 0 if previous iteration */
        /* converged but there was a status change somewhere).            */
        if (concheck == 0) concheck = 1;
        else if (*relerr <= Hacc)
        {
            /* We have convergence. Quit if we are into extra iterations. */
            if (*iter > MaxIter) break;

/*** Updated 9/7/00 ***/
            /* Quit if no status changes occur. */
            if (!valvechange
                &&  !linkstatus (sim)
                &&  !pswitch (sim)) break;

            /* We have a status change so force another iteration. */
            concheck = 0;
            nextcheck = *iter + CheckFreq;
        }

        /* See if its time for a periodic status check; update iteration count. */
        if (*iter <= MaxCheck && *iter == nextcheck)
        {

/*** Updated 9/7/00 **/
            /*valvestatus();*/

            linkstatus (sim);
            nextcheck += CheckFreq;
        }
        (*iter)++;
    } while (*iter <= maxtrials);

    /* Iterations ended. Report any errors. */
    if (errcode < 0) errcode = 101;      /* Memory allocation error */
    else if (errcode > 0)
    {
        writehyderr(sim, errcode);      /* Ill-conditioned eqns. error */
        errcode = 110;
    }

    /* Add any emitter flows to junction demands */
    for (i=1; i<=Njuncs; i++) 
        sim->D[i] += sim->E[i];
    return(errcode);
}                        /* End of netsolve */


static int  badvalve(ENsimulation_t sim, int n)
/*
**-----------------------------------------------------------------
**  Input:   n = node index
**  Output:  returns 1 if node n belongs to an active control valve,
**           0 otherwise
**  Purpose: determines if a node belongs to an active control valve
**           whose setting causes an inconsistent set of eqns. If so,
**           the valve status is fixed open and a warning condition
**           is generated.
**-----------------------------------------------------------------
*/
{
    int i,k,n1,n2;
    for (i=1; i<=Nvalves; i++)
    {
        k = Valve[i].Link;
        n1 = Link[k].N1;
        n2 = Link[k].N2;
        if (n == n1 || n == n2)
        {
            if (Link[k].Type == PRV ||
                Link[k].Type == PSV ||
                Link[k].Type == FCV)
            {
                if (sim->S[k] == ACTIVE)
                {
#if 0
                    if (Statflag == FULL)
                    {
                        ENclockstring_t atime;
                        sprintf(Msg,FMT61, clocktime (&atime, sim->Htime),Link[k].ID);
                        writeline(Msg);
                    }
#endif
                    if (Link[k].Type == FCV) sim->S[k] = XFCV;
                    else                     sim->S[k] = XPRESSURE;
                    return(1);
                }
            }
            return(0);
        }
    }
    return(0);
}


static int  valvestatus(ENsimulation_t sim)
/*
**-----------------------------------------------------------------
**  Input:   none
**  Output:  returns 1 if any pressure control valve changes status,
**           0 otherwise
**  Purpose: updates status for PRVs & PSVs whose status
**           is not fixed to OPEN/CLOSED
**-----------------------------------------------------------------
*/
{
    int   change = FALSE,            /* Status change flag      */
        i,k,                       /* Valve & link indexes    */
        n1,n2;                     /* Start & end nodes       */
    char  s;                         /* Valve status settings   */
    float hset;                      /* Valve head setting      */

    for (i=1; i<=Nvalves; i++)                   /* Examine each valve   */
    {
        k = Valve[i].Link;                        /* Link index of valve  */
        if (sim->K[k] == MISSING) continue;/* Valve status fixed   */
        n1 = Link[k].N1;                          /* Start & end nodes    */
        n2 = Link[k].N2;
        s  = sim->S[k];                    /* Save current status  */
        if (s != CLOSED                    /* No change if flow is */
            && ABS (sim->Q[k]) < Qtol)     /* negligible.          */
            continue;            
        switch (Link[k].Type)                     /* Evaluate new status: */
        {
        case PRV:  
            hset = Node[n2].El + sim->K[k];
            sim->S[k] = prvstatus (sim, k, s, hset, sim->H[n1], sim->H[n2]);
            break;
        case PSV:  
            hset = Node[n1].El + sim->K[k];
            sim->S[k] = psvstatus (sim, k, s, hset, sim->H[n1], sim->H[n2]);
            break;
        default: continue;
        }

/*** Updated 9/7/00 ***/
        /* Do not reset flow in valve if its status changes. */
        /* This strategy improves convergence. */

        /* Check for status change */
        if (s != sim->S[k])
        {
#if 0
            if (Statflag == FULL) writestatchange (k, s, sim->S[k]);
#endif
            change = TRUE;
        }
    }
    return(change);
}                       /* End of valvestatus() */


/*** Updated 9/7/00 ***/
static int  linkstatus(ENsimulation_t sim)
/*
**--------------------------------------------------------------
**  Input:   none
**  Output:  returns 1 if any link changes status, 0 otherwise
**  Purpose: determines new status for pumps, CVs, FCVs & pipes
**           to tanks.
**--------------------------------------------------------------
*/
{
    int   change = FALSE,             /* Status change flag      */
        k,                          /* Link index              */
        n1,                         /* Start node index        */
        n2;                         /* End node index          */
    float dh;                         /* Head difference         */
    char  status;                     /* Current status          */

    /* Examine each link */
    for (k=1; k<=Nlinks; k++)
    {
        n1 = Link[k].N1;
        n2 = Link[k].N2;
        dh = sim->H[n1] - sim->H[n2];

        /* Re-open temporarily closed links (status = XHEAD or TEMPCLOSED) */
        status = sim->S[k];
        if (status == XHEAD || status == TEMPCLOSED) sim->S[k] = OPEN;

        /* Check for status changes in CVs and pumps */
        if (Link[k].Type == CV) 
            sim->S[k] = cvstatus (sim->S[k], dh, sim->Q[k]);

        if (Link[k].Type == PUMP && sim->S[k] >= OPEN)
            sim->S[k] = pumpstatus (sim, k, -dh);

        /* Check for status changes in non-fixed FCVs */
        if (Link[k].Type == FCV && sim->K[k] != MISSING)
            sim->S[k] = fcvstatus (sim, k, status, sim->H[n1], sim->H[n2]);

        /* Check for flow into (out of) full (empty) tanks */
        if (n1 > Njuncs || n2 > Njuncs) tankstatus(sim, k,n1,n2);

        /* Note change in link status and revise link flow */
        if (status != sim->S[k])
        {
            change = TRUE;
#if 0
            if (Statflag == FULL) writestatchange(k,status, sim->S[k]);
#endif
            if (sim->S[k] <= CLOSED) sim->Q[k] = QZERO;
            else setlinkflow (sim, k, dh);
        }
    }
    return(change);
}                        /* End of linkstatus */


static char  cvstatus(char s, float dh, float q)
/*
**--------------------------------------------------
**  Input:   s  = current status
**           dh = headloss
**           q  = flow
**  Output:  returns new link status
**  Purpose: updates status of a check valve.
**--------------------------------------------------
*/
{
    /* Prevent reverse flow through CVs */
    if (ABS(dh) > Htol)
    {
        if (dh < -Htol)     return(CLOSED);
        else if (q < -Qtol) return(CLOSED);
        else                return(OPEN);
    }
    else
    {
        if (q < -Qtol) return(CLOSED);
        else           return(s);
    }
}


static char  pumpstatus(ENsimulation_t sim, int k, float dh)
/*
**--------------------------------------------------
**  Input:   k  = link index
**           dh = head gain
**  Output:  returns new pump status
**  Purpose: updates status of a pump.
**--------------------------------------------------
*/
{
    int   p;
    float hmax;

    /* Prevent reverse flow through pump */
    p = PUMPINDEX(k);
    if (Pump[p].Ptype == CONST_HP) hmax = BIG;
    else hmax = SQR (sim->K[k]) * Pump[p].Hmax;
    if (dh > hmax + Htol) return(XHEAD);

    /* Check if pump cannot deliver flow */
    if (sim->Q[k] > sim->K[k] * Pump[p].Qmax + Qtol) 
        return XFLOW;
    return(OPEN);
}


static char  prvstatus(ENsimulation_t sim, int k, char s, float hset, 
                       float h1, float h2)
/*
**-----------------------------------------------------------
**  Input:   k    = link index
**           s    = current status
**           hset = valve head setting
**           h1   = head at upstream node
**           h2   = head at downstream node
**  Output:  returns new valve status
**  Purpose: updates status of a pressure reducing valve.
**-----------------------------------------------------------
*/
{
    char  status;     /* New valve status */
    float hml;        /* Minor headloss   */

    status = s;
    if (sim->K[k] == MISSING)           /* Status fixed by user */
        return status;
    hml = Link[k].Km * SQR (sim->Q[k]); /* Head loss when open  */
    switch (s)
    {
    case ACTIVE:
        if (sim->Q[k] < -Qtol) status = CLOSED;
        else if (h1 < hset+hml-Htol) status = OPEN;
        else                         status = ACTIVE;
        break;
    case OPEN:
        if (sim->Q[k] < -Qtol) status = CLOSED;
        else if (h1 > hset+hml+Htol) status = ACTIVE;
        else                         status = OPEN;
        break;
    case CLOSED:
        if (h1 > h2+Htol
            &&  h1 < hset-Htol)          status = OPEN;
        else if (h1 > h2+Htol
                 &&  h2 < hset-Htol)     status = ACTIVE;
        else                         status = CLOSED;
        break;
    case XPRESSURE:
        if (sim->Q[k] < -Qtol)            status = CLOSED;
        break;
    }
    return(status);
}


static char  psvstatus(ENsimulation_t sim, int k, char s, float hset, 
                       float h1, float h2)
/*
**-----------------------------------------------------------
**  Input:   k    = link index
**           s    = current status
**           hset = valve head setting
**           h1   = head at upstream node
**           h2   = head at downstream node
**  Output:  returns new valve status
**  Purpose: updates status of a pressure sustaining valve.
**-----------------------------------------------------------
*/
{
    char  status;       /* New valve status */
    float hml;          /* Minor headloss   */

    status = s;
    if (sim->K[k] == MISSING)           /* Status fixed by user */
        return status;       
    hml = Link[k].Km * SQR (sim->Q[k]); /* Head loss when open  */
    switch (s)
    {
    case ACTIVE:
        if (sim->Q[k] < -Qtol)            status = CLOSED;
        else if (h2+hml > hset+Htol) status = OPEN;
        else                         status = ACTIVE;
        break;
    case OPEN:
        if (sim->Q[k] < -Qtol)            status = CLOSED;
        else if (h2+hml < hset-Htol) status = ACTIVE;
        else                         status = OPEN;
        break;
    case CLOSED:
        if (h1 > h2+Htol
            &&  h2 > hset+Htol)          status = OPEN;
        else if (h1 > h2+Htol
                 &&  h1 > hset+Htol)          status = ACTIVE;
        else                         status = CLOSED;
        break;
    case XPRESSURE:
        if (sim->Q[k] < -Qtol)            status = CLOSED;
        break;
    }
    return(status);
}


static char  fcvstatus(ENsimulation_t sim, int k, char s, float h1, float h2)
/*
**-----------------------------------------------------------
**  Input:   k    = link index
**           s    = current status
**           h1   = head at upstream node
**           h2   = head at downstream node
**  Output:  returns new valve status
**  Purpose: updates status of a flow control valve.
**
**    Valve status changes to XFCV if flow reversal.
**    If current status is XFCV and current flow is
**    above setting, then valve becomes active.
**    If current status is XFCV, and current flow
**    positive but still below valve setting, then
**    status remains same.
**-----------------------------------------------------------
*/
{
    char  status;        /* New valve status */
    status = s;
    if (h1 - h2 < -Htol)                status = XFCV;
    else if (s == XFCV && sim->Q[k] >= sim->K[k]) 
        status = ACTIVE;
    return(status);
}


/*** Updated 9/7/00 ***/
/*** Updated 11/19/01 ***/
static void  tankstatus(ENsimulation_t sim, int k, int n1, int n2)
/*
**----------------------------------------------------------------
**  Input:   k  = link index
**           n1 = start node of link
**           n2 = end node of link
**  Output:  none
**  Purpose: closes link flowing into full or out of empty tank
**----------------------------------------------------------------
*/
{
    int   i,n;
    float h,q;

    /* Make node n1 be the tank */
    q = sim->Q[k];
    i = n1 - Njuncs;
    if (i <= 0)
    {
        i = n2 - Njuncs;
        if (i <= 0) return;
        n = n1;
        n1 = n2;
        n2 = n;
        q = -q;
    }
    h = sim->H[n1] - sim->H[n2];

    /* Skip reservoirs & closed links */
    if (Tank[i].A == 0.0 || sim->S[k] <= CLOSED) return;

    /* If tank full, then prevent flow into it */
    if (sim->H[n1] >= Tank[i].Hmax - Htol)
    {

        /* Case 1: Link is a pump discharging into tank */
        if ( Link[k].Type == PUMP )
        {
            if (Link[k].N2 == n1) sim->S[k] = TEMPCLOSED;
        }

        /* Case 2: Downstream head > tank head */
        /* (i.e., an open outflow check valve would close) */
        else if (cvstatus(OPEN, h, q) == CLOSED) sim->S[k] = TEMPCLOSED;
    }

    /* If tank empty, then prevent flow out of it */
    if (sim->H[n1] <= Tank[i].Hmin + Htol)
    {

        /* Case 1: Link is a pump discharging from tank */
        if ( Link[k].Type == PUMP)
        {
            if (Link[k].N1 == n1) sim->S[k] = TEMPCLOSED;
        }

        /* Case 2: Tank head > downstream head */
        /* (i.e., a closed outflow check valve would open) */
        else if (cvstatus(CLOSED, h, q) == OPEN) sim->S[k] = TEMPCLOSED;
    }
}                        /* End of tankstatus */


static int  pswitch(ENsimulation_t sim)
/*
**--------------------------------------------------------------
**  Input:   none
**  Output:  returns 1 if status of any link changes, 0 if not
**  Purpose: adjusts settings of links controlled by junction
**           pressures after a hydraulic solution is found
**--------------------------------------------------------------
*/
{
    int   i,                 /* Control statement index */
        k,                 /* Link being controlled */
        n,                 /* Node controlling link */
        reset,             /* Flag on control conditions */
        change,            /* Flag for status or setting change */
        anychange = 0;     /* Flag for 1 or more changes */
    char  s;                 /* Current link status */

    /* Check each control statement */
    for (i = 1; i <= sim->Ncontrols; i++)
    {
        reset = 0;
        if ( (k = sim->Control[i].Link) <= 0) continue;

        /* Determine if control based on a junction, not a tank */
        if ( (n = sim->Control[i].Node) > 0 && n <= Njuncs)
        {
            /* Determine if control conditions are satisfied */
            if (sim->Control[i].Type == LOWLEVEL
                && sim->H[n] <= sim->Control[i].Grade + Htol )
                reset = 1;
            if (sim->Control[i].Type == HILEVEL
                && sim->H[n] >= sim->Control[i].Grade - Htol )
                reset = 1;
        }

        /* Determine if control forces a status or setting change */
        if (reset == 1)
        {
            change = 0;
            s = sim->S[k];
            if (Link[k].Type == PIPE)
            {
                if (s != sim->Control[i].Status) change = 1;
            }
            if (Link[k].Type == PUMP)
            {
                if (sim->K[k] != sim->Control[i].Setting) change = 1;
            }
            if (Link[k].Type >= PRV)
            {
                if (sim->K[k] != sim->Control[i].Setting) change = 1;
                else if (sim->K[k] == MISSING 
                         && s != sim->Control[i].Status) change = 1;
            }

            /* If a change occurs, update status & setting */
            if (change)
            {
                sim->S[k] = sim->Control[i].Status;
                if (Link[k].Type > PIPE) sim->K[k] = sim->Control[i].Setting;
#if 0
                if (Statflag == FULL) writestatchange (k, s, sim->S[k]);
#endif
                /* Re-set flow if status has changed */
                if (sim->S[k] != s) {
                    if (Link[k].Type == PUMP) 
                        pumpswitch (sim, k, sim->S[k]);
                    sim->Q[k] = initlinkflow (k, sim->S[k], sim->K[k]);
                }
                anychange = 1;
            }
        }
    }
    return(anychange);
}                        /* End of pswitch */


static REAL newflows(ENsimulation_t sim)
/*
**----------------------------------------------------------------
**  Input:   none
**  Output:  returns solution convergence error
**  Purpose: updates link flows after new nodal heads computed
**----------------------------------------------------------------
*/
{
    REAL  dh,                    /* Link head loss       */
        dq;                    /* Link flow change     */
    REAL  dqsum,                 /* Network flow change  */
        qsum;                  /* Network total flow   */
    int   k, n, n1, n2;

    /* Initialize net inflows (i.e., demands) at tanks */
    for (n=Njuncs+1; n <= Nnodes; n++) 
        sim->D[n] = 0.0;

    /* Initialize sum of flows & corrections */
    qsum  = 0.0;
    dqsum = 0.0;

    /* Update flows in all links */
    for (k=1; k<=Nlinks; k++)
    {

        /*
        ** Apply flow update formula:
        **   dq = Y - P*(new head loss)
        **    P = 1/(dh/dq)
        **    Y = P*(head loss based on current flow)
        ** where P & Y were computed in newcoeffs().
        */

        n1 = Link[k].N1;
        n2 = Link[k].N2;
        dh = sim->H[n1] - sim->H[n2];
        dq = sim->Y[k] - sim->P[k] * dh;

        /* Prevent flow in constant HP pumps from going negative */
        if (Link[k].Type == PUMP)
        {
            n = PUMPINDEX(k);
            if (Pump[n].Ptype == CONST_HP && dq > sim->Q[k]) 
                dq = sim->Q[k]/2.0;
        }
        sim->Q[k] -= dq;

        /* Update sum of absolute flows & flow corrections */
        qsum += ABS (sim->Q[k]);
        dqsum += ABS(dq);

        /* Update net flows to tanks */
        if (n1 > Njuncs) sim->D[n1] -= sim->Q[k];
        if (n2 > Njuncs) sim->D[n2] += sim->Q[k];

    }

    /* Update emitter flows */
    for (k=1; k<=Njuncs; k++)
    {
        if (Node[k].Ke == 0.0) continue;
        dq = emitflowchange(sim, k);
        sim->E[k] -= dq;
        qsum += ABS (sim->E[k]);
        dqsum += ABS(dq);
    }

    /* Return ratio of total flow corrections to total flow */
    if (qsum > Hacc) return(dqsum/qsum);
    else return(dqsum);

}                        /* End of newflows */


static void   newcoeffs(ENsimulation_t sim)
/*
**--------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Purpose: computes coefficients of linearized network eqns.
**--------------------------------------------------------------
*/
{
    /* Reset coeffs. to 0 */
    memset (sim->Aii, 0, (Nnodes + 1) * sizeof(REAL));
    memset (sim->Aij, 0, (smatrix_non_zero_coeffs (sim->smatrix) + 1) 
            * sizeof(REAL));
    memset (sim->F, 0, (Nnodes+1) * sizeof(REAL));
    memset (sim->X, 0, (Nnodes+1) * sizeof(float));
    memset (sim->P, 0, (Nlinks+1) * sizeof(float));
    memset (sim->Y, 0, (Nlinks+1) * sizeof(float));
    linkcoeffs (sim);                            /* Compute link coeffs.  */
    emittercoeffs (sim);                         /* Compute emitter coeffs.*/
    nodecoeffs (sim);                            /* Compute node coeffs.  */
    valvecoeffs (sim);                           /* Compute valve coeffs. */
}                        /* End of newcoeffs */


static void  linkcoeffs(ENsimulation_t sim)
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  none
**   Purpose: computes solution matrix coefficients for links
**--------------------------------------------------------------
*/
{
    int   k,n1,n2;

    /* Examine each link of network */
    for (k=1; k<=Nlinks; k++)
    {
        n1 = Link[k].N1;           /* Start node of link */
        n2 = Link[k].N2;           /* End node of link   */

        /* Compute P[k] = 1 / (dh/dQ) and Y[k] = h * P[k]   */
        /* for each link k (where h = link head loss).      */
        /* FCVs, PRVs, and PSVs with non-fixed status       */
        /* are analyzed later.                              */

        switch (Link[k].Type)
        {
        case CV:
        case PIPE:  pipecoeff (sim, k, Link[k].Km); break;
        case PUMP:  pumpcoeff (sim, k); break;
        case PBV:   pbvcoeff (sim, k);  break;
        case TCV:   tcvcoeff (sim, k);  break;
        case GPV:   gpvcoeff (sim, k);  break;
        case FCV:
        case PRV:
        case PSV:   /* If valve status fixed then treat as pipe */
            /* otherwise ignore the valve for now. */
            if (sim->K[k] == MISSING) pipecoeff(sim, k, Link[k].Km);
            else continue;
            break;
        default:    continue;
        }

        /* Update net nodal inflows (X), solution matrix (A) and RHS
           array (F). Use covention that flow out of node is (-), flow
           into node is (+) */
        sim->X[n1] -= sim->Q[k];
        sim->X[n2] += sim->Q[k];
        
        /* Off-diagonal coeff. */
        sim->Aij[sim->smatrix->Ndx[k]] -= sim->P[k];  

        if (n1 <= Njuncs)                 /* Node n1 is junction */
        {
            sim->Aii[sim->smatrix->Row[n1]] += sim->P[k]; /* Diagonal coeff. */
            sim->F[sim->smatrix->Row[n1]] += sim->Y[k];   /* RHS coeff.      */
        }
        else   /* Node n1 is a tank   */
            sim->F[sim->smatrix->Row[n2]] += (sim->P[k] * sim->H[n1]);

        if (n2 <= Njuncs)                 /* Node n2 is junction */
        {
            sim->Aii[sim->smatrix->Row[n2]] += sim->P[k]; /* Diagonal coeff. */
            sim->F[sim->smatrix->Row[n2]] -= sim->Y[k];   /* RHS coeff.      */
        }
        else   /* Node n2 is a tank   */
            sim->F[sim->smatrix->Row[n1]] += (sim->P[k] * sim->H[n2]);
    }
}                        /* End of linkcoeffs */


static void  nodecoeffs(ENsimulation_t sim)
/*
**----------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Purpose: completes calculation of nodal flow imbalance (X)
**           & flow correction (F) arrays
**----------------------------------------------------------------
*/
{
    int   i;

    /* For junction nodes, subtract demand flow from net */
    /* flow imbalance & add imbalance to RHS array F.    */
    for (i=1; i<=Njuncs; i++)
    {
        sim->X[i] -= sim->D[i];
        sim->F[sim->smatrix->Row[i]] += sim->X[i];
    }
}                        /* End of nodecoeffs */


static void  valvecoeffs(ENsimulation_t sim)
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  none
**   Purpose: computes matrix coeffs. for PRVs, PSVs & FCVs
**            whose status is not fixed to OPEN/CLOSED
**--------------------------------------------------------------
*/
{
    int i,k,n1,n2;

    for (i=1; i<=Nvalves; i++)                   /* Examine each valve   */
    {
        k = Valve[i].Link;                        /* Link index of valve  */
        if (sim->K[k] == MISSING) continue;/* Valve status fixed   */
        n1 = Link[k].N1;                          /* Start & end nodes    */
        n2 = Link[k].N2;
        switch (Link[k].Type)                     /* Call valve-specific  */
        {                                         /*   function           */
        case PRV:  prvcoeff (sim, k, n1, n2); break;
        case PSV:  psvcoeff (sim, k, n1, n2); break;
        case FCV:  fcvcoeff (sim, k, n1, n2); break;
        default:   continue;
        }
    }
}                        /* End of valvecoeffs */


static void  emittercoeffs(ENsimulation_t sim)
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  none
**   Purpose: computes matrix coeffs. for emitters
**
**   Note: Emitters consist of a fictitious pipe connected to
**         a fictitious reservoir whose elevation equals that
**         of the junction. The headloss through this pipe is
**         Ke*(Flow)^Qexp, where Ke = emitter headloss coeff.
**--------------------------------------------------------------
*/
{
    int   i;
    REAL  ke;
    REAL  p;
    REAL  q;
    REAL  y;
    REAL  z;
    for (i=1; i<=Njuncs; i++)
    {
        if (Node[i].Ke == 0.0) continue;
        ke = MAX(CSMALL, Node[i].Ke);
        q = sim->E[i];
        z = ke*pow(ABS(q),Qexp);
        p = Qexp*z/ABS(q);
        if (p < RQtol) p = 1.0/RQtol;
        else p = 1.0/p;
        y = SGN(q)*z*p;
        sim->Aii[sim->smatrix->Row[i]] += p;
        sim->F[sim->smatrix->Row[i]] += y + p*Node[i].El;
        sim->X[i] -= q;
    }
}


static REAL  emitflowchange(ENsimulation_t sim, int i)
/*
**--------------------------------------------------------------
**   Input:   i = node index
**   Output:  returns change in flow at an emitter node
**   Purpose: computes flow change at an emitter node
**--------------------------------------------------------------
*/
{
    REAL ke, p;
    ke = MAX(CSMALL, Node[i].Ke);
    p = Qexp * ke * pow (ABS (sim->E[i]), (Qexp - 1.0));
    if (p < RQtol)
        p = 1/RQtol;
    else
        p = 1.0/p;
    return sim->E[i] / Qexp - p * (sim->H[i] - Node[i].El);
}


static inline
void  pipecoeff(ENsimulation_t sim, int k, REAL ml)
/*
**--------------------------------------------------------------
**   Input:   k = link index
              ml = Minor loss coeff.
**   Output:  none
**  Purpose:  computes P & Y coefficients for pipe k
**
**    P = inverse head loss gradient = 1/(dh/dQ)
**    Y = flow correction term = h*P
**--------------------------------------------------------------
*/
{
    REAL  hpipe,     /* Normal head loss          */
        hml,       /* Minor head loss           */
        p,         /* q*(dh/dq)                 */
        q,         /* Abs. value of flow        */
        r,         /* Resistance coeff.         */
        r1,        /* Total resistance factor   */
        f,         /* D-W friction factor       */
        dfdq;      /* Derivative of fric. fact. */

    /* For closed pipe use headloss formula: h = CBIG*q */
    if (sim->S[k] <= CLOSED)
    {
        sim->P[k] = INVCBIG;
        sim->Y[k] = sim->Q[k];
        return;
    }

    /* Evaluate headloss coefficients */
    q = fabs (sim->Q[k]);           /* Absolute flow       */
    r = Link[k].R;                         /* Resistance coeff.   */
    f = 1.0;                               /* D-W friction factor */
    if (Formflag == DW) 
        f = DWcoeff (Link[k].Type, Link[k].Diam, sim->Q[k], Link[k].Kc, &dfdq);

    r1 = f*r+ml;

    /* Use large P coefficient for small flow resistance product */
    if (r1*q < RQtol)
    {
        sim->P[k] = 1.0/RQtol;
        sim->Y[k] = sim->Q[k] / Hexp;
        return;
    }

    /* Compute P and Y coefficients */
    if (Formflag == DW)                  /* D-W eqn. */
    {
        hpipe = r1*SQR(q);                /* Total head loss */
        p = 2.0*r1*q;                     /* |dh/dQ| */
        /* + dfdq*r*q*q;*/                 /* Ignore df/dQ term */
        p = 1.0/p;
        sim->P[k] = p;
        sim->Y[k] = SGN (sim->Q[k]) * hpipe * p;
    }
    else                                 /* H-W or C-M eqn.   */
    {
        hpipe = r*pow(q,Hexp);            /* Friction head loss  */
        p = Hexp*hpipe;                   /* Q*dh(friction)/dQ   */
        if (ml > 0.0)
        {
            hml = ml*q*q;                  /* Minor head loss   */
            p += 2.0*hml;                  /* Q*dh(Total)/dQ    */
        }
        else  hml = 0.0;
        p = sim->Q[k] / p;          /* 1 / (dh/dQ) */
        sim->P[k] = fabs (p);
        sim->Y[k] = p * (hpipe + hml);
    }
}                        /* End of pipecoeff */


static REAL 
DWcoeff(char type, float diam, float flow, float roughness, REAL *dfdq)
/*
**--------------------------------------------------------------
**   Input:   k = link index
**   Output:  returns Darcy-Weisbach friction factor
**   Purpose: computes Darcy-Weisbach friction factor
**
**    Uses interpolating polynomials developed by
**    E. Dunlop for transition flow from 2000 < Re < 4000.
**
**   df/dq term is ignored as it slows convergence rate.
**--------------------------------------------------------------
*/
{
    REAL q,             /* Abs. value of flow */
        f;             /* Friction factor    */
    REAL x1,x2,x3,x4,
        y1,y2,y3,
        fa,fb,r;
    REAL s,w;

    *dfdq = 0.0;
    if (type > PIPE) return(1.0); /* Only apply to pipes */
    q = ABS (flow);
    s = Viscos * diam;
    w = q/s;                       /* w = Re(Pi/4) */
    if (w >= A1)                   /* Re >= 4000; Colebrook Formula */
    {
        y1 = A8/pow(w,0.9);
        y2 = roughness / (3.7 * diam) + y1;
        y3 = A9*log(y2);
        f = 1.0/SQR(y3);
        /*  *dfdq = (2.0+AA*y1/(y2*y3))*f; */   /* df/dq */
    }
    else if (w > A2)              /* Re > 2000; Interpolation formula */
    {
        y2 = roughness / (3.7 * diam) + AB;
        y3 = A9*log(y2);
        fa = 1.0/SQR(y3);
        fb = (2.0+AC/(y2*y3))*fa;
        r = w/A2;
        x1 = 7.0*fa - fb;
        x2 = 0.128 - 17.0*fa + 2.5*fb;
        x3 = -0.128 + 13.0*fa - (fb+fb);
        x4 = r*(0.032 - 3.0*fa + 0.5*fb);
        f = x1 + r*(x2 + r*(x3+x4));
        /*  *dfdq = (x1 + x1 + r*(3.0*x2 + r*(4.0*x3 + 5.0*x4)));  */
    }
    else if (w > A4)              /* Laminar flow: Hagen-Poiseuille Formula */
    {
        f = A3*s/q;
        /*  *dfdq = A3*s; */
    }
    else
    {
        f = 8.0;
        *dfdq = 0.0;
    }
    return(f);
}                        /* End of DWcoeff */


/*** Updated 10/25/00 ***/
/*** Updated 12/29/00 ***/
static void  pumpcoeff(ENsimulation_t sim, int k)
/*
**--------------------------------------------------------------
**   Input:   k = link index
**   Output:  none
**   Purpose: computes P & Y coeffs. for pump in link k
**--------------------------------------------------------------
*/
{
    int   p;         /* Pump index             */
    float h0,        /* Shutoff head           */
        q,         /* Abs. value of flow     */
        r,         /* Flow resistance coeff. */
        n;         /* Flow exponent coeff.   */

    /* Use high resistance pipe if pump closed or cannot deliver head */
    if (sim->S[k] <= CLOSED || sim->K[k] == 0.0)
    {
        pipecoeff(sim, k, Link[k].Km);
        return;
    }

    q = fabs (sim->Q[k]);
    q = MAX(q,TINY);
    p = PUMPINDEX(k);

    /* Get pump curve coefficients for custom pump curve. */
    if (Pump[p].Ptype == CUSTOM)
    {
        /* Find intercept (h0) & slope (r) of pump curve    */
        /* line segment which contains speed-adjusted flow. */
        curvecoeff (Pump[p].Hcurve, q / sim->K[k], &h0, &r);

        /* Determine head loss coefficients. */
        Pump[p].H0 = -h0;
        Pump[p].R  = -r;
        Pump[p].N  = 1.0;
    }

    /* Adjust head loss coefficients for pump speed. */
    h0 = SQR (sim->K[k]) * Pump[p].H0;
    n  = Pump[p].N;
    r  = Pump[p].R * pow (sim->K[k], 2.0 - n);
    if (n != 1.0) r = n*r*pow(q,n-1.0);

    /* Compute inverse headloss gradient (P) and flow correction factor (Y) */
    sim->P[k] = 1.0/MAX(r,RQtol);
    sim->Y[k] = sim->Q[k]/n + sim->P[k] * h0;
}                        /* End of pumpcoeff */


/*** Updated 10/25/00 ***/
/*** Updated 12/29/00 ***/
static void  curvecoeff(int i, float q, float *h0, float *r)
/*
**-------------------------------------------------------------------
**   Input:   i   = curve index
**            q   = flow rate
**   Output:  *h0 = head at zero flow (y-intercept)
**            *r  = dHead/dFlow (slope)
**   Purpose: computes intercept and slope of head v. flow curve
**            at current flow.
**-------------------------------------------------------------------
*/
{
    int   k1, k2, npts;
    float *x, *y;

    /* Remember that curve is stored in untransformed units */
    q *= Ucf[FLOW];
    x = Curve[i].X;           /* x = flow */
    y = Curve[i].Y;           /* y = head */
    npts = Curve[i].Npts;

    /* Find linear segment of curve that brackets flow q */
    k2 = 0;
    while (k2 < npts && x[k2] < q) k2++;
    if      (k2 == 0)    k2++;
    else if (k2 == npts) k2--;
    k1  = k2 - 1;

    /* Compute slope and intercept of this segment */
    *r  = (y[k2]-y[k1])/(x[k2]-x[k1]);
    *h0 = y[k1] - (*r)*x[k1];

    /* Convert units */
    *h0 = (*h0)/Ucf[HEAD];
    *r  = (*r)*Ucf[FLOW]/Ucf[HEAD];
}                       /* End of curvecoeff */


static void  gpvcoeff(ENsimulation_t sim, int k)
/*
**--------------------------------------------------------------
**   Input:   k = link index
**   Output:  none
**   Purpose: computes P & Y coeffs. for general purpose valve
**--------------------------------------------------------------
*/
{
    float h0,        /* Headloss curve intercept */
        q,         /* Abs. value of flow       */
        r;         /* Flow resistance coeff.   */

/*** Updated 9/7/00 ***/
    /* Treat as a pipe if valve closed */
    if (sim->S[k] == CLOSED) pipecoeff (sim, k, Link[k].Km);

    /* Otherwise utilize headloss curve   */
    /* whose index is stored in K */
    else
    {
        /* Find slope & intercept of headloss curve. */
        q = ABS (sim->Q[k]);
        q = MAX (q, TINY);

/*** Updated 10/25/00 ***/
/*** Updated 12/29/00 ***/
        curvecoeff ((int)ROUND (sim->K[k]), q, &h0, &r);

        /* Compute inverse headloss gradient (P) and flow  */
        /* correction factor (Y) using formulas for pumps. */
        sim->P[k] = 1.0 / MAX (r, RQtol);
        sim->Y[k] = sim->Q[k] + sim->P[k] * h0;
    }
}


static void  pbvcoeff(ENsimulation_t sim, int k)
/*
**--------------------------------------------------------------
**   Input:   k = link index
**   Output:  none
**   Purpose: computes P & Y coeffs. for pressure breaker valve
**--------------------------------------------------------------
*/
{
    /* If valve fixed OPEN or CLOSED then treat as a pipe */
    if (sim->K[k] == MISSING || sim->K[k] == 0.0) 
        pipecoeff (sim, k, Link[k].Km);

    /* If valve is active */
    else
    {
        /* Treat as a pipe if minor loss > valve setting */
        if (Link[k].Km * SQR (sim->Q[k]) > sim->K[k]) 
            pipecoeff (sim, k, Link[k].Km);

        /* Otherwise force headloss across valve to be equal to setting */
        else
        {
            sim->P[k] = CBIG;
            sim->Y[k] = sim->K[k] * CBIG;
        }
    }
}                        /* End of pbvcoeff */


static void  tcvcoeff(ENsimulation_t sim, int k)
/*
**--------------------------------------------------------------
**   Input:   k = link index
**   Output:  none
**   Purpose: computes P & Y coeffs. for throttle control valve
**--------------------------------------------------------------
*/
{
    float km;

    /* If valve not fixed OPEN or CLOSED, compute its loss coeff. */
    if (sim->K[k] != MISSING)
        km = 0.02517 * sim->K[k] 
            / (SQR (Link[k].Diam) * SQR (Link[k].Diam));

    /* Then apply usual pipe formulas */
    pipecoeff (sim, k, km);
    
}                        /* End of tcvcoeff */


static void  prvcoeff(ENsimulation_t sim, int k, int n1, int n2)
/*
**--------------------------------------------------------------
**   Input:   k    = link index
**            n1   = upstream node of valve
**            n2   = downstream node of valve
**   Output:  none
**   Purpose: computes solution matrix coeffs. for pressure
**            reducing valves
**--------------------------------------------------------------
*/
{
    int   i,j;                       /* Rows of solution matrix */
    float hset;                      /* Valve head setting      */
    i  = sim->smatrix->Row[n1]; /* Matrix rows of nodes    */
    j  = sim->smatrix->Row[n2];
    hset   = Node[n2].El + sim->K[k];     /* Valve setting           */

    if (sim->S[k] == ACTIVE)
    {
        /*
          Set coeffs. to force head at downstream
          node equal to valve setting & force flow (when updated in
          newflows()) equal to flow imbalance at downstream node.
        */
        sim->P[k] = 0.0;
        /* Force flow balance   */
        sim->Y[k] = sim->Q[k] + sim->X[n2]; 
        sim->F[j] += (hset*CBIG);       /* Force head = hset    */
        sim->Aii[j] += CBIG;            /*   at downstream node */
        if (sim->X[n2] < 0.0) 
            sim->F[i] +=  sim->X[n2];
        return;
    }
    /*
      For OPEN, CLOSED, or XPRESSURE valve
      compute matrix coeffs. using the pipecoeff() function.
    */
    pipecoeff (sim, k, Link[k].Km);
    sim->Aij[sim->smatrix->Ndx[k]] -= sim->P[k];
    sim->Aii[i] += sim->P[k];
    sim->Aii[j] += sim->P[k];
    sim->F[i] += (sim->Y[k] - sim->Q[k]);
    sim->F[j] -= (sim->Y[k] - sim->Q[k]);
}                        /* End of prvcoeff */


static void  psvcoeff(ENsimulation_t sim, int k, int n1, int n2)
/*
**--------------------------------------------------------------
**   Input:   k    = link index
**            n1   = upstream node of valve
**            n2   = downstream node of valve
**   Output:  none
**   Purpose: computes solution matrix coeffs. for pressure
**            sustaining valve
**--------------------------------------------------------------
*/
{
    int   i,j;                       /* Rows of solution matrix */
    float hset;                      /* Valve head setting      */
    i  = sim->smatrix->Row[n1];      /* Matrix rows of nodes    */
    j  = sim->smatrix->Row[n2];
    hset   = Node[n1].El + sim->K[k];     /* Valve setting           */

    if (sim->S[k] == ACTIVE)
    {
        /*
          Set coeffs. to force head at upstream
          node equal to valve setting & force flow (when updated in
          newflows()) equal to flow imbalance at upstream node.
        */
        sim->P[k] = 0.0;
        /* Force flow balance   */
        sim->Y[k] = sim->Q[k] - sim->X[n1];
        sim->F[i] += (hset*CBIG);  /* Force head = hset    */
        sim->Aii[i] += CBIG;       /*   at upstream node   */
        if (sim->X[n1] > 0.0) 
            sim->F[j] += sim->X[n1];
        return;
    }

    /*
      For OPEN, CLOSED, or XPRESSURE valve
      compute matrix coeffs. using the pipecoeff() function.
    */
    pipecoeff (sim, k, Link[k].Km);
    sim->Aij[sim->smatrix->Ndx[k]] -= sim->P[k];
    sim->Aii[i] += sim->P[k];
    sim->Aii[j] += sim->P[k];
    sim->F[i] += (sim->Y[k] - sim->Q[k]);
    sim->F[j] -= (sim->Y[k] - sim->Q[k]);
}                        /* End of psvcoeff */


static void  fcvcoeff(ENsimulation_t sim, int k, int n1, int n2)
/*
**--------------------------------------------------------------
**   Input:   k    = link index
**            n1   = upstream node of valve
**            n2   = downstream node of valve
**   Output:  none
**   Purpose: computes solution matrix coeffs. for flow control
**            valve
**--------------------------------------------------------------
*/
{
    int   i,j;                   /* Rows in solution matrix */
    float q;                     /* Valve flow setting      */
    q = sim->K[k];
    i = sim->smatrix->Row[n1];
    j = sim->smatrix->Row[n2];

    /*
      If valve active, break network at valve and treat
      flow setting as external demand at upstream node
      and external supply at downstream node.
    */
    if (sim->S[k] == ACTIVE)
    {
        sim->X[n1] -= q;
        sim->F[i] -= q;
        sim->X[n2] += q;
        sim->F[j] += q;
        sim->P[k] = 0.0;
        sim->Y[k] = sim->Q[k] - q;
    }
    /*
      Otherwise treat valve as an open pipe
    */
    else
    {
        pipecoeff (sim, k, Link[k].Km);
        sim->Aij[sim->smatrix->Ndx[k]] -= sim->P[k];
        sim->Aii[i] += sim->P[k];
        sim->Aii[j] += sim->P[k];
        sim->F[i] += (sim->Y[k] - sim->Q[k]);
        sim->F[j] -= (sim->Y[k] - sim->Q[k]);
    }
}                        /* End of fcvcoeff */

/* 
   This function updates the number of switches and the schedule
   of a pump. Invoke this function whenever a pump changes status.
   
   index : link index of the pump.
   status : new status of the pump.
*/
static void
pumpswitch(ENsimulation_t sim, int index, char status)
{
    int p = PUMPINDEX(index);

    /* Htime == 0 is ignored because if the final status is the same
       as the initial status, there is no pump switch.  */
    if (sim->Htime > 0 && status != Link[index].Stat) { 
        sim->Pump[p].Nswitches++; 
    }

    if (sim->Htime < 24 * 3600 && sim->Pump[p].Schedule_idx < 24) {
        sim->Pump[p].Schedule[sim->Pump[p].Schedule_idx] = sim->Htime;
        sim->Pump[p].Schedule_idx++;
    }
}

/****************  END OF HYDRAUL.C  ***************/
