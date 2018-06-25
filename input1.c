/*
*********************************************************************

INPUT1.C -- Input Functions for EPANET Program

AUTHOR:     M. Lopez-Ibanez (M. L\'opez-Ib\'a\~nez)
            Napier University, Edinburgh, UK.

$LastChangedBy: manu $ $Revision: 158 $
$LastChangedDate: 2007-12-19 13:01:00 +0100 (Wed, 19 Dec 2007) $

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
DATE:       5/30/00
            9/7/00
            11/19/01
            6/24/02
AUTHOR:     L. Rossman
            US EPA - NRMRL

  This module initializes, retrieves, and adjusts the input
  data for a network simulation.

  The entry point for this module is:
     getdata() -- called from ENopen() in EPANET.C.

*********************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
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

/*
  --------------------- Module Global Variables  ----------------------
*/

#define   MAXITER   40       /* Default max. # hydraulic iterations    */
#define   HACC      0.001    /* Default hydraulics convergence ratio   */
#define   HTOL      0.0005   /* Default hydraulic head tolerance (ft)  */

/*** Updated 11/19/01 ***/
#define   QTOL      0.0001   /* Default flow rate tolerance (cfs)      */

#define   AGETOL    0.01     /* Default water age tolerance (hrs)      */
#define   CHEMTOL   0.01     /* Default concentration tolerance        */
#define   PAGESIZE  0        /* Default uses no page breaks            */
#define   SPGRAV    1.0      /* Default specific gravity               */
#define   EPUMP     75       /* Default pump efficiency                */
#define   DEFPATID  "1"      /* Default demand pattern ID              */

/*
  These next three parameters are used in the hydraulics solver:
*/

#define   RQTOL     1E-7     /* Default low flow resistance tolerance  */
#define   CHECKFREQ 2        /* Default status check frequency         */
#define   MAXCHECK  10       /* Default # iterations for status checks */

extern char *Fldname[];      /* Defined in enumstxt.h in EPANET.C      */
extern char *RptFlowUnitsTxt[];


static void  convertunits(void);             /* Converts data to std. units  */


int  getdata()
/*
**----------------------------------------------------------------
**  Input:   none
**  Output:  returns error code
**  Purpose: reads in network data from disk file
**----------------------------------------------------------------
*/
{
    int errcode = 0;
    setdefaults();                /* Assign default data values     */
#if 0
    initreport();                 /* Initialize reporting options   */
#endif
    rewind(InFile);               /* Rewind input file              */
    ERRCODE(readdata());          /* Read in network data           */
    if (!errcode) adjustdata();   /* Adjust data for default values */
    if (!errcode) initunits();    /* Initialize units on input data */
    ERRCODE(inittanks());         /* Initialize tank volumes        */
    if (!errcode) convertunits(); /* Convert units on input data    */
    return(errcode);
}                       /*  End of getdata  */


void  setdefaults()
/*
**----------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Purpose: assigns default values to global variables
**----------------------------------------------------------------
*/
{
    strncpy(Title[0],"",MAXMSG);
    strncpy(Title[1],"",MAXMSG);
    strncpy(Title[2],"",MAXMSG);
    strncpy(HydFname,"",MAXFNAME);
    strncpy(MapFname,"",MAXFNAME);
    strncpy(ChemName,t_CHEMICAL,MAXID);
    strncpy(ChemUnits,u_MGperL,MAXID);
    strncpy(DefPatID,DEFPATID,MAXID);
#if 0
    Hydflag   = SCRATCH;         /* No external hydraulics file    */
    Qualflag  = NONE;            /* No quality simulation          */
#endif
    Formflag  = HW;              /* Use Hazen-Williams formula     */
    Unitsflag = US;              /* US unit system                 */
    Flowflag  = GPM;             /* Flow units are gpm             */
    Pressflag = PSI;             /* Pressure units are psi         */
    Tstatflag = SERIES;          /* Generate time series output    */
    Htol      = HTOL;            /* Default head tolerance         */
    Qtol      = QTOL;            /* Default flow tolerance         */
    Hacc      = HACC;            /* Default hydraulic accuracy     */
#if 0
    Ctol      = MISSING;         /* No pre-set quality tolerance   */
#endif
    MaxIter   = MAXITER;         /* Default max. hydraulic trials  */
    ExtraIter = -1;              /* Stop if network unbalanced     */
    Dur       = 0;               /* 0 sec duration (steady state)  */
    Tstart    = 0;               /* Starting time of day           */
    Pstart    = 0;               /* Starting pattern period        */
    Hstep     = 3600;            /* 1 hr hydraulic time step       */
#if 0
    Qstep     = 0;               /* No pre-set quality time step   */
#endif
    Pstep     = 3600;            /* 1 hr time pattern period       */
    Rstep     = 3600;            /* 1 hr reporting period          */
    Rulestep  = 0;               /* No pre-set rule time step      */
    Rstart    = 0;               /* Start reporting at time 0      */
#if 0
    TraceNode = 0;               /* No source tracing              */
    BulkOrder = 1.0;             /* 1st-order bulk reaction rate   */
    WallOrder = 1.0;             /* 1st-order wall reaction rate   */
    TankOrder = 1.0;             /* 1st-order tank reaction rate   */
    Kbulk     = 0.0;             /* No global bulk reaction        */
    Kwall     = 0.0;             /* No global wall reaction        */
    Climit    = 0.0;             /* No limiting potential quality  */
    Diffus    = MISSING;         /* Temporary diffusivity          */
    Rfactor   = 0.0;             /* No roughness-reaction factor   */
#endif
    Viscos    = MISSING;         /* Temporary viscosity            */
    SpGrav    = SPGRAV;          /* Default specific gravity       */
    DefPat    = 0;               /* Default demand pattern index   */
    Epat      = 0;               /* No energy price pattern        */
    Ecost     = 0.0;             /* Zero unit energy cost          */
    Dcost     = 0.0;             /* Zero energy demand charge      */
    Epump     = EPUMP;           /* Default pump efficiency        */
    Qexp      = 2.0;             /* Flow exponent for emitters     */
    Dmult     = 1.0;             /* Demand multiplier              */
    RQtol     = RQTOL;           /* Default hydraulics parameters  */
    CheckFreq = CHECKFREQ;
    MaxCheck  = MAXCHECK;
}                       /*  End of setdefaults  */

#if 0
void  initreport()
/*
**----------------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Purpose: initializes reporting options
**----------------------------------------------------------------------
*/
{
    int i;
    strncpy(Rpt2Fname,"",MAXFNAME);
    PageSize    = PAGESIZE;      /* Default page size for report   */
    Summaryflag = TRUE;          /* Write summary report           */
    Messageflag = TRUE;          /* Report error/warning messages  */
    Statflag    = FALSE;         /* No hydraulic status reports    */
    Energyflag  = FALSE;         /* No energy usage report         */
    Nodeflag    = 0;             /* No reporting on nodes          */
    Linkflag    = 0;             /* No reporting on links          */
    for (i=0; i<MAXVAR; i++)     /* For each reporting variable:   */
    {
        strncpy(Field[i].Name,Fldname[i],MAXID);
        Field[i].Enabled = FALSE;        /* Not included in report  */
        Field[i].Precision = 2;          /* 2 decimal precision     */

/*** Updated 6/24/02 ***/
        Field[i].RptLim[LOW] =   SQR(BIG); /* No reporting limits   */
        Field[i].RptLim[HI]  =  -SQR(BIG);
    }
    Field[FRICTION].Precision = 3;
    for (i=DEMAND; i<=QUALITY; i++) Field[i].Enabled = TRUE;
    for (i=FLOW; i<=HEADLOSS; i++) Field[i].Enabled = TRUE;
}
#endif

void  adjustdata()
/*
**----------------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Purpose: adjusts data after input file has been processed
**----------------------------------------------------------------------
*/
{
    int   i;
    float ucf;                   /* Unit conversion factor */
    Pdemand demand;              /* Pointer to demand record */

/* Use 1 hr pattern & report time step if none specified */
    if (Pstep <= 0) Pstep = 3600;
    if (Rstep == 0) Rstep = Pstep;

/* Hydraulic time step cannot be greater than pattern or report time step */
    if (Hstep <=  0)   Hstep = 3600;
    if (Hstep > Pstep) Hstep = Pstep;
    if (Hstep > Rstep) Hstep = Rstep;

/* Report start time cannot be greater than simulation duration */
    if (Rstart > Dur) Rstart = 0;

#if 0
/* No water quality analysis for steady state run */
    if (Dur == 0) Qualflag = NONE;
/* If no quality timestep, then make it 1/10 of hydraulic timestep */
    if (Qstep == 0) Qstep = Hstep/10;
#endif
/* If no rule time step then make it 1/10 of hydraulic time step; */
/* Rule time step cannot be greater than hydraulic time step */
    if (Rulestep == 0) Rulestep = Hstep/10;
    Rulestep = MIN(Rulestep, Hstep);
#if 0
/* Quality timestep cannot exceed hydraulic timestep */
    Qstep = MIN(Qstep, Hstep);

/* If no quality tolerance, then use default values */
    if (Ctol == MISSING)
    {
        if (Qualflag == AGE) Ctol = AGETOL;
        else                 Ctol = CHEMTOL;
    }
#endif
/* Determine unit system based on flow units */
    switch (Flowflag)
    {
    case LPS:          /* Liters/sec */
    case LPM:          /* Liters/min */
    case MLD:          /* megaliters/day  */
    case CMH:          /* cubic meters/hr */
    case CMD:          /* cubic meters/day */
        Unitsflag = SI;
        break;
    default:
        Unitsflag = US;
    }

/* Revise pressure units depending on flow units */
    if (Unitsflag != SI) Pressflag = PSI;
    else if (Pressflag == PSI) Pressflag = METERS;

/* Store value of viscosity & diffusivity */
    ucf = 1.0;
    if (Unitsflag == SI) ucf = SQR(MperFT);

    if (Viscos == MISSING)     /* No value supplied */
        Viscos = VISCOS;
    else if (Viscos > 1.e-3)   /* Multiplier supplied */
        Viscos = Viscos*VISCOS;
    else                       /* Actual value supplied */
        Viscos = Viscos/ucf;
#if 0
    if (Diffus == MISSING)
        Diffus = DIFFUS;
    else if (Diffus > 1.e-4)
        Diffus = Diffus*DIFFUS;
    else
        Diffus = Diffus/ucf;
#endif
/*
  Set exponent in head loss equation and adjust flow-resistance tolerance.
*/
    if (Formflag == HW) Hexp = 1.852;
    else                Hexp = 2.0;

/*** Updated 9/7/00 ***/
/*** No adjustment made to flow-resistance tolerance ***/
    /*RQtol = RQtol/Hexp;*/

/* See if default reaction coeffs. apply */
    for (i=1; i<=Nlinks; i++)
    {
        if (Link[i].Type > PIPE) continue;
#if 0
        if (Link[i].Kb == MISSING) Link[i].Kb = Kbulk;  /* Bulk coeff. */
        if (Link[i].Kw == MISSING)                      /* Wall coeff. */
        {
            /* Rfactor is the pipe roughness correlation factor */
            if (Rfactor == 0.0)   Link[i].Kw = Kwall;
            else if ((Link[i].Kc > 0.0) && (Link[i].Diam > 0.0))
            {
                if (Formflag == HW) Link[i].Kw = Rfactor/Link[i].Kc;
                if (Formflag == DW) Link[i].Kw = Rfactor/ABS(log(Link[i].Kc/Link[i].Diam));
                if (Formflag == CM) Link[i].Kw = Rfactor*Link[i].Kc;
            }
            else Link[i].Kw = 0.0;
        }
#endif
    }
#if 0
    for (i=1; i<=Ntanks; i++)
        if (Tank[i].Kb == MISSING) Tank[i].Kb = Kbulk;
#endif

/* Use default pattern if none assigned to a demand */
    for (i=1; i<=Nnodes; i++)
    {
        for (demand = Node[i].D; demand != NULL; demand = demand->next)
            if (demand->Pat == 0) demand->Pat = DefPat;
    }
#if 0
/* Remove QUALITY as a reporting variable if no WQ analysis */
    if (Qualflag == NONE) Field[QUALITY].Enabled = FALSE;
#endif
}                       /*  End of adjustdata  */


int  inittanks()
/*
**---------------------------------------------------------------
**  Input:   none
**  Output:  returns error code
**  Purpose: initializes volumes in non-cylindrical tanks
**---------------------------------------------------------------
*/
{
    int   i,j,n;
    float a;
    int   errcode = 0,
        levelerr;

    for (j=1; j<=Ntanks; j++)
    {

        /* Skip reservoirs */
        if (Tank[j].A == 0.0) continue;

        levelerr = 0;
        i = Tank[j].Vcurve;

        /* Check for valid lower/upper tank levels */
        if (Tank[j].H0   > Tank[j].Hmax ||
            Tank[j].Hmin > Tank[j].Hmax ||
            Tank[j].H0   < Tank[j].Hmin )  {

            levelerr = 1;
        }
        /* Else if tank has a volume curve, */
        else if(i > 0) {

            /* Check that tank heights are within volume curve */
            n = Curve[i].Npts - 1;
            if (Tank[j].Hmin < Curve[i].X[0] ||
                Tank[j].Hmax > Curve[i].X[n] )  {

                levelerr = 1;
            }
            else {
                /* Find min., max., and initial volumes from curve */
                Tank[j].Vmin = interp(Curve[i].Npts,Curve[i].X,
                                      Curve[i].Y,Tank[j].Hmin);
                Tank[j].Vmax = interp(Curve[i].Npts,Curve[i].X,
                                      Curve[i].Y,Tank[j].Hmax);
                Tank[j].V0   = interp(Curve[i].Npts,Curve[i].X,
                                      Curve[i].Y,Tank[j].H0);

                /* Find a "nominal" diameter for tank */
                a = (Curve[i].Y[n] - Curve[i].Y[0])/
                    (Curve[i].X[n] - Curve[i].X[0]);
                Tank[j].A = sqrt(4.0*a/PI);
            }
        }
        /* Report error in levels if found */
        if (levelerr)
        {
#if 0
            sprintf(Msg,ERR225,Node[Tank[j].Node].ID);
            writeline(Msg);
#endif
            errcode = 200;
        }
    }
    return(errcode);
}                       /* End of inittanks */


void  initunits()
/*
**--------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Purpose: determines unit conversion factors
**--------------------------------------------------------------
*/
{
    float  dcf,  /* distance conversion factor      */
#if 0
        ccf,  /* concentration conversion factor */
#endif
        qcf,  /* flow conversion factor          */
        hcf,  /* head conversion factor          */
        pcf,  /* pressure conversion factor      */
        wcf;  /* energy conversion factor        */

    if (Unitsflag == SI)                            /* SI units */
    {
        strcpy(Field[DEMAND].Units,RptFlowUnitsTxt[(int) Flowflag]);
        strcpy(Field[ELEV].Units,u_METERS);
        strcpy(Field[HEAD].Units,u_METERS);
        if (Pressflag == METERS) strcpy(Field[PRESSURE].Units,u_METERS);
        else strcpy(Field[PRESSURE].Units,u_KPA);
        strcpy(Field[LENGTH].Units,u_METERS);
        strcpy(Field[DIAM].Units,u_MMETERS);
        strcpy(Field[FLOW].Units,RptFlowUnitsTxt[(int) Flowflag]);
        strcpy(Field[VELOCITY].Units,u_MperSEC);
        strcpy(Field[HEADLOSS].Units,u_per1000M);
        strcpy(Field[FRICTION].Units,"");
        strcpy(Field[POWER].Units,u_KW);
        dcf = 1000.0*MperFT;
        qcf = LPSperCFS;
        if (Flowflag == LPM) qcf = LPMperCFS;
        if (Flowflag == MLD) qcf = MLDperCFS;
        if (Flowflag == CMH) qcf = CMHperCFS;
        if (Flowflag == CMD) qcf = CMDperCFS;
        hcf = MperFT;
        if (Pressflag == METERS) pcf = MperFT*SpGrav;
        else pcf = KPAperPSI*PSIperFT*SpGrav;
        wcf = KWperHP;
    }
    else                                         /* US units */
    {
        strcpy(Field[DEMAND].Units,RptFlowUnitsTxt[(int) Flowflag]);
        strcpy(Field[ELEV].Units,u_FEET);
        strcpy(Field[HEAD].Units,u_FEET);
        strcpy(Field[PRESSURE].Units,u_PSI);
        strcpy(Field[LENGTH].Units,u_FEET);
        strcpy(Field[DIAM].Units,u_INCHES);
        strcpy(Field[FLOW].Units,RptFlowUnitsTxt[(int) Flowflag]);
        strcpy(Field[VELOCITY].Units,u_FTperSEC);
        strcpy(Field[HEADLOSS].Units,u_per1000FT);
        strcpy(Field[FRICTION].Units,"");
        strcpy(Field[POWER].Units,u_HP);
        dcf = 12.0;
        qcf = 1.0;
        if (Flowflag == GPM) qcf = GPMperCFS;
        if (Flowflag == MGD) qcf = MGDperCFS;
        if (Flowflag == IMGD)qcf = IMGDperCFS;
        if (Flowflag == AFD) qcf = AFDperCFS;
        hcf = 1.0;
        pcf = PSIperFT*SpGrav;
        wcf = 1.0;
    }
#if 0
    strcpy(Field[QUALITY].Units,"");
    ccf = 1.0;
    if (Qualflag == CHEM)
    {
        ccf = 1.0/LperFT3;
        strncpy(Field[QUALITY].Units,ChemUnits,MAXID);
        strncpy(Field[REACTRATE].Units,ChemUnits,MAXID);
        strcat(Field[REACTRATE].Units,t_PERDAY);
    }
    else if (Qualflag == AGE) strcpy(Field[QUALITY].Units,u_HOURS);
    else if (Qualflag == TRACE) strcpy(Field[QUALITY].Units,u_PERCENT);
#endif
    Ucf[DEMAND]    = qcf;
    Ucf[ELEV]      = hcf;
    Ucf[HEAD]      = hcf;
    Ucf[PRESSURE]  = pcf;
#if 0
    Ucf[QUALITY]   = ccf;
#endif
    Ucf[LENGTH]    = hcf;
    Ucf[DIAM]      = dcf;
    Ucf[FLOW]      = qcf;
    Ucf[VELOCITY]  = hcf;
    Ucf[HEADLOSS]  = hcf;
#if 0
    Ucf[LINKQUAL]  = ccf;
    Ucf[REACTRATE] = ccf;
#endif
    Ucf[FRICTION]  = 1.0;
    Ucf[POWER]     = wcf;
    Ucf[VOLUME]    = hcf*hcf*hcf;
    if (Hstep < 1800)                    /* Report time in mins.    */
    {                                    /* if hydraulic time step  */
        Ucf[TIME] = 1.0/60.0;             /* is less than 1/2 hour.  */
        strcpy(Field[TIME].Units,u_MINUTES);
    }
    else
    {
        Ucf[TIME] = 1.0/3600.0;
        strcpy(Field[TIME].Units,u_HOURS);
    }

}                       /*  End of initunits  */


static void  convertunits()
/*
**--------------------------------------------------------------
**  Input:   none
**  Output:  none
**  Purpose: converts units of input data
**--------------------------------------------------------------
*/
{
    int   i,j,k;
    float ucf;        /* Unit conversion factor */
    Pdemand demand;   /* Pointer to demand record */

/* Convert nodal elevations & initial WQ */
/* (WQ source units are converted in QUALITY.C */
    for (i=1; i<=Nnodes; i++)
    {
        Node[i].El /= Ucf[ELEV];
#if 0
        Node[i].C0 /= Ucf[QUALITY];
#endif
    }

/* Convert demands */
    for (i=1; i<=Njuncs; i++)
    {
        for (demand = Node[i].D; demand != NULL; demand = demand->next)
            demand->Base /= Ucf[DEMAND];
    }

/* Convert emitter discharge coeffs. to head loss coeff. */
    ucf = pow(Ucf[FLOW],Qexp)/Ucf[PRESSURE];
    for (i=1; i<=Njuncs; i++)
        if (Node[i].Ke > 0.0) Node[i].Ke = ucf/pow(Node[i].Ke,Qexp);

/* Initialize tank variables (convert tank levels to elevations) */
    for (j=1; j<=Ntanks; j++)
    {
        i = Tank[j].Node;
        Tank[j].H0 = Node[i].El + Tank[j].H0/Ucf[ELEV];
        Tank[j].Hmin = Node[i].El + Tank[j].Hmin/Ucf[ELEV];
        Tank[j].Hmax = Node[i].El + Tank[j].Hmax/Ucf[ELEV];
        Tank[j].A = PI*SQR(Tank[j].A/Ucf[ELEV])/4.0;
        Tank[j].V0 /= Ucf[VOLUME];
        Tank[j].Vmin /= Ucf[VOLUME];
        Tank[j].Vmax /= Ucf[VOLUME];
#if 0
        Tank[j].Kb /= SECperDAY;
/*        Tank[j].V = Tank[j].V0;*/
        Tank[j].C = Node[i].C0;
        Tank[j].V1max *= Tank[j].Vmax;
#endif
    }
#if 0
/* Convert WQ option concentration units */
    Climit /= Ucf[QUALITY];
    Ctol   /= Ucf[QUALITY];

/* Convert global reaction coeffs. */
    Kbulk /= SECperDAY;
    Kwall /= SECperDAY;
#endif
    /* Convert units of link parameters */
    for (k=1; k<=Nlinks; k++)
    {
        if (Link[k].Type <= PIPE)
        {
            /* Convert pipe parameter units:                         */
            /*    - for Darcy-Weisbach formula, convert roughness    */
            /*      from millifeet (or mm) to ft (or m)              */
            /*    - for US units, convert diameter from inches to ft */
            if (Formflag  == DW) Link[k].Kc /= (1000.0*Ucf[ELEV]);
            Link[k].Diam /= Ucf[DIAM];
            Link[k].Len /= Ucf[LENGTH];

            /* Convert minor loss coeff. from V^2/2g basis to Q^2 basis */
            Link[k].Km = 0.02517*Link[k].Km/SQR(Link[k].Diam)/SQR(Link[k].Diam);
#if 0
            /* Convert units on reaction coeffs. */
            Link[k].Kb /= SECperDAY;
            Link[k].Kw /= SECperDAY;
#endif
        }

        else if (Link[k].Type == PUMP )
        {
            /* Convert units for pump curve parameters */
            i = PUMPINDEX(k);
            if (Pump[i].Ptype == CONST_HP)
            {
                /* For constant hp pump, convert kw to hp */
                if (Unitsflag == SI) Pump[i].R /= Ucf[POWER];
            }
            else
            {
                /* For power curve pumps, convert     */
                /* shutoff head and flow coefficient  */
                if (Pump[i].Ptype == POWER_FUNC)
                {
                    Pump[i].H0 /= Ucf[HEAD];
                    Pump[i].R  *= (pow(Ucf[FLOW],Pump[i].N)/Ucf[HEAD]);
                }
                /* Convert flow range & max. head units */
                Pump[i].Q0   /= Ucf[FLOW];
                Pump[i].Qmax /= Ucf[FLOW];
                Pump[i].Hmax /= Ucf[HEAD];
            }
        }

        else
        {
            /* For flow control valves, convert flow setting    */
            /* while for other valves convert pressure setting  */
            Link[k].Diam /= Ucf[DIAM];
            Link[k].Km = 0.02517*Link[k].Km/SQR(Link[k].Diam)/SQR(Link[k].Diam);
            if (Link[k].Kc != MISSING) switch (Link[k].Type)
            {
            case FCV: Link[k].Kc /= Ucf[FLOW]; break;
            case PRV:
            case PSV:
            case PBV: Link[k].Kc /= Ucf[PRESSURE]; break;
            }
        }

        /* Compute flow resistances */
        resistance(k);
    }

    /* Convert units on control settings */
    for (i = 1; i <= MaxControls; i++)
    {
        if ( (k = InpControl[i].Link) == 0) continue;
        if ( (j = InpControl[i].Node) > 0)
        {
            /* j = index of controlling node, and if           */
            /* j > Njuncs, then control is based on tank level */
            /* otherwise control is based on nodal pressure    */
            if (j > Njuncs)
                InpControl[i].Grade = Node[j].El + InpControl[i].Grade/Ucf[ELEV];
            else InpControl[i].Grade = Node[j].El + InpControl[i].Grade/Ucf[PRESSURE];
        }

        /* Convert units on valve settings */
        if (InpControl[i].Setting != MISSING) switch (Link[k].Type)
        {
        case PRV:
        case PSV:
        case PBV:
            InpControl[i].Setting /= Ucf[PRESSURE];
            break;
        case FCV:
            InpControl[i].Setting /= Ucf[FLOW];
        }
    }
}                       /*  End of convertunits  */

/************************ END OF INPUT1.C ************************/
