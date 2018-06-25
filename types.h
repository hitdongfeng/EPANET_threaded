/*
***********************************************************************

TYPES.H -- Global constants and data types for EPANET program

AUTHOR:     M. Lopez-Ibanez (M. L\'opez-Ib\'a\~nez)
            Napier University, Edinburgh, UK.

$LastChangedBy: manu $ $Revision: 151 $
$LastChangedDate: 2007-12-12 17:26:42 +0100 (Wed, 12 Dec 2007) $

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
            9/7/00
            10/25/00
            3/1/01
            12/6/01
            6/24/02
AUTHOR:     L. Rossman
            US EPA - NRMRL

---------------------------------------------------------------------

References:

 [LA94a]  Kevin E. Lansey and K. Awumah. Optimal pump operations
          considering pump switches. Journal of Water Resources
          Planning and Management, ASCE, 120(1):17--35, January /
          February 1994.

**********************************************************************
*/
/*
-----------------------------
   Global Constants
-----------------------------
*/
/*** Updated 6/24/02 ***/
#define   CODEVERSION        20010

#define   MAGICNUMBER        516114521
#define   VERSION            200
#define   EOFMARK            0x1A  /* Use 0x04 for UNIX systems */
#define   MAXTITLE  3        /* Max. # title lines                     */
#define   MAXID     EN_MAX_ID_LEN  // Max. # characters in ID name
#define   MAXMSG    EN_MAX_MSG_LEN // Max. # characters in message text
#define   MAXLINE   255            // Max. # characters read from input line
#define   MAXFNAME  EN_MAX_FILENAME_LEN // Max. # characters in file name
#define   MAXTOKS   40       /* Max. items per line of input           */
#define   TZERO     1.E-4    /* Zero time tolerance                    */
#define   TRUE      1
#define   FALSE     0
#define   FULL      2
#define   BIG       1.E10
#define   TINY      1.E-6
#define   MISSING   -1.E10
#define   PI        3.141592654

/*** Updated 9/7/00 ***/
/* Various conversion factors */
#define   LperM3        1.0e-03
#define   M3perL        1.0e+03
#define   MperFT      0.3048 // exact value
#define   LperFT3     28.317 // 28.316846592
#define   M3perFT3    0.028317 // (LperFT3 * LperM3)
#define   GPMperCFS   448.831
#define   AFDperCFS   1.9837
#define   MGDperCFS   0.64632
#define   IMGDperCFS  0.5382
#define   LPSperCFS   LperFT3
#define   LPMperCFS   1699.0
#define   CMHperCFS   101.94
#define   CMDperCFS   2446.6
#define   MLDperCFS   2.4466
#define   PSIperFT    0.4333
#define   KPAperPSI   6.895
#define   KWperHP     0.7457
#define   SECperDAY   86400

#define   DIFFUS    1.3E-8   /* Diffusivity of chlorine                */
                             /* @ 20 deg C (sq ft/sec)                 */
#define   VISCOS    1.1E-5   /* Kinematic viscosity of water           */
                             /* @ 20 deg C (sq ft/sec)                 */

#define   SEPSTR    " \t\n\r"  /* Token separator characters */

/*
---------------------------------------------------------------------
   Macro to test for successful allocation of memory
---------------------------------------------------------------------
*/
#define  MEMCHECK(x)  (((x) == NULL) ? 101 : 0 )
#define  FREE(x)      (free((x)))

/*
---------------------------------------------------------------------
   Conversion macros to be used in place of functions
---------------------------------------------------------------------
*/
#define INT(x)   ((int)(x))                   /* integer portion of x  */
#define FRAC(x)  ((x)-(int)(x))               /* fractional part of x  */
#define ABS(x)   (((x)<0) ? -(x) : (x))       /* absolute value of x   */
#define MIN(x,y) (((x)<=(y)) ? (x) : (y))     /* minimum of x and y    */
#define MAX(x,y) (((x)>=(y)) ? (x) : (y))     /* maximum of x and y    */
#define ROUND(x) (((x)>=0) ? (int)((x)+.5) : (int)((x)-.5))
                                              /* round-off of x        */
#define MOD(x,y) ((x)%(y))                    /* x modulus y           */
#define SQR(x)   ((x)*(x))                    /* x-squared             */
#define SGN(x)   (((x)<0) ? (-1) : (1))       /* sign of x             */
#define UCHAR(x) (((x) >= 'a' && (x) <= 'z') ? ((x)&~32) : (x))
                                              /* uppercase char of x   */
/*
------------------------------------------------------
   Macro to evaluate function x with error checking
   (Fatal errors are numbered higher than 100)
------------------------------------------------------
*/
#define ERRCODE(x) (errcode = ((errcode>100) ? (errcode) : (x)))

/*
------------------------------------------------------
   Macro to find Pump index of Link[x]
   (Diameter = pump index for pump links)
------------------------------------------------------
*/
#define PUMPINDEX(x) (ROUND(Link[(x)].Diam))

/*
------------------------------------------------------
   Global Data Structures
------------------------------------------------------
*/
typedef  double        REAL;
typedef  long          INT4;

struct IDstring    /* Holds component ID labels */
{
    char ID[MAXID+1];
};

struct  Floatlist  /* Element of list of floats */
{
    float   value;
    struct  Floatlist *next;
};
typedef struct Floatlist SFloatlist;

struct  Tmplist    /* Element of temp list for Pattern & Curve data */
{
    int        i;
    char       ID[MAXID+1];
    SFloatlist *x;
    SFloatlist *y;
    struct     Tmplist  *next;
};
typedef struct Tmplist STmplist;

typedef struct        /* TIME PATTERN OBJECT */
{
    char  ID[MAXID+1]; /* Pattern ID       */
    int   Length;      /* Pattern length   */
    float *F;          /* Pattern factors  */
}  Spattern;

typedef struct        /* CURVE OBJECT */
{
    char  ID[MAXID+1]; /* Curve ID         */
    int   Type;        /* Curve type       */
    int   Npts;        /* Number of points */
    float *X;          /* X-values         */
    float *Y;          /* Y-values         */
}  Scurve;

struct Sdemand            /* DEMAND CATEGORY OBJECT */
{
    float Base;            /* Baseline demand  */
    int   Pat;             /* Pattern index    */
    struct Sdemand *next;  /* Next record      */
};
typedef struct Sdemand *Pdemand; /* Pointer to demand object */

#if 0
struct Ssource     /* WQ SOURCE OBJECT */
{
    /*int   Node;*/     /* Node index of source     */
    float C0;       /* Base concentration/mass  */
    int   Pat;      /* Pattern index            */
    float Smass;    /* Actual mass flow rate    */
    char  Type;     /* SourceType (see below)   */
};
typedef struct Ssource *Psource; /* Pointer to WQ source object */
#endif

typedef struct            /* NODE OBJECT */
{
    char   ID[MAXID+1];    /* Node ID          */
    float  El;             /* Elevation        */
    Pdemand D;             /* Demand pointer   */
    char    Ddefault;       /* Demand comes from [JUNCTION] section */
#if 0
    Psource S;             /* Source pointer   */
    float  C0;             /* Initial quality  */
#endif
    float  Ke;             /* Emitter coeff.   */
    char   Rpt;            /* Reporting flag   */
    float  x;              // x coordinate
    float  y;              // y coordinate
}  Snode;

typedef struct            /* LINK OBJECT */
{
    char   ID[MAXID+1];    /* Link ID           */
    int    N1;             /* Start node index  */
    int    N2;             /* End node index    */
    float  Diam;           // Diameter (or Pump index if Type == PUMP)
    float  Len;            /* Length            */
    float  Kc;             /* Roughness         */
    float  Km;             /* Minor loss coeff. */
#if 0
    float  Kb;             /* Bulk react. coeff */
    float  Kw;             /* Wall react. coeff */
#endif
    float  R;              /* Flow resistance   */
    char   Type;           /* Link type         */
    char   Stat;           /* Initial status    */
    char   Rpt;            /* Reporting flag    */
}  Slink;

typedef struct     /* TANK OBJECT */
{
    int   Node;     /* Node index of tank       */
    float A;        /* Tank area                */
    float Hmin;     /* Minimum water elev       */
    float Hmax;     /* Maximum water elev       */
    float H0;       /* Initial water elev       */
    float Vmin;     /* Minimum volume           */
    float Vmax;     /* Maximum volume           */
    float V0;       /* Initial volume           */
#if 0
    float Kb;       /* Reaction coeff. (1/days) */
    float C;        /* Concentration            */
#endif
    int   Pat;      /* Fixed grade time pattern */
    int   Vcurve;   /* Vol.- elev. curve index  */
#if 0
    char  MixModel; /* Type of mixing model     */
                    /* (see MixType below)      */
    float V1max;    /* Mixing compartment size  */
#endif
}  Stank;

typedef struct     /* PUMP OBJECT */
{
    int   Link;     /* Link index of pump          */
    int   Ptype;    /* Pump curve type             */
                    /* (see PumpType below)        */
    float Q0;       /* Initial flow                */
    float Qmax;     /* Maximum flow                */
    float Hmax;     /* Maximum head                */
    float H0;       /* Shutoff head                */
    float R;        /* Flow coeffic.               */
    float N;        /* Flow exponent               */
    int   Hcurve;   /* Head v. flow curve index    */
    int   Ecurve;   /* Effic. v. flow curve index  */
    int   Upat;     /* Utilization pattern index   */
    int   Epat;     /* Energy cost pattern index   */
    float Ecost;    /* Unit energy cost            */
}  Spump;

typedef struct     /* VALVE OBJECT */
{
    int   Link;     /* Link index of valve */
}  Svalve;

typedef struct     /* CONTROL STATEMENT */
{
    int   Link;     /* Link index         */
    int   Node;     /* Control node index */
    long  Time;     /* Control time       */
    float Grade;    /* Control grade      */
    float Setting;  /* New link setting   */
    char  Status;   /* New link status    */
    char  Type;     /* Control type       */
                    /* (see ControlType below) */
}  Scontrol;

struct  Sseg              /* PIPE SEGMENT record used */
{                         /*   for WQ routing         */
    float  v;              /* Segment volume      */
    float  c;              /* Water quality value */
    struct Sseg *prev;     /* Record for previous segment */
};
typedef struct Sseg *Pseg;    /* Pointer to pipe segment */

typedef struct            /* FIELD OBJECT of report table */
{
    char  Name[MAXID+1];   /* Name of reported variable  */
    char  Units[MAXID+1];  /* Units of reported variable */
    char  Enabled;         /* Enabled if in table        */
    int   Precision;       /* Number of decimal places   */
    float RptLim[2];       /* Lower/upper report limits  */
} SField;

typedef char ENclockstring_t[13]; /* String for Clock time (hrs:min:sec)  */

typedef struct {
    int   Nswitches; /* turning on a pump that was not operating in the
                       previous period [LA94a] */
    int   Schedule[24]; /* Time where a status change occurred */
    int   Schedule_idx; /* Index of the last status change */

    float Energy[6];  /* Energy usage statistics:  */
                      /* 0 = pump utilization      */
                      /* 1 = avg. efficiency       */
                      /* 2 = avg. kW/flow          */
                      /* 3 = avg. kwatts           */
                      /* 4 = peak kwatts           */
                      /* 5 = cost/day              */
} SimPump_t;

#include "smatrix.h" /* Fox smatrix_t */
/* A simulation instance. This is an internal object.  */
struct _ENsimulation_t {

    smatrix_t smatrix;

    int Haltflag; /* Flag used to halt taking further time steps */

    char
                Warnflag,              /* Warnings occurred            */
                SaveHflag,             /* Hydraul. results saved flag  */
                Saveflag;              /* Save results                 */

    char
                *S,                    /* Link status                  */
                *OldStat;              /* Previous link/tank status    */

    long
                Rtime,                 /* Next reporting time          */
                Htime;                 /* Current hyd. time (sec)      */


    float 
#if 0
                *C,                    /* Node actual quality          */
#endif
                *D,                    /* Node actual demand           */
                *E,                    /* Emitter flows                */
                *K,                    /* Link settings                */
#if 0
                *R,                    /* Pipe reaction rate           */
#endif
                *Q,                    /* Link flows                   */
                *V,                    /* Tank volumes              */
                *X,                    /* General purpose array        */

                Emax,                  /* Peak energy usage            */
                Dsystem,               /* Total system demand          */

                TotalSystemDemand,     /* Simulation total system demand  */
                TotalSystemInflow,     /* Simulation total system inflow  */
                TotalSystemLeakage;    /* Simulation total system leakage */

    int

             Ncontrols,             /* Number of simple controls    */
             Npats,                 /* Number of time patterns      */
             Nwarnings;             /* Number of warnings from the simulator */

    REAL     *H;                    /* Node heads                   */

    SimPump_t *Pump;           /* Pump characteristics during simulation */
    Spattern *Pattern;           /* Time patterns during simulation */
    Scontrol *Control;           /* Control data during simulation  */

/*
** NOTE: Hydraulic analysis of the pipe network at a given point in time
**       is done by repeatedly solving a linearized version of the
**       equations for conservation of flow & energy:
**
**           A*H = F
**
**       where H = vector of heads (unknowns) at each node,
**             F = vector of right-hand side coeffs.
**             A = square matrix of coeffs.
**       and both A and F are updated at each iteration until there is
**       negligible change in pipe flows.
**
**       Each row (or column) of A corresponds to a junction in the pipe
**       network. Each link (pipe, pump or valve) in the network has a
**       non-zero entry in the row-column of A that corresponds to its
**       end points. This results in A being symmetric and very sparse.
**       The following arrays are used to efficiently manage this sparsity:
*/
    REAL     *Aii,        /* Diagonal coeffs. of A               */
             *Aij,        /* Non-zero, off-diagonal coeffs. of A */
             *F;          /* Right hand side coeffs.             */
    float    *P,          /* Inverse headloss derivatives        */
             *Y;          /* Flow correction factors             */

};


/*
----------------------------------------------
   Global Enumeration Variables
----------------------------------------------
*/
 enum Hydtype                   /* Hydraulics solution option:         */
                {USE,           /*    use from previous run            */
                 SAVE,          /*    save after current run           */
                 SCRATCH};      /*    use temporary file               */

 enum QualType                  /* Water quality analysis option:      */
                {NONE,          /*    no quality analysis              */
                 CHEM,          /*    analyze a chemical               */
                 AGE,           /*    analyze water age                */
                 TRACE};        /*    trace % of flow from a source    */

 enum NodeType                  /* Type of node:                       */
                {JUNC,          /*    junction                         */
                 RESERV,        /*    reservoir                        */
                 TANK};         /*    tank                             */

 enum LinkType                  /* Type of link:                       */
                 {CV,           /*    pipe with check valve            */
                  PIPE,         /*    regular pipe                     */
                  PUMP,         /*    pump                             */
                  PRV,          /*    pressure reducing valve          */
                  PSV,          /*    pressure sustaining valve        */
                  PBV,          /*    pressure breaker valve           */
                  FCV,          /*    flow control valve               */
                  TCV,          /*    throttle control valve           */
                  GPV};         /*    general purpose valve            */

 enum CurveType                /* Type of curve:                       */
                 {V_CURVE,     /*    volume curve                      */
                  P_CURVE,     /*    pump curve                        */
                  E_CURVE,     /*    efficiency curve                  */
                  H_CURVE};    /*    head loss curve                   */

 enum PumpType                  /* Type of pump curve:                 */
                {CONST_HP,      /*    constant horsepower              */
                 POWER_FUNC,    /*    power function                   */
                 CUSTOM,        /*    user-defined custom curve        */
                 NOCURVE};

 enum SourceType                /* Type of source quality input        */
                {CONCEN,        /*    inflow concentration             */
                 MASS,          /*    mass inflow booster              */
                 SETPOINT,      /*    setpoint booster                 */
                 FLOWPACED};    /*    flow paced booster               */

 enum ControlType               /* Control condition type:             */
                {LOWLEVEL,      /*    act when grade below set level   */
                 HILEVEL,       /*    act when grade above set level   */
                 TIMER,         /*    act when set time reached        */
                 TIMEOFDAY};    /*    act when time of day occurs      */

 enum StatType                  /* Link/Tank status:                   */
                 {XHEAD,        /*   pump cannot deliver head (closed) */
                  TEMPCLOSED,   /*   temporarily closed                */
                  CLOSED,       /*   closed                            */
                  OPEN,         /*   open                              */
                  ACTIVE,       /*   valve active (partially open)     */
                  XFLOW,        /*   pump exceeds maximum flow         */
                  XFCV,         /*   FCV cannot supply flow            */
                  XPRESSURE,    /*   valve cannot supply pressure      */
                  FILLING,      /*   tank filling                      */
                  EMPTYING};    /*   tank emptying                     */

 enum FormType                  /* Head loss formula:                  */
                 {HW,           /*   Hazen-Williams                    */
                  DW,           /*   Darcy-Weisbach                    */
                  CM};          /*   Chezy-Manning                     */

 enum UnitsType                 /* Unit system:                        */
                 {US,           /*   US                                */
                  SI};          /*   SI (metric)                       */

 enum FlowUnitsType             /* Flow units:                         */
                 {CFS,          /*   cubic feet per second             */
                  GPM,          /*   gallons per minute                */
                  MGD,          /*   million gallons per day           */
                  IMGD,         /*   imperial million gal. per day     */
                  AFD,          /*   acre-feet per day                 */
                  LPS,          /*   liters per second                 */
                  LPM,          /*   liters per minute                 */
                  MLD,          /*   megaliters per day                */
                  CMH,          /*   cubic meters per hour             */
                  CMD};         /*   cubic meters per day              */

 enum PressUnitsType            /* Pressure units:                     */
                 {PSI,          /*   pounds per square inch            */
                  KPA,          /*   kiloPascals                       */
                  METERS};      /*   meters                            */

 enum RangeType                 /* Range limits:                       */
                 {LOW,          /*   lower limit                       */
                  HI,           /*   upper limit                       */
                  PREC};        /*   precision                         */

 enum MixType                   /* Tank mixing regimes                 */
                 {MIX1,         /*   1-compartment model               */
                  MIX2,         /*   2-compartment model               */
                  FIFO,         /*   First in, first out model         */
                  LIFO};        /*   Last in, first out model          */

 enum TstatType                 /* Time series statistics              */
                 {SERIES,       /*   none                              */
                  AVG,          /*   time-averages                     */
                  MIN,          /*   minimum values                    */
                  MAX,          /*   maximum values                    */
                  RANGE};       /*   max - min values                  */

#define MAXVAR   21             /* Max. # types of network variables   */
                                /* (equals # items enumed below)       */
 enum FieldType                 /* Network variables:                  */
                 {ELEV,         /*   nodal elevation                   */
                  DEMAND,       /*   nodal demand flow                 */
                  HEAD,         /*   nodal hydraulic head              */
                  PRESSURE,     /*   nodal pressure                    */
                  QUALITY,      /*   nodal water quality               */

                  LENGTH,       /*   link length                       */
                  DIAM,         /*   link diameter                     */
                  FLOW,         /*   link flow rate                    */
                  VELOCITY,     /*   link flow velocity                */
                  HEADLOSS,     /*   link head loss                    */
                  LINKQUAL,     /*   avg. water quality in link        */
                  STATUS,       /*   link status                       */
                  SETTING,      /*   pump/valve setting                */
                  REACTRATE,    /*   avg. reaction rate in link        */
                  FRICTION,     /*   link friction factor              */

                  POWER,        /*   pump power output                 */
                  TIME,         /*   simulation time                   */
                  VOLUME,       /*   tank volume                       */
                  CLOCKTIME,    /*   simulation time of day            */
                  FILLTIME,     /*   time to fill a tank               */
                  DRAINTIME};   /*   time to drain a tank              */

enum SectType    {_TITLE,_JUNCTIONS,_RESERVOIRS,_TANKS,_PIPES,_PUMPS,
                  _VALVES,_CONTROLS,_RULES,_DEMANDS,_SOURCES,_EMITTERS,
                  _PATTERNS,_CURVES,_QUALITY,_STATUS,_ROUGHNESS,_ENERGY,
                  _REACTIONS,_MIXING,_REPORT,_TIMES,_OPTIONS,
                  _COORDS,_VERTICES,_LABELS,_BACKDROP,_TAGS,_END};

enum HdrType                    /* Type of table heading   */
                 {STATHDR,      /*  Hydraulic Status       */
                  ENERHDR,      /*  Energy Usage           */
                  NODEHDR,      /*  Node Results           */
                  LINKHDR};     /*  Link Results           */
