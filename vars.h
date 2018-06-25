/*
************************************************************************
            Global Variables for EPANET Program

AUTHOR:     M. Lopez-Ibanez (M. L\'opez-Ib\'a\~nez)
            Napier University, Edinburgh, UK.

$LastChangedBy: manu $ $Revision: 147 $
$LastChangedDate: 2007-11-20 15:14:06 +0100 (Tue, 20 Nov 2007) $

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
            6/24/02
AUTHOR:     L. Rossman
            US EPA - NRMRL

************************************************************************
*/

EXTERN FILE     *InFile;               /* Input file pointer           */
#if 0
                *OutFile,              /* Output file pointer          */
                *RptFile,              /* Report file pointer          */
                *HydFile,              /* Hydraulics file pointer      */
                *TmpOutFile;           /* Temporary file handle        */

EXTERN long     HydOffset,             /* Hydraulics file byte offset  */
                OutOffset1,            /* 1st output file byte offset  */
                OutOffset2;            /* 2nd output file byte offset  */
#endif
EXTERN char     
#if 0
                Msg[MAXMSG+1],         /* Text of output message       */
#endif
                InpFname[MAXFNAME+1],  /* Input file name              */
                Rpt1Fname[MAXFNAME+1], /* Primary report file name     */
                Rpt2Fname[MAXFNAME+1], /* Secondary report file name   */
                HydFname[MAXFNAME+1],  /* Hydraulics file name         */
                OutFname[MAXFNAME+1],  /* Binary output file name      */
                MapFname[MAXFNAME+1],  /* Map file name                */
                Title[MAXTITLE][MAXMSG+1], /* Problem title            */
                ChemName[MAXID+1],     /* Name of chemical             */
                ChemUnits[MAXID+1],    /* Units of chemical            */
                DefPatID[MAXID+1],     /* Default demand pattern ID    */
#if 0
                Hydflag,               /* Hydraulics flag              */
                Qualflag,              /* Water quality flag           */
#endif
                Unitsflag,             /* Unit system flag             */
                Flowflag,              /* Flow units flag              */
                Pressflag,             /* Pressure units flag          */
                Formflag,              /* Hydraulic formula flag       */
#if 0
                Rptflag,               /* Report flag                  */
                Summaryflag,           /* Report summary flag          */
                Messageflag,           /* Error/warning message flag   */
                Statflag,              /* Status report flag           */
                Energyflag,            /* Energy report flag           */
                Nodeflag,              /* Node report flag             */
                Linkflag,              /* Link report flag             */
#endif
                Tstatflag,             /* Time statistics flag         */
                Openflag;              /* Input processed flag         */
#if 0
                OpenQflag,             /* Quality system opened flag   */
                SaveQflag;             /* Quality results saved flag   */
#endif
EXTERN int      MaxNodes,              /* Node count from input file   */
                MaxLinks,              /* Link count from input file   */
                MaxJuncs,              /* Junction count               */
                MaxPipes,              /* Pipe count                   */
                MaxTanks,              /* Tank count                   */
                MaxPumps,              /* Pump count                   */
                MaxValves,             /* Valve count                  */
                MaxControls,           /* Control count                */
                MaxRules,              /* Rule count                   */
                MaxPats,               /* Pattern count                */
                MaxCurves,             /* Curve count                  */
                Nnodes,                /* Number of network nodes      */
                Ntanks,                // Number of tanks and reservoirs
                Nreservs,              // Number of reservoirs
                Njuncs,                /* Number of junction nodes     */
                Nlinks,                /* Number of network links      */
                Npipes,                /* Number of pipes              */
                Npumps,                /* Number of pumps              */
                Nvalves,               /* Number of valves             */
                Nrules,                /* Number of control rules      */
                Ncurves,               /* Number of data curves        */
                Nperiods,              /* Number of reporting periods  */
                DefPat,                /* Default demand pattern       */
                Epat,                  /* Energy cost time pattern     */
                MaxIter,               /* Max. hydraulic trials        */
                ExtraIter,             /* Extra hydraulic trials       */
#if 0
                TraceNode,             /* Source node for flow tracing */
                PageSize,              /* Lines/page in output report  */
#endif
                CheckFreq,             /* Hydraulics solver parameter  */
                MaxCheck;              /* Hydraulics solver parameter  */

EXTERN float    Ucf[MAXVAR],           /* Unit conversion factors      */
#if 0
                Ctol,                  /* Water quality tolerance      */
#endif
                Htol,                  /* Hydraulic head tolerance     */
                Qtol,                  /* Flow rate tolerance          */
                RQtol,                 /* Flow resistance tolerance    */
                Hexp,                  /* Exponent in headloss formula */
                Qexp,                  /* Exponent in orifice formula  */
                Dmult,                 /* Demand multiplier            */
                Hacc,                  /* Hydraulics solution accuracy */
#if 0
                BulkOrder,             /* Bulk flow reaction order     */
                WallOrder,             /* Pipe wall reaction order     */
                TankOrder,             /* Tank reaction order          */
                Kbulk,                 /* Global bulk reaction coeff.  */
                Kwall,                 /* Global wall reaction coeff.  */
                Climit,                /* Limiting potential quality   */
                Rfactor,               /* Roughness-reaction factor    */
                Diffus,                /* Diffusivity (sq ft/sec)      */
#endif
                Viscos,                /* Kin. viscosity (sq ft/sec)   */
                SpGrav,                /* Specific gravity             */
                Ecost,                 /* Base energy cost per kwh     */
                Dcost,                 /* Energy demand charge/kw/day  */
                Epump;                 /* Global pump efficiency       */
#if 0
                Wbulk,                 /* Avg. bulk reaction rate      */
                Wwall,                 /* Avg. wall reaction rate      */
                Wtank,                 /* Avg. tank reaction rate      */
                Wsource;               /* Avg. mass inflow             */
#endif
EXTERN long     Tstart,                /* Starting time of day (sec)   */
                Hstep,                 /* Nominal hyd. time step (sec) */
#if 0
                Qstep,                 /* Quality time step (sec)      */
#endif
                Pstep,                 /* Time pattern time step (sec) */
                Pstart,                /* Starting pattern time (sec)  */
                Rstep,                 /* Reporting time step (sec)    */
                Rstart,                /* Time when reporting starts   */
#if 0
                Qtime,                 /* Current quality time (sec)   */
#endif
                Rulestep,              /* Rule evaluation time step    */
                Dur;                   /* Duration of simulation (sec) */

EXTERN SField   Field[MAXVAR];         /* Output reporting fields      */

/* Array pointers not allocated and freed in same routine */
EXTERN STmplist *Patlist;              /* Temporary time pattern list  */
EXTERN STmplist *Curvelist;            /* Temporary list of curves     */
EXTERN Spattern *InpPattern;           /* Time patterns from input file */
EXTERN Scurve   *Curve;                /* Curve data                   */
EXTERN Snode    *Node;                 /* Node data                    */
EXTERN Slink    *Link;                 /* Link data                    */
EXTERN Stank    *Tank;                 /* Tank data                    */
EXTERN Spump    *Pump;                 /* Pump data                    */
EXTERN Svalve   *Valve;                /* Valve data                   */
EXTERN Scontrol *InpControl;           /* Control data from input file */
EXTERN HTtable  *Nht, *Lht;            /* Hash tables for ID labels    */


