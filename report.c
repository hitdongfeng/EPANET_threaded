/*
*********************************************************************

REPORT.C -- Reporting Routines for EPANET Program

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
            6/24/02
AUTHOR:     L. Rossman
            US EPA - NRMRL

This module contains various procedures (all beginning with
'write') that are called from other modules to write formatted
output to a report file.

It also contains function disconnected(), called from writehydwarn()
and writehyderr(), that checks if a hydraulic solution causes a
network to become disconnected.

The function writeline(S) is used throughout to write a
formatted string S to the report file.

********************************************************************
*/

#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <time.h>
#include "toolkit.h"
#include "hash.h"
#include "text.h"
#include "types.h"
#include "funcs.h"

#define  EXTERN  extern
#include "vars.h"

#define   MAXCOUNT 10     /* Max. # of disconnected nodes listed */
long      LineNum;        /* Current line number     */
long      PageNum;        /* Current page number     */
char      DateStamp[26];  /* Current date & time     */
char      Fprinterr;      /* File write error flag   */

/* Defined in enumstxt.h in EPANET.C */
extern char *NodeTxt[];
extern char *LinkTxt[];
extern char *StatTxt[];
extern char *TstatTxt[];
extern char *LogoTxt[];
extern char *RptFormTxt[];

typedef   float *Pfloat;
#if 0
static void      writenodetable(Pfloat *);
static void      writelinktable(Pfloat *);
#endif

static int   disconnected(ENsimulation_t);     /* Checks for disconnections  */
static void  marknodes(ENsimulation_t,         /* Identifies connected nodes */
                       int, int *, char *);
static void  getclosedlink(ENsimulation_t,     /* Finds a disconnecting link */
                           int, char *);

#if 0
int  writereport()
/*
**------------------------------------------------------
**   Input:   none
**   Output:  returns error code
**   Purpose: writes formatted output report to file
**
**   Calls strcomp() from the EPANET.C module.
**------------------------------------------------------
*/
{
    char  tflag;
    FILE  *tfile;
    int   errcode = 0;

    /* If no secondary report file specified then    */
    /* write formatted output to primary report file. */
    Fprinterr = FALSE;
    if (Rptflag && strlen(Rpt2Fname) == 0 && RptFile != NULL)
    {
        writecon(FMT17);
        writecon(Rpt1Fname);
        if (Energyflag) writeenergy();
        errcode = writeresults();
    }

    /* A secondary report file was specified */
    else if (strlen(Rpt2Fname) > 0)
    {

        /* If secondary report file has same name as either input */
        /* or primary report file then use primary report file.   */
        if (strcomp(Rpt2Fname,InpFname) ||
            strcomp(Rpt2Fname,Rpt1Fname))
        {
            writecon(FMT17);
            writecon(Rpt1Fname);
            if (Energyflag) writeenergy();
            errcode = writeresults();
        }

        /* Otherwise write report to secondary report file. */
        else
        {

            /* Try to open file */
            tfile = RptFile;
            tflag = Rptflag;
            if ((RptFile = fopen(Rpt2Fname,"wt")) == NULL)
            {
                RptFile = tfile;
                Rptflag = tflag;
                errcode = 303;
            }

            /* Write full formatted report to file */
            else
            {
                Rptflag = 1;
                writecon(FMT17);
                writecon(Rpt2Fname);
                writelogo();
                if (Summaryflag) writesummary();
                if (Energyflag)  writeenergy();
                errcode = writeresults();
                fclose(RptFile);
                RptFile = tfile;
                Rptflag = tflag;
            }
        }
    }

    /* Special error handler for write-to-file error */
    if (Fprinterr) errmsg(309);
    return(errcode);
}                        /* End of writereport */
#endif
#if 0
void  writelogo()
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  none
**   Purpose: writes program logo to report file.
**--------------------------------------------------------------
*/
{
    int    i;
    time_t timer;         /* time_t structure & functions time() & */
    /* ctime() are defined in time.h         */
    time(&timer);
    strcpy(DateStamp,ctime(&timer));
    PageNum = 1;
    LineNum = 2;
    fprintf(RptFile,FMT18);
    fprintf(RptFile,"%s",DateStamp);
    for (i=0; LogoTxt[i] != NULL; i++) writeline(LogoTxt[i]);
    writeline("");
}                        /* End of writelogo */
#endif

#if 0
void  writesummary()
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  none
**   Purpose: writes summary system information to report file
**--------------------------------------------------------------
*/
{
    char s[MAXFNAME+1];
    int  i;
    int  nres = 0;

    for (i=0; i<3; i++)
    {
        if (strlen(Title[i]) > 0)
        {
            sprintf(s,"%-.70s",Title[i]);
            writeline(s);
        }
    }
    writeline(" ");
    sprintf(s,FMT19,InpFname);
    writeline(s);
    sprintf(s,FMT20,Njuncs);
    writeline(s);
    for (i=1; i<=Ntanks; i++)
        if (Tank[i].A == 0.0) nres++;
    sprintf(s,FMT21a,nres);
    writeline(s);
    sprintf(s,FMT21b,Ntanks-nres);
    writeline(s);
    sprintf(s,FMT22,Npipes);
    writeline(s);
    sprintf(s,FMT23,Npumps);
    writeline(s);
    sprintf(s,FMT24,Nvalves);
    writeline(s);
    sprintf(s,FMT25,RptFormTxt[(int) Formflag]);
    writeline(s);
    sprintf(s,FMT26,Hstep*Ucf[TIME],Field[TIME].Units);
    writeline(s);
    sprintf(s,FMT27,Hacc);
    writeline(s);
    sprintf(s,FMT28,MaxIter);
    writeline(s);
    if (Qualflag == NONE || Dur == 0.0)
        sprintf(s,FMT29);
    else if (Qualflag == CHEM)
        sprintf(s,FMT30,ChemName);
    else if (Qualflag == TRACE)
        sprintf(s,FMT31,Node[TraceNode].ID);
    else if (Qualflag == AGE)
        sprintf(s,FMT32);
    writeline(s);
    if (Qualflag != NONE && Dur > 0)
    {
        sprintf(s,FMT33,(float)Qstep/60.0);
        writeline(s);
        sprintf(s,FMT34,Ctol*Ucf[QUALITY],Field[QUALITY].Units);
        writeline(s);
    }
    sprintf(s,FMT36,SpGrav);
    writeline(s);
    sprintf(s,FMT37a,Viscos/VISCOS);
    writeline(s);
    sprintf(s,FMT37b,Diffus/DIFFUS);
    writeline(s);
    sprintf(s,FMT38,Dmult);
    writeline(s);
    sprintf(s,FMT39,Dur*Ucf[TIME],Field[TIME].Units);
    writeline(s);
    if (Rptflag)
    {
        sprintf(s,FMT40);
        writeline(s);
        if (Nodeflag == 0)  writeline(FMT41);
        if (Nodeflag == 1)  writeline(FMT42);
        if (Nodeflag == 2)  writeline(FMT43);
        writelimits(DEMAND,QUALITY);
        if (Linkflag == 0)  writeline(FMT44);
        if (Linkflag == 1)  writeline(FMT45);
        if (Linkflag == 2)  writeline(FMT46);
        writelimits(DIAM,HEADLOSS);
    }
    writeline(" ");
}                        /* End of writesummary */
#endif

void  writehydstat(ENsimulation_t sim, int iter, float relerr)
/*
**--------------------------------------------------------------
**   Input:   iter   = # iterations to find hydraulic solution
**            relerr = convergence error in hydraulic solution
**   Output:  none
**   Purpose: writes hydraulic status report for solution found
**            at current time period to report file
**--------------------------------------------------------------
*/
{
    int    i,n;
    int    newstat;
#if 0
    char   s1[MAXLINE+1];

/*** Updated 6/24/02 ***/
    ENclockstring_t   atime;

    /* Display system status */
    clocktime (&atime, simulation->Htime);
    if (iter > 0)
    {
        if (relerr <= Hacc)
            sprintf(s1,FMT58,atime,iter);
        else
            sprintf(s1,FMT59,atime,iter,relerr);
        writeline(s1);
    }
#endif
    /*
      Display status changes for tanks.
      D[n] is net inflow to tank at node n.
      Old tank status is stored in OldStat[]
      at indexes Nlinks+1 to Nlinks+Ntanks.
    */
    for (i=1; i<=Ntanks; i++)
    {
        n = Tank[i].Node;
        if (ABS (sim->D[n]) < 0.001) newstat = CLOSED;
        else if (sim->D[n] >  0.0)  newstat = FILLING;
        else if (sim->D[n] <  0.0)  newstat = EMPTYING;
        else newstat = sim->OldStat[Nlinks+i];

        if (newstat != sim->OldStat[Nlinks+i])
        {
#if 0
            if (Tank[i].A > 0.0)
                sprintf(s1,FMT50,atime,Node[n].ID,StatTxt[newstat],
                        (simulation->H[n] - Node[n].El)
                        * Ucf[HEAD], Field[HEAD].Units);
            else sprintf(s1,FMT51,atime,Node[n].ID,StatTxt[newstat]);
            writeline(s1);
#endif
            sim->OldStat[Nlinks+i] = newstat;
        }
    }

    /* Display status changes for links */
    for (i=1; i<=Nlinks; i++)
    {
        if (sim->S[i] != sim->OldStat[i])
        {
#if 0
            if (simulation->Htime == 0)
                sprintf(s1,FMT52,atime,LinkTxt[(int) Link[i].Type],Link[i].ID,
                        StatTxt[(int) simulation->S[i]]);
            else sprintf(s1,FMT53,atime,LinkTxt[(int) Link[i].Type],Link[i].ID,
                         StatTxt[(int) simulation->OldStat[i]],
                         StatTxt[(int) simulation->S[i]]);
            writeline(s1);
#endif
            sim->OldStat[i] = sim->S[i];
        }
    }
#if 0
    writeline(" ");
#endif
}                        /* End of writehydstat */

#if 0
void  writeenergy()
/*
**-------------------------------------------------------------
**   Input:   none
**   Output:  none
**   Purpose: writes energy usage report to report file
**-------------------------------------------------------------
*/
{
    int   j;
    float csum;
    char  s[MAXLINE+1];
    if (Npumps == 0) return;
    writeline(" ");
    writeheader(ENERHDR,0);
    csum = 0.0;
    for (j=1; j<=Npumps; j++)
    {
        csum += Pump[j].Energy[5];
        if (LineNum == (long)PageSize) writeheader(ENERHDR,1);
        sprintf(s,"%-8s  %6.2f %6.2f %9.2f %9.2f %9.2f %9.2f",
                Link[Pump[j].Link].ID,Pump[j].Energy[0],Pump[j].Energy[1],
                Pump[j].Energy[2],Pump[j].Energy[3],Pump[j].Energy[4],
                Pump[j].Energy[5]);
        writeline(s);
    }
    fillstr(s,'-',63);
    writeline(s);

/*** Updated 6/24/02 ***/
    sprintf(s,FMT74,"",Emax*Dcost);
    writeline(s);
    sprintf(s,FMT75,"",csum+Emax*Dcost);
/*** End of update ***/

    writeline(s);
    writeline(" ");
}                       /* End of writeenergy */
#endif
#if 0
int  writeresults()
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  returns error code
**   Purpose: writes simulation results to report file
**--------------------------------------------------------------
*/
{
    Pfloat *x;                /* Array of pointers to floats */
    int    j,m,n,np,nnv,nlv;
    int    errcode = 0;

    /*
    **-----------------------------------------------------------
    **  NOTE:  The OutFile contains results for 4 node variables
    **         (demand, head, pressure, & quality) and 8 link
    **         variables (flow, velocity, headloss, quality,
    **         status, setting, reaction rate & friction factor)
    **         at each reporting time.
    **-----------------------------------------------------------
    */

    /* Return if no output file */
    if (OutFile == NULL) return(106);

    /* Return if no nodes or links selected for reporting */
    /* or if no node or link report variables enabled.    */
    if (!Nodeflag && !Linkflag) return(errcode);
    nnv = 0;
    for (j=ELEV; j<=QUALITY; j++) nnv += Field[j].Enabled;
    nlv = 0;
    for (j=LENGTH; j<=FRICTION; j++) nlv += Field[j].Enabled;
    if (nnv == 0 && nlv == 0) return(errcode);

    /* Allocate memory for output variables. */
    /* m = larger of # node variables & # link variables */
    /* n = larger of # nodes & # links */
    m = MAX( (QUALITY-DEMAND+1), (FRICTION-FLOW+1) );
    n = MAX( (Nnodes+1), (Nlinks+1));
    x = (Pfloat *) calloc(m, sizeof(Pfloat));
    ERRCODE( MEMCHECK(x) );
    if (errcode) return(errcode);
    for (j=0; j<m; j++)
    {
        x[j] = (float *) calloc(n, sizeof(float));
        ERRCODE( MEMCHECK(x[j]) );
    }
    if (errcode) return(errcode);

    /* Re-position output file & initialize report time. */
    fseek(OutFile,OutOffset2,SEEK_SET);
    simulation->Htime = Rstart;

    /* For each reporting time: */
    for (np=1; np<=Nperiods; np++)
    {

        /* Read in node results & write node table. */
        /* (Remember to offset x[j] by 1 because array is zero-based). */
        for (j=DEMAND; j<=QUALITY; j++)
            fread((x[j-DEMAND])+1,sizeof(float),Nnodes,OutFile);
        if (nnv > 0 && Nodeflag > 0) writenodetable(x);

        /* Read in link results & write link table. */
        for (j=FLOW; j<=FRICTION; j++)
            fread((x[j-FLOW])+1,sizeof(float),Nlinks,OutFile);
        if (nlv > 0 && Linkflag > 0) writelinktable(x);
        simulation->Htime += Rstep;
    }

    /* Free allocated memory */
    for (j=0; j<m; j++) free(x[j]);
    free(x);
    return(errcode);
}                        /* End of writereport */
#endif

#if 0
static void  writenodetable(Pfloat *x)
/*
**---------------------------------------------------------------
**   Input:   x = pointer to node results for current time
**   Output:  none
**   Purpose: writes node results for current time to report file
**---------------------------------------------------------------
*/
{
    int i,j;
    char s[MAXLINE+1],s1[16];
    float y[MAXVAR];

    /* Write table header */
    writeheader(NODEHDR,0);

    /* For each node: */
    for (i=1; i<=Nnodes; i++)
    {

        /* Place results for each node variable in y */
        y[ELEV] = Node[i].El*Ucf[ELEV];
        for (j=DEMAND; j<=QUALITY; j++) y[j] = *((x[j-DEMAND])+i);

        /* Check if node gets reported on */
        if ((Nodeflag == 1 || Node[i].Rpt) && checklimits(y,ELEV,QUALITY))
        {

            /* Check if new page needed */
            if (LineNum == (long)PageSize) writeheader(NODEHDR,1);

            /* Add node ID and each reported field to string s */
            sprintf(s,"%-15s",Node[i].ID);
            for (j=ELEV; j<=QUALITY; j++)
            {
                if (Field[j].Enabled == TRUE)
                {

/*** Updated 6/24/02 ***/
                    if (fabs(y[j]) > 1.e6) sprintf(s1, "%10.2e", y[j]);
                    else                   sprintf(s1, "%10.*f", Field[j].Precision, y[j]);
/*** End of update ***/

                    strcat(s, s1);
                }
            }

            /* Note if node is a reservoir/tank */
            if (i > Njuncs)
            {
                strcat(s, "  ");
                strcat(s, NodeTxt[getnodetype(i)]);
            }

            /* Write results for node */
            writeline(s);
        }
    }
    writeline(" ");
}
#endif
#if 0
static void  writelinktable(Pfloat *x)
/*
**---------------------------------------------------------------
**   Input:   x = pointer to link results for current time
**   Output:  none
**   Purpose: writes link results for current time to report file
**---------------------------------------------------------------
*/
{
    int i,j,k;
    char s[MAXLINE+1],s1[16];
    float y[MAXVAR];

    /* Write table header */
    writeheader(LINKHDR,0);

    /* For each link: */
    for (i=1; i<=Nlinks; i++)
    {

        /* Place results for each link variable in y */
        y[LENGTH] = Link[i].Len*Ucf[LENGTH];
        y[DIAM] = Link[i].Diam*Ucf[DIAM];
        for (j=FLOW; j<=FRICTION; j++) y[j] = *((x[j-FLOW])+i);

        /* Check if link gets reported on */
        if ((Linkflag == 1 || Link[i].Rpt) && checklimits(y,DIAM,FRICTION))
        {

            /* Check if new page needed */
            if (LineNum == (long)PageSize) writeheader(LINKHDR,1);

            /* Add link ID and each reported field to string s */
            sprintf(s,"%-15s",Link[i].ID);
            for (j=LENGTH; j<=FRICTION; j++)
            {
                if (Field[j].Enabled == TRUE)
                {
                    if (j == STATUS)
                    {
                        if      (y[j] <= CLOSED) k = CLOSED;
                        else if (y[j] == ACTIVE) k = ACTIVE;
                        else                     k = OPEN;
                        sprintf(s1, "%10s", StatTxt[k]);
                    }

/*** Updated 6/24/02 ***/
                    else
                    {
                        if (fabs(y[j]) > 1.e6) sprintf(s1, "%10.2e", y[j]);
                        else                   sprintf(s1, "%10.*f", Field[j].Precision, y[j]);
                    }
/*** End of update ***/

                    strcat(s, s1);
                }
            }

            /* Note if link is a pump or valve */
            if ( (j = Link[i].Type) > PIPE)
            {
                strcat(s, "  ");
                strcat(s, LinkTxt[j]);
            }

            /* Write results for link */
            writeline(s);
        }
    }
    writeline(" ");
}
#endif

#if 0
void  writeheader(int type, int contin)
/*
**--------------------------------------------------------------
**   Input:   type   = table type
**            contin = table continuation flag
**   Output:  none
**   Purpose: writes column headings for output report tables
**--------------------------------------------------------------
*/
{
    char   s[MAXLINE+1],s1[MAXLINE+1],s2[MAXLINE+1],s3[MAXLINE+1];
    int    i,n;

    ENclockstring_t atime;

    /* Move to next page if < 11 lines remain on current page. */
    if (Rptflag && LineNum+11 > (long)PageSize)
    {
        while (LineNum < (long)PageSize) writeline(" ");
    }
    writeline(" ");

    /* Hydraulic Status Table */
    if (type == STATHDR)
    {
        sprintf(s,FMT49);
        if (contin) strcat(s,t_CONTINUED);
        writeline(s);
        fillstr(s,'-',70);
        writeline(s);
    }

    /* Energy Usage Table */
    if (type == ENERHDR)
    {
        if (Unitsflag == SI) strcpy(s1,t_perM3);
        else                 strcpy(s1,t_perMGAL);
        sprintf(s,FMT71);
        if (contin) strcat(s,t_CONTINUED);
        writeline(s);
        fillstr(s,'-',63);
        writeline(s);
        sprintf(s,FMT72);
        writeline(s);
        sprintf(s,FMT73,s1);
        writeline(s);
        fillstr(s,'-',63);
        writeline(s);
    }

    /* Node Results Table */
    if (type == NODEHDR)
    {
        if      (Tstatflag == RANGE)  sprintf(s,FMT76,t_DIFFER);
        else if (Tstatflag != SERIES) sprintf(s,FMT76,TstatTxt[(int) Tstatflag]);
        else if (Dur == 0)            sprintf(s,FMT77);
        else                          
            sprintf (s, FMT78, clocktime (&atime, simulation->Htime));

        if (contin) strcat(s,t_CONTINUED);
        writeline(s);
        n = 15;
        sprintf(s2,"%15s","");
        strcpy(s,t_NODEID);
        sprintf(s3,"%-15s",s);
        for (i=ELEV; i<QUALITY; i++) if (Field[i].Enabled == TRUE)
        {
            n += 10;
            sprintf(s,"%10s",Field[i].Name);
            strcat(s2,s);
            sprintf(s,"%10s",Field[i].Units);
            strcat(s3,s);
        }
        if (Field[QUALITY].Enabled == TRUE)
        {
            n += 10;
            sprintf(s,"%10s",ChemName);
            strcat(s2,s);
            sprintf(s,"%10s",ChemUnits);
            strcat(s3,s);
        }
        fillstr(s1,'-',n);
        writeline(s1);
        writeline(s2);
        writeline(s3);
        writeline(s1);
    }

    /* Link Results Table */
    if (type == LINKHDR)
    {
        if      (Tstatflag == RANGE)  sprintf(s,FMT79,t_DIFFER);
        else if (Tstatflag != SERIES) sprintf(s,FMT79,TstatTxt[(int) Tstatflag]);
        else if (Dur == 0)            sprintf(s,FMT80);
        else                          
            sprintf (s, FMT81, clocktime (&atime, simulation->Htime));

        if (contin) strcat(s,t_CONTINUED);
        writeline(s);
        n = 15;
        sprintf(s2,"%15s","");
        strcpy(s,t_LINKID);
        sprintf(s3,"%-15s",s);
        for (i=LENGTH; i<=FRICTION; i++) if (Field[i].Enabled == TRUE)
        {
            n += 10;
            sprintf(s,"%10s",Field[i].Name);
            strcat(s2,s);
            sprintf(s,"%10s",Field[i].Units);
            strcat(s3,s);
        }
        fillstr(s1,'-',n);
        writeline(s1);
        writeline(s2);
        writeline(s3);
        writeline(s1);
    }
}                        /* End of writeheader */
#endif

#if 0
void  writeline(char *s)
/*
**--------------------------------------------------------------
**   Input:   *s = text string
**   Output:  none
**   Purpose: writes a line of output to report file
**--------------------------------------------------------------
*/
{
    if (RptFile == NULL) return;
    if (Rptflag)
    {
        if (LineNum == (long)PageSize)
        {
            PageNum++;
            if (fprintf(RptFile,FMT82,PageNum,Title[0]) == EOF)
                Fprinterr = TRUE;
            LineNum = 3;
        }
    }
    if (fprintf(RptFile,"\n  %s",s) == EOF) Fprinterr = TRUE;
    LineNum++;
}                        /* End of writeline */
#endif

#if 0
void  writerelerr(int iter, float relerr)
/*
**-----------------------------------------------------------------
**   Input:   iter   = current iteration of hydraulic solution
**            relerr = current convergence error
**   Output:  none
**   Purpose: writes out convergence status of hydraulic solution
**-----------------------------------------------------------------
*/
{
    ENclockstring_t atime;

    if (iter == 0)
    {
        sprintf (Msg, FMT64, clocktime (&atime, simulation->Htime));
        writeline(Msg);
    }
    else
    {
        sprintf(Msg,FMT65,iter,relerr);
        writeline(Msg);
    }
}                        /* End of writerelerr */
#endif

#if 0
void  writestatchange(int k, char s1, char s2)
/*
**--------------------------------------------------------------
**   Input:   k  = link index
**            s1 = old link status
**            s2 = new link status
**   Output:  none
**   Purpose: writes change in link status to output report
**--------------------------------------------------------------
*/
{
    int   j1,j2;
    float setting;

/* We have a pump/valve setting change instead of a status change */
    if (s1 == s2)
    {

/*** Updated 10/25/00 ***/
        setting = simulation->K[k];   /* Link[k].Kc; */

        switch (Link[k].Type)
        {
        case PRV:
        case PSV:
        case PBV: setting *= Ucf[PRESSURE]; break;
        case FCV: setting *= Ucf[FLOW];
        }
        sprintf(Msg,FMT56,LinkTxt[(int) Link[k].Type],Link[k].ID,setting);
        writeline(Msg);
        return;
    }

/* We have a status change. Write the old & new status types. */
    if      (s1 == ACTIVE) j1 = ACTIVE;
    else if (s1 <= CLOSED) j1 = CLOSED;
    else                   j1 = OPEN;
    if      (s2 == ACTIVE) j2 = ACTIVE;
    else if (s2 <= CLOSED) j2 = CLOSED;
    else                   j2 = OPEN;
    if (j1 != j2)
    {
        sprintf(Msg,FMT57,LinkTxt[(int) Link[k].Type],Link[k].ID,
                StatTxt[j1],StatTxt[j2]);
        writeline(Msg);
    }
}                        /* End of writestatchange */
#endif

#if 0
void writecontrolaction(int k, int i)
/*
  ----------------------------------------------------------------
  **   Input:   k  = link index
  **            i  = control index
  **   Output:  none
  **   Purpose: writes control action taken to status report
  **--------------------------------------------------------------
  */
{
    int n;
    ENclockstring_t atime;

    switch (Control[i].Type)
    {
    case LOWLEVEL:
    case HILEVEL:
        n = Control[i].Node;
        sprintf (Msg, FMT54, clocktime (&atime, simulation->Htime),
                 LinkTxt[(int) Link[k].Type], Link[k].ID,
                 NodeTxt[getnodetype(n)], Node[n].ID);
        break;
    case TIMER:
    case TIMEOFDAY:
        sprintf (Msg, FMT55, clocktime (&atime, simulation->Htime),
                 LinkTxt[(int) Link[k].Type], Link[k].ID);
        break;
    default: return;
    }
    writeline(Msg);
}
#endif

#if 0
void writeruleaction(int k, char *ruleID)
/*
**--------------------------------------------------------------
**   Input:   k  = link index
**            *ruleID  = rule ID
**   Output:  none
**   Purpose: writes rule action taken to status report
**--------------------------------------------------------------
*/
{
    ENclockstring_t atime;
    sprintf (Msg, FMT63, clocktime (&atime, simulation->Htime), 
             LinkTxt[(int) Link[k].Type], Link[k].ID, ruleID);
    writeline(Msg);
}
#endif


int  writehydwarn(ENsimulation_t sim, int iter, float relerr)
/*
**--------------------------------------------------------------
**   Input:   iter = # iterations to find hydraulic solution
**   Output:  warning flag code
**   Purpose: writes hydraulic warning message to report file
**
**   Note: Warning conditions checked in following order:
**         1. System balanced but unstable
**         2. Negative pressures
**         3. FCV cannot supply flow or PRV/PSV cannot maintain pressure
**         4. Pump out of range
**         5. Network disconnected
**         6. System unbalanced
**--------------------------------------------------------------
*/
{
    int  i,j;
    char flag = 0;
#if 0
    ENclockstring_t atime;
#endif

    /* Check if system unstable */
    if (iter > MaxIter && relerr <= Hacc)
    {
#if 0
        sprintf (Msg, WARN02, clocktime (&atime, sim->Htime));
        if (Messageflag) writeline(Msg);
#endif
        flag = 2;
        sim->Nwarnings++;
    }

    /* Check for negative pressures */
    for (i=1; i<=Njuncs; i++)
    {
        if (sim->H[i] < Node[i].El && sim->D[i] > 0.0)
        {
#if 0
            sprintf (Msg, WARN06, clocktime (&atime, sim->Htime));
            if (Messageflag) writeline(Msg);
#endif
            flag = 6;
            sim->Nwarnings++;
            break;
        }
    }

    /* Check for abnormal valve condition */
    for (i=1; i<=Nvalves; i++)
    {
        j = Valve[i].Link;
        if (sim->S[j] >= XFCV)
        {
#if 0
            sprintf (Msg, WARN05, LinkTxt[(int) Link[j].Type],Link[j].ID,
                     StatTxt[(int) sim->S[j]], 
                     clocktime (&atime, sim->Htime));
            if (Messageflag) writeline(Msg);
#endif
            flag = 5;
            sim->Nwarnings++;
        }
    }

    /* Check for abnormal pump condition */
    for (i=1; i<=Npumps; i++)
    {
        j = Pump[i].Link;
        if (sim->S[j] == XHEAD || sim->S[j] == XFLOW)
        {
#if 0
            sprintf(Msg,WARN04,Link[j].ID,StatTxt[(int) sim->S[j]],
                    clocktime (&atime, sim->Htime));
            if (Messageflag) writeline(Msg);
#endif
            flag = 4;
            sim->Nwarnings++;
        }
    }

    /* Check if system is unbalanced */
    if (iter > MaxIter && relerr > Hacc)
    {
#if 0
        sprintf (Msg, WARN01, clocktime (&atime, sim->Htime));
        if (ExtraIter == -1) strcat(Msg,t_HALTED);
        if (Messageflag) writeline(Msg);
#endif
        flag = 1;
        sim->Nwarnings++;
    }

    /* Check for disconnected network */
    /* & update global warning flag   */
    if (flag > 0)
    {
        disconnected (sim);
        sim->Warnflag = flag;
    }
    return(flag);
}                        /* End of writehydwarn */


void  writehyderr(ENsimulation_t sim, int errnode)
/*
**-----------------------------------------------------------
**   Input:   none
**   Output:  none
**   Purpose: outputs status & checks connectivity when
**            network hydraulic equations cannot be solved.
**-----------------------------------------------------------
*/
{
#if 0
    ENclockstring_t atime;
    sprintf (Msg, FMT62, clocktime (&atime, simulation->Htime),
             Node[errnode].ID);
    if (Messageflag) writeline(Msg);
#endif
    writehydstat (sim, 0,0);
    disconnected (sim);
}                        /* End of writehyderr */


static int  disconnected(ENsimulation_t sim)
/*
**-------------------------------------------------------------------
**   Input:   None
**   Output:  Returns number of disconnected nodes
**   Purpose: Tests current hydraulic solution to see if any closed
**            links have caused the network to become disconnected.
**-------------------------------------------------------------------
*/
{
    int  i, j;
    int  count, mcount;
    int  errcode = 0;
    int  *nodelist;
    char *marked;
#if 0
    ENclockstring_t atime;
#endif

    /* Allocate memory for node list & marked list */
    nodelist = (int *)  calloc(Nnodes+1,sizeof(int));
    marked   = (char *) calloc(Nnodes+1,sizeof(char));
    ERRCODE(MEMCHECK(nodelist));
    ERRCODE(MEMCHECK(marked));
    if (errcode) return(0);

    /* Place tanks on node list and marked list */
    for (i=1; i<=Ntanks; i++)
    {
        j = Njuncs + i;
        nodelist[i] = j;
        marked[j] = 1;
    }

    /* Place junctions with negative demands on the lists */
    mcount = Ntanks;
    for (i=1; i<=Njuncs; i++)
    {
        if (sim->D[i] < 0.0)
        {
            mcount++;
            nodelist[mcount] = i;
            marked[i] = 1;
        }
    }

    /* Mark all nodes that can be connected to tanks */
    /* and count number of nodes remaining unmarked. */
    marknodes (sim, mcount, nodelist, marked);
    j = 0;
    count = 0;
    for (i=1; i<=Njuncs; i++)
    {
        if (!marked[i] && sim->D[i] != 0.0)  /* Skip if no demand */
        {
            count++;
#if 0
            if (count <= MAXCOUNT && Messageflag)
            {
                sprintf (Msg, WARN03a, Node[i].ID, 
                         clocktime (&atime, simulation->Htime));
                writeline(Msg);
            }
#endif
            j = i;                       /* Last unmarked node */
        }
    }

    /* Report number of unmarked nodes and find closed link */
    /* on path from node j back to a tank.                  */
    if (count > 0) /* FIXME: && Messageflag) */
    {
        if (count > MAXCOUNT)
        {
#if 0
            sprintf (Msg, WARN03b, count - MAXCOUNT, 
                     clocktime (&atime, simulation->Htime));
            writeline(Msg);
#endif
        }
        getclosedlink (sim, j, marked);
    }

    /* Free allocated memory */
    free(nodelist);
    free(marked);
    return(count);
}                   /* End of disconnected() */


static void  marknodes(ENsimulation_t sim, int m, int *nodelist, char *marked)
/*
**----------------------------------------------------------------
**   Input:   m = number of source nodes
**            nodelist[] = list of nodes to be traced from
**            marked[]   = TRUE if node connected to source
**   Output:  None.
**   Purpose: Marks all junction nodes connected to tanks.
**----------------------------------------------------------------
*/
{
    int   i, j, k, n;
    Padjlist alink;

    /* Scan each successive entry of node list */
    n = 1;
    while (n <= m )
    {

        /* Scan all nodes connected to current node */
        i = nodelist[n];
        for (alink = smatrix_adjlist_item (sim->smatrix, i); alink != NULL; alink = alink->next)
        {

            /* Get indexes of connecting link and node */
            k = alink->link;
            j = alink->node;
            if (marked[j]) continue;

            /* Check if valve connection is in correct direction */
            switch (Link[k].Type)
            {
            case CV:
            case PRV:
            case PSV: if (j == Link[k].N1) continue;
            }

            /* Mark connection node if link not closed */
            if (sim->S[k] > CLOSED)
            {
                marked[j] = 1;
                m++;
                nodelist[m] = j;
            }
        }
        n++;
    }
}                   /* End of marknodes() */


static void getclosedlink(ENsimulation_t sim, int i, char *marked)
/*
**----------------------------------------------------------------
**   Input:   i = junction index
**            marked[] = marks nodes already examined
**   Output:  None.
**   Purpose: Determines if a closed link connects to junction i.
**----------------------------------------------------------------
*/
{
    int j,k;
    Padjlist alink;
    marked[i] = 2;
    for (alink = smatrix_adjlist_item (sim->smatrix, i); 
         alink != NULL; alink = alink->next)
    {
        k = alink->link;
        j = alink->node;
        if (marked[j] == 2) continue;
        if (marked[j] == 1)
        {
#if 0
            sprintf(Msg, WARN03c, Link[k].ID);
            writeline(Msg);
#endif
            return;
        }
        else getclosedlink (sim, j, marked);
    }
}

#if 0
void  writelimits(int j1, int j2)
/*
**--------------------------------------------------------------
**   Input:   j1 = index of first output variable
**            j2 = index of last output variable
**   Output:  none
**   Purpose: writes reporting criteria to output report
**--------------------------------------------------------------
*/
{
    int  j;
    for (j=j1; j<=j2; j++)
    {
        if (Field[j].RptLim[LOW] < BIG)
        {
            sprintf(Msg,FMT47,
                    Field[j].Name,Field[j].RptLim[LOW],Field[j].Units);
            writeline(Msg);
        }
        if (Field[j].RptLim[HI] > -BIG)
        {
            sprintf(Msg,FMT48,
                    Field[j].Name,Field[j].RptLim[HI],Field[j].Units);
            writeline(Msg);
        }
    }
}                        /* End of writelimits */
#endif

int  checklimits(float *y, int j1, int j2)
/*
**--------------------------------------------------------------
**   Input:   *y = array of output results
**            j1 = index of first output variable
**            j2 = index of last output variable
**   Output:  returns 1 if criteria met, 0 otherwise
**   Purpose: checks if output reporting criteria is met
**--------------------------------------------------------------
*/
{
    int j;
    for (j=j1; j<=j2; j++)
    {
        if (y[j] > Field[j].RptLim[LOW]
            ||  y[j] < Field[j].RptLim[HI]) return(0);
    }
    return(1);
}                        /* End of checklim */

#if 0
void writetime(char *fmt)
/*
**----------------------------------------------------------------
**   Input:   fmt = format string
**   Output:  none
**   Purpose: writes starting/ending time of a run to report file
**----------------------------------------------------------------
*/
{
    time_t timer;
    time(&timer);
    sprintf(Msg, fmt, ctime(&timer));
    writeline(Msg);
}
#endif

char *clocktime(ENclockstring_t *atime, long seconds)
/*
**--------------------------------------------------------------
**   Input:   seconds = time in seconds
**   Output:  atime = time in hrs:min
**            (returns pointer to atime)
**   Purpose: converts time in seconds to hours:minutes format
**--------------------------------------------------------------
*/
{
/*** Updated 6/24/02 ***/
    int h,m,s;
    h = seconds/3600;
    m = (seconds % 3600) / 60;
    s = seconds - 3600*h - 60*m;
    sprintf (*atime, "%01d:%02d:%02d", h,m,s);
    return *atime;
}                        /* End of clocktime */


char *fillstr(char *s, char ch, int n)
/*
**---------------------------------------------------------
**  Fills n bytes of s to character ch.
**  NOTE: does not check for overwriting s.
**---------------------------------------------------------
*/
{
    int i;
    for (i=0; i<=n; i++) s[i] = ch;
    s[n+1] = '\0';
    return(s);
}


int  getnodetype(int i)
/*
**---------------------------------------------------------
**  Determines type of node with index i
**  (junction = 0, reservoir = 1, tank = 2).
**---------------------------------------------------------
*/
{
    if (i <= Njuncs) return(0);
    if (Tank[i-Njuncs].A == 0.0) return(1);
    return(2);
}


/********************* END OF REPORT.C ********************/
