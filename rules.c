/*
**********************************************************************

RULES.C -- Rule processor module for EPANET

AUTHOR:     M. Lopez-Ibanez (M. L\'opez-Ib\'a\~nez)
            Napier University, Edinburgh, UK.

$LastChangedBy: manu $ $Revision: 173 $
$LastChangedDate: 2008-09-10 15:30:14 +0200 (Wed, 10 Sep 2008) $

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
AUTHOR:     L. Rossman
            US EPA - NRMRL

  The entry points for this module are:
     initrules()  -- called from ENopen() in EPANET.C
     addrule()    -- called from netsize() in INPUT2.C
     allocrules() -- called from allocdata() in EPANET.C
     ruledata()   -- called from newline() in INPUT2.C
     freerules()  -- called from freedata() in EPANET.C
     checkrules() -- called from ruletimestep() in HYDRAUL.C

**********************************************************************
*/
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include "toolkit.h"
#include "hash.h"
#include "text.h"
#include "types.h"
#include "funcs.h"
#define  EXTERN  extern
#include "vars.h"

struct      Premise         /* Rule Premise Clause */
{
    int      logop;          /* Logical operator */
    int      object;         /* Node or link */
    int      index;          /* Object's index */
    int      variable;       /* Pressure, flow, etc. */
    int      relop;          /* Relational operator */
    int      status;         /* Variable's status */
    float    value;          /* Variable's value */
    struct   Premise *next;
};

struct     Action           /* Rule Action Clause */
{
    int     link;            /* Link index */
    int     status;          /* Link's status */
    float   setting;         /* Link's setting */
    struct  Action *next;
};

struct      aRule           /* Control Rule Structure */
{
    char     label[MAXID+1];    /* Rule character label */
    float    priority;          /* Priority level */
    struct   Premise  *Pchain;  /* Linked list of premises */
    struct   Action   *Tchain;  /* Linked list of actions if true */
    struct   Action   *Fchain;  /* Linked list of actions if false */
    struct   aRule    *next;
};

struct      ActItem         /* Action list item */
{
    int      ruleindex;        /* Index of rule action belongs to */
    struct   Action   *action; /* An action structure */
    struct   ActItem  *next;
};

struct  aRule *Rule;        /* Array of rules */
int     RuleState;          /* State of rule interpreter */
struct  Premise *Plast;     /* Previous premise clause */

enum    Rulewords      {r_RULE,r_IF,r_AND,r_OR,r_THEN,r_ELSE,r_PRIORITY,r_ERROR};
char    *Ruleword[]  = {w_RULE,w_IF,w_AND,w_OR,w_THEN,w_ELSE,w_PRIORITY,NULL};

enum    Varwords       {r_DEMAND, r_HEAD, r_GRADE, r_LEVEL, r_PRESSURE,
                        r_FLOW, r_STATUS, r_SETTING, r_POWER, r_TIME,
                        r_CLOCKTIME, r_FILLTIME, r_DRAINTIME};
char    *Varword[]   = {w_DEMAND, w_HEAD, w_GRADE, w_LEVEL, w_PRESSURE,
                        w_FLOW, w_STATUS, w_SETTING, w_POWER,w_TIME,
                        w_CLOCKTIME,w_FILLTIME,w_DRAINTIME, NULL};

enum    Objects        {r_JUNC,r_RESERV,r_TANK,r_PIPE,r_PUMP,r_VALVE,
                        r_NODE,r_LINK,r_SYSTEM};
char    *Object[]    = {w_JUNC,w_RESERV,w_TANK,w_PIPE,w_PUMP,w_VALVE,
                        w_NODE,w_LINK,w_SYSTEM,NULL};

/* NOTE: place "<=" & ">=" before "<" & ">" so that findmatch() works correctly. */
enum    Operators      { EQ, NE,  LE,  GE,  LT, GT, IS,  NOT,  BELOW,  ABOVE};
char    *Operator[]  = {"=","<>","<=",">=","<",">",w_IS,w_NOT,w_BELOW,w_ABOVE,NULL};

enum    Values         {IS_NUMBER,IS_OPEN,IS_CLOSED,IS_ACTIVE};
char    *Value[]     = {"XXXX",   w_OPEN, w_CLOSED, w_ACTIVE,NULL};

/* External variables declared in INPUT2.C */
extern char      *Tok[MAXTOKS];
extern int       Ntokens;

/*
**   Local function prototypes are defined here and not in FUNCS.H
**   because some of them utilize the Premise and Action structures
**   defined locally in this module.
*/
static void newrule(void);
static int  newpremise(int);
static int  newaction(void);
static int  newpriority(void);
static int  evalpremises(ENsimulation_t, int, long);
static void updateactlist(struct ActItem **, int, struct Action *);
static int  checkaction(struct ActItem *, int, struct Action *);
static int  checkpremise(ENsimulation_t, struct Premise *, long dt);
static int  checktime(struct Premise *, long Htime, long dt);
static int  checkstatus(struct Premise *, char *S);
static int  checkvalue(ENsimulation_t, struct Premise *);
static int  takeactions(ENsimulation_t sim, struct ActItem *item);
static void clearactlist(struct  ActItem *a);
#if 0
static void ruleerrmsg(int);
#endif

void initrules()
/*
**--------------------------------------------------------------
**    Initializes rule base.
**    Called by ENopen() in EPANET.C module
**--------------------------------------------------------------
*/
{
    RuleState = r_PRIORITY;
    Rule = NULL;
}


void addrule(char *tok)
/*
**--------------------------------------------------------------
**    Updates rule count if RULE keyword found in line of input.
**    Called by netsize() in INPUT2.C module.
**--------------------------------------------------------------
*/
{
    if (match(tok,w_RULE)) MaxRules++;
}

/*----------------------------------------------------------------
** Adds rule with format:

 IF SYSTEM CLOCKTIME >= start_time (in seconds)
 AND SYSTEM CLOCKTIME <  stop_time  (in seconds)
 AND TANK id(tank_index) LEVEL [BELOW|ABOVE] level
 THEN PUMP id(pump_index) STATUS IS status

**----------------------------------------------------------------
*/
int addlevelbasedtrigger(int pump_index, int tank_index,
                         int start_time, int stop_time,
                         float level, int status)
{
    struct Premise * p;
    struct Action *a;

    if (status < 0) return(202);  /* Illegal status */
    /* Illegal time interval */
    if (start_time >= stop_time) return(202);

    Rule = (struct aRule *) realloc(Rule, (MaxRules+2)*sizeof(struct aRule));
    if (Rule == NULL) return(101);

    Nrules++;
    snprintf(Rule[Nrules].label,  MAXID, "leveltrigger%d", Nrules);
    Rule[Nrules].Pchain = NULL;
    Rule[Nrules].Tchain = NULL;
    Rule[Nrules].Fchain = NULL;
    Rule[Nrules].priority = 0.0;
    Plast = NULL;

    /* Create new premise structure
       AND SYSTEM CLOCKTIME >= start_time */
    p = (struct Premise *) malloc(sizeof(struct Premise));
    if (p == NULL) return(101);

    p->object = r_SYSTEM;
    p->index =  0;
    p->variable = r_CLOCKTIME;
    p->relop = GE;
    p->logop = r_AND;
    p->status   = 0;
    p->value    = start_time;

    /* Add premise to current rule's premise list */
    p->next = NULL;
    Rule[Nrules].Pchain = p;
    Plast = p;

    /* Create new premise structure
       SYSTEM CLOCKTIME < stop_time */
    p = (struct Premise *) malloc(sizeof(struct Premise));
    if (p == NULL) return(101);

    p->object = r_SYSTEM;
    p->index =  0;
    p->variable = r_CLOCKTIME;
    p->relop = LT;
    p->logop = r_AND;
    p->status   = 0;
    p->value    = stop_time;

    /* Add premise to current rule's premise list */
    p->next = NULL;
    Plast->next = p;
    Plast = p;

    /* Create new premise structure
       AND TANK id(tank_index) LEVEL [BELOW|ABOVE] level */
    p = (struct Premise *) malloc(sizeof(struct Premise));
    if (p == NULL) return(101);

    p->object   = r_NODE;
    p->index    = tank_index;
    p->variable = r_LEVEL;
    p->logop    = r_AND;
    p->status   = 0;
    p->value    = level;
    /* if CLOSED then ABOVE else BELOW */
    p->relop    = ( (status == 0) ? (GT) : (LT));

    /* Add premise to current rule's premise list */
    p->next = NULL;
    Plast->next = p;
    Plast = p;

    /* Create a new action structure:
       THEN PUMP id(pump_index) STATUS IS status */
    a = (struct Action *) malloc(sizeof(struct Action));
    if (a == NULL) return(101);
    a->link = pump_index;
    a->status = -1;
    a->setting = status;

    /* Add action to current rule's action list */
    a->next = Rule[Nrules].Tchain;
    Rule[Nrules].Tchain = a;

    MaxRules++;

    return(0);

}

int  allocrules()
/*
**--------------------------------------------------------------
**    Allocates memory for rule-based controls.
**    Called by allocdata() in EPANET.C module.
**--------------------------------------------------------------
*/
{
    Rule = (struct aRule *) calloc(MaxRules+1,sizeof(struct aRule));
    if (Rule == NULL) return(101);
    else return(0);
}


void freerules()
/*
**--------------------------------------------------------------
**    Frees memory used for rule-based controls.
**    Called by freedata() in EPANET.C module.
**--------------------------------------------------------------
*/
{
    clearrules();
    free(Rule);
}


int checkrules(ENsimulation_t sim, long dt)
/*
**-----------------------------------------------------
**    Checks which rules should fire at current time.
**    Called by ruletimestep() in HYDRAUL.C.
**-----------------------------------------------------
*/
{
    struct  ActItem *actlist = NULL;   /* Linked list of action items */

    int i,
        r;    /* Number of actions actually taken */

    /* Iterate through each rule */
    r = 0;
    for (i=1; i<=Nrules; i++)
    {
        /* If premises true, add THEN clauses to action list. */
        if (evalpremises(sim, i, dt) == TRUE) 
            updateactlist (&actlist, i, Rule[i].Tchain);

        /* If premises false, add ELSE actions to list. */
        else if (Rule[i].Fchain != NULL) 
                updateactlist (&actlist, i, Rule[i].Fchain);
    }

    /* Execute actions then clear list. */
    if (actlist != NULL) r = takeactions (sim, actlist);
    clearactlist (actlist);
    return(r);
}


int  ruledata()
/*
**--------------------------------------------------------------
**    Parses a line from [RULES] section of input.
**    Called by newline() in INPUT2.C module.
**    Tok[] is global array of tokens parsed from input line.
**--------------------------------------------------------------
*/
{
    int    key,                      /* Keyword code */
        err;

    /* Exit if current rule has an error */
    if (RuleState == r_ERROR) return(0);

    /* Find the key word that begins the rule statement */
    err = 0;
    key = findmatch(Tok[0],Ruleword);
    switch (key)
    {
    case -1:     err = 201;      /* Unrecognized keyword */
        break;
    case r_RULE: Nrules++;
        newrule();
        RuleState = r_RULE;
        break;
    case r_IF:   if (RuleState != r_RULE)
    {
        err = 221;   /* Mis-placed IF clause */
        break;
    }
        RuleState = r_IF;
        err = newpremise(r_AND);
        break;
    case r_AND:  if (RuleState == r_IF) err = newpremise(r_AND);
    else if (RuleState == r_THEN || RuleState == r_ELSE)
        err = newaction();
    else err = 221;
        break;
    case r_OR:   if (RuleState == r_IF) err = newpremise(r_OR);
    else err = 221;
        break;
    case r_THEN: if (RuleState != r_IF)
    {
        err = 221;   /* Mis-placed THEN clause */
        break;
    }
        RuleState = r_THEN;
        err = newaction();
        break;
    case r_ELSE: if (RuleState != r_THEN)
    {
        err = 221;   /* Mis-placed ELSE clause */
        break;
    }
        RuleState = r_ELSE;
        err = newaction();
        break;
    case r_PRIORITY: if (RuleState != r_THEN && RuleState != r_ELSE)
    {
        err = 221;
        break;
    }
        RuleState = r_PRIORITY;
        err = newpriority();
        break;
    default:         err = 201;
    }

    /* Set RuleState to r_ERROR if errors found */
    if (err)
    {
        RuleState = r_ERROR;
#if 0
        ruleerrmsg(err);
#endif
        err = 200;
    }
    return(err);
}


static void  clearactlist(struct  ActItem *a)
/*
**----------------------------------------------------------
**    Clears memory used for action list
**----------------------------------------------------------
*/
{
    struct ActItem *anext;
    while (a != NULL)
    {
        anext = a->next;
        free(a);
        a = anext;
    }
}


void  clearrules()
/*
**-----------------------------------------------------------
**    Clears memory used for premises & actions for all rules
**-----------------------------------------------------------
*/
{
    struct Premise *p;
    struct Premise *pnext;
    struct Action  *a;
    struct Action  *anext;
    int i;
    for (i=1; i<=Nrules; i++)
    {
        p = Rule[i].Pchain;
        while (p != NULL)
        {
            pnext = p->next;
            free(p);
            p = pnext;
        }
        a = Rule[i].Tchain;
        while (a != NULL)
        {
            anext = a->next;
            free(a);
            a = anext;
        }
        a = Rule[i].Fchain;
        while (a != NULL)
        {
            anext = a->next;
            free(a);
            a = anext;
        }
    }
}


static void  newrule(void)
/*
**----------------------------------------------------------
**    Adds new rule to rule base
**----------------------------------------------------------
*/
{
    strncpy(Rule[Nrules].label, Tok[1], MAXID);
    Rule[Nrules].Pchain = NULL;
    Rule[Nrules].Tchain = NULL;
    Rule[Nrules].Fchain = NULL;
    Rule[Nrules].priority = 0.0;
    Plast = NULL;
}


static int  newpremise(int logop)
/*
**--------------------------------------------------------------------
**   Adds new premise to current rule.
**   Formats are:
**     IF/AND/OR <object> <id> <variable> <operator> <value>
**     IF/AND/OR  SYSTEM <variable> <operator> <value> (units)
**
**   Calls findmatch() and hour() in INPUT2.C.
**   Calls findnode() and findlink() in EPANET.C.
**---------------------------------------------------------------------
*/
{
    int   i,j,k,m,r,s,v;
    float x;
    struct Premise *p;

    /* Check for correct number of tokens */
    if (Ntokens != 5 && Ntokens != 6) return(201);

    /* Find network object & id if present */
    i = findmatch(Tok[1],Object);
    if (i == r_SYSTEM)
    {
        j = 0;
        v = findmatch(Tok[2],Varword);
        if (v != r_DEMAND && v != r_TIME && v != r_CLOCKTIME) return(201);
    }
    else
    {
        v = findmatch(Tok[3],Varword);
        if (v < 0) return(201);
        switch (i)
        {
        case r_NODE:
        case r_JUNC:
        case r_RESERV:
        case r_TANK:   k = r_NODE; break;
        case r_LINK:
        case r_PIPE:
        case r_PUMP:
        case r_VALVE:  k = r_LINK; break;
        default: return(201);
        }
        i = k;
        if (i == r_NODE)
        {
            j = findnode(Tok[2]);
            if (j == 0) return(203);
            switch (v)
            {
            case r_DEMAND:
            case r_HEAD:
            case r_GRADE:
            case r_LEVEL:
            case r_PRESSURE: break;

/*** Updated 9/7/00 ***/
            case r_FILLTIME:
            case r_DRAINTIME: if (j <= Njuncs) return(201); break;

            default: return(201);
            }
        }
        else
        {
            j = findlink(Tok[2]);
            if (j == 0) return(204);
            switch (v)
            {
            case r_FLOW:
            case r_STATUS:
            case r_SETTING: break;
            default: return(201);
            }
        }
    }

    /* Parse relational operator (r) and check for synonyms */
    if (i == r_SYSTEM) m = 3;
    else m = 4;
    k = findmatch(Tok[m],Operator);
    if (k < 0) return(201);
    switch(k)
    {
    case IS:    r = EQ; break;
    case NOT:   r = NE; break;
    case BELOW: r = LT; break;
    case ABOVE: r = GT; break;
    default:    r = k;
    }

    /* Parse for status (s) or numerical value (x) */
    s = 0;
    x = MISSING;
    if (v == r_TIME || v == r_CLOCKTIME)
    {
        if (Ntokens == 6)
            x = hour(Tok[4],Tok[5])*3600.;
        else
            x = hour(Tok[4],"")*3600.;
        if (x < 0.0) return(202);
    }
    else if ((k = findmatch(Tok[Ntokens-1],Value)) > IS_NUMBER) s = k;
    else
    {
        if (!getfloat(Tok[Ntokens-1],&x)) return(202);
    }

    /* Create new premise structure */
    p = (struct Premise *) malloc(sizeof(struct Premise));
    if (p == NULL) return(101);
    p->object = i;
    p->index =  j;
    p->variable = v;
    p->relop = r;
    p->logop = logop;
    p->status   = s;
    p->value    = x;

    /* Add premise to current rule's premise list */
    p->next = NULL;
    if (Plast == NULL) Rule[Nrules].Pchain = p;
    else Plast->next = p;
    Plast = p;
    return(0);
}


static int  newaction()
/*
**----------------------------------------------------------
**   Adds new action to current rule.
**   Format is:
**      THEN/ELSE/AND LINK <id> <variable> IS <value>
**
**   Calls findlink() from EPANET.C.
**   Calls getfloat() and findmatch() from INPUT2.C.
**----------------------------------------------------------
*/
{
    int   j,k,s;
    float x;
    struct Action *a;

    /* Check for correct number of tokens */
    if (Ntokens != 6) return(201);

    /* Check that link exists */
    j = findlink(Tok[2]);
    if (j == 0) return(204);

/***  Updated 9/7/00  ***/
    /* Cannot control a CV */
    if (Link[j].Type == CV) return(207);

    /* Find value for status or setting */
    s = -1;
    x = MISSING;
    if ((k = findmatch(Tok[5],Value)) > IS_NUMBER) s = k;
    else
    {
        if (!getfloat(Tok[5],&x)) return(202);
        if (x < 0.0) return(202);
    }

/*** Updated 9/7/00 ***/
    /* Cannot change setting for a GPV ***/
    if (x != MISSING && Link[j].Type == GPV) return(202);

/*** Updated 3/1/01 ***/
    /* Set status for pipe in case setting was specified */
    if (x != MISSING && Link[j].Type == PIPE)
    {
        if (x == 0.0) s = IS_CLOSED;
        else          s = IS_OPEN;
        x = MISSING;
    }

    /* Create a new action structure */
    a = (struct Action *) malloc(sizeof(struct Action));
    if (a == NULL) return(101);
    a->link = j;
    a->status = s;
    a->setting = x;

    /* Add action to current rule's action list */
    if (RuleState == r_THEN)
    {
        a->next = Rule[Nrules].Tchain;
        Rule[Nrules].Tchain = a;
    }
    else
    {
        a->next = Rule[Nrules].Fchain;
        Rule[Nrules].Fchain = a;
    }
    return(0);
}


static int  newpriority()
/*
**---------------------------------------------------
**    Adds priority rating to current rule
**---------------------------------------------------
*/
{
    float x;
    if (!getfloat(Tok[1],&x)) return(202);
    Rule[Nrules].priority = x;
    return(0);
}


static int  evalpremises(ENsimulation_t sim, int i, long dt)
/*
**----------------------------------------------------------
**    Checks if premises to rule i are true
**----------------------------------------------------------
*/
{
    int result;
    struct Premise *p;

    result = TRUE;
    p = Rule[i].Pchain;
    while (p != NULL)
    {
        if (p->logop == r_OR)
        {
            if (result == FALSE)
            {
                result = checkpremise (sim, p, dt);
            }
        }
        else
        {
            if (result == FALSE) return(FALSE);
            result = checkpremise (sim, p, dt);
        }
        p = p->next;
    }
    return(result);
}


static int  checkpremise(ENsimulation_t sim, struct Premise *p, long dt)
/*
**----------------------------------------------------------
**    Checks if a particular premise is true
**----------------------------------------------------------
*/
{
    if (p->variable == r_TIME || p->variable == r_CLOCKTIME)
        return(checktime (p, sim->Htime, dt));
    else if (p->status > IS_NUMBER)
        return checkstatus (p, sim->S);
    else
        return(checkvalue (sim, p));
}


static int  checktime(struct Premise *p, long Htime, long dt)
/*
**------------------------------------------------------------
**    Checks if condition on system time holds
**------------------------------------------------------------
*/
{
    char  flag;
    long  t1,t2,x;

    /* Get start and end of rule evaluation time interval */
    if (p->variable == r_TIME)
    {
        /* Start of rule evaluation time interval */
        t1 = Htime - dt + 1;
        t2 = Htime;
    }
    else if (p->variable == r_CLOCKTIME)
    {
        t1 = ((Htime - dt + 1) + Tstart) % SECperDAY;
        t2 = (Htime + Tstart) % SECperDAY;
    }
    else return(0);

    /* Test premise's time */
    x = (long)(p->value);
    switch (p->relop)
    {
        /* For inequality, test against current time */
    case LT: if (t2 >= x) return(0); break;
    case LE: if (t2 >  x) return(0); break;
    case GT: if (t2 <= x) return(0); break;
    case GE: if (t2 <  x) return(0); break;

        /* For equality, test if within interval */
    case EQ:
    case NE:
        flag = FALSE;
        if (t2 < t1)     /* E.g., 11:00 am to 1:00 am */
        {
            if (x >= t1 || x <= t2) flag = TRUE;
        }
        else
        {
            if (x >= t1 && x <= t2) flag = TRUE;
        }
        if (p->relop == EQ && flag == FALSE) return(0);
        if (p->relop == NE && flag == TRUE)  return(0);
        break;
    }

    /* If we get to here then premise was satisfied */
    return(1);
}


static int  checkstatus(struct Premise *p, char *S)
/*
**------------------------------------------------------------
**    Checks if condition on link status holds
**------------------------------------------------------------
*/
{
    char i;
    int  j;
    switch (p->status)
    {
    case IS_OPEN:
    case IS_CLOSED:
    case IS_ACTIVE:
        i = S[p->index];
        if      (i <= CLOSED) j = IS_CLOSED;
        else if (i == ACTIVE) j = IS_ACTIVE;
        else                  j = IS_OPEN;
        if (j == p->status &&
            p->relop == EQ)     return(1);
        if (j != p->status &&
            p->relop == NE) return(1);
    }
    return(0);
}


static int  checkvalue(ENsimulation_t sim, struct Premise *p)
/*
**----------------------------------------------------------
**    Checks if numerical condition on a variable is true.
**    Uses tolerance of 0.001 when testing conditions.
**----------------------------------------------------------
*/
{
    int   i,j,v;
    float x,
        tol = 1.e-3;

    i = p->index;
    v = p->variable;
    switch (v)
    {

/*** Updated 10/25/00 ***/
    case r_DEMAND:    
        if (p->object == r_SYSTEM) x = sim->Dsystem * Ucf[DEMAND];
        else x = sim->D[i] * Ucf[DEMAND];
        break;

    case r_HEAD:
    case r_GRADE:     x = sim->H[i] * Ucf[HEAD];
        break;
    case r_PRESSURE:  x = (sim->H[i] - Node[i].El) * Ucf[PRESSURE];
        break;
    case r_LEVEL:     x = (sim->H[i] - Node[i].El) * Ucf[HEAD];
        break;
    case r_FLOW:      x = ABS (sim->Q[i]) * Ucf[FLOW];
        break;
    case r_SETTING:   
        if (sim->K[i] == MISSING) return(0);
        x = sim->K[i];
        switch (Link[i].Type)
        {
        case PRV:
        case PSV:
        case PBV:  x = x*Ucf[PRESSURE]; break;
        case FCV:  x = x*Ucf[FLOW];     break;
        }
        break;
    case r_FILLTIME:  if (i <= Njuncs) return(0);
        j = i - Njuncs;
        if (Tank[j].A == 0.0) return(0);
        if (sim->D[i] <= TINY) return(0);
        x = (Tank[j].Vmax - sim->V[j]) / sim->D[i];
        break;
    case r_DRAINTIME: if (i <= Njuncs) return(0);
        j = i - Njuncs;
        if (Tank[j].A == 0.0) return(0);
        if (sim->D[i] >= -TINY) return(0);
        x = (Tank[j].Vmin - sim->V[j]) / sim->D[i];
        break;
    default:          return(0);
    }
    switch (p->relop)
    {
    case EQ:        if (ABS(x - p->value) > tol) return(0);
        break;
    case NE:        if (ABS(x - p->value) < tol) return(0);
        break;
    case LT:        if (x > p->value + tol) return(0); break;
    case LE:        if (x > p->value - tol) return(0); break;
    case GT:        if (x < p->value - tol) return(0); break;
    case GE:        if (x < p->value + tol) return(0); break;
    }
    return(1);
}


static void
updateactlist(struct ActItem **actlist, int i, struct Action *actions)
/*
**---------------------------------------------------
**    Adds rule's actions to action list
**---------------------------------------------------
*/
{
    struct ActItem *item;
    struct Action *a;

    /* Iterate through each action of Rule i */
    a = actions;
    while (a != NULL)
    {
        /* Add action to list if not already on it */
        if (!checkaction (*actlist, i, a))
        {
            item = (struct ActItem *) malloc(sizeof(struct ActItem));
            if (item != NULL)
            {
                item->action = a;
                item->ruleindex = i;
                item->next = *actlist;
                *actlist = item;
            }
        }
        a = a->next;
    }
}


static int  checkaction(struct ActItem *item, int i, struct Action *a)
/*
**-----------------------------------------------------------
**    Checks if an action is already on the Action List
**-----------------------------------------------------------
*/
{
    int i1,k,k1;
    struct Action *a1;

    /* Search action list for link named in action */
    k = a->link;                 /* Action applies to link k */
    while (item != NULL)
    {
        a1 = item->action;
        i1 = item->ruleindex;
        k1 = a1->link;

        /* If link on list then replace action if rule has higher priority. */
        if (k1 == k)
        {
            if (Rule[i].priority > Rule[i1].priority)
            {
                item->action = a;
                item->ruleindex = i;
            }
            return(1);
        }
        item = item->next;
    }
    return(0);
}


static int  takeactions(ENsimulation_t sim, struct ActItem *item)
/*
**-----------------------------------------------------------
**    Implements actions on action list
**-----------------------------------------------------------
*/
{
    struct Action *a;

    char   flag;
    int    k, s, n;
    float  tol = 1.e-3,
        v, x;

    n = 0;
    while (item != NULL)
    {
        flag = FALSE;
        a = item->action;
        k = a->link;
        s = sim->S[k];
        v = sim->K[k];
        x = a->setting;

        /* Switch link from closed to open */
        if (a->status == IS_OPEN && s <= CLOSED)
        {
            setlinkstatus (sim, k, 1, &(sim->S[k]), &(sim->K[k]));
            flag = TRUE;
        }

        /* Switch link from not closed to closed */
        else if (a->status == IS_CLOSED && s > CLOSED)
        {
            setlinkstatus (sim, k, 0, &(sim->S[k]), &(sim->K[k]));
            flag = TRUE;
        }

        /* Change link's setting */
        else if (x != MISSING)
        {
            switch(Link[k].Type)
            {
            case PRV:
            case PSV:
            case PBV:    x = x/Ucf[PRESSURE];  break;
            case FCV:    x = x/Ucf[FLOW];      break;
            }
            if (ABS(x-v) > tol)
            {
                setlinksetting (sim, k, x, &(sim->S[k]), &(sim->K[k]));
                flag = TRUE;
            }
        }

        /* Report rule action */
        if (flag == TRUE)
        {
            n++;
#if 0
            if (Statflag) writeruleaction(k,Rule[item->ruleindex].label);
#endif
        }

        /* Move to next action on list */
        item = item->next;
    }
    return(n);
}

#if 0
static void  ruleerrmsg(int err)
/*
**-----------------------------------------------------------
**    Reports error message
**-----------------------------------------------------------
*/
{
    int    i;
    char   label[81];
    char   fmt[256];
    switch (err)
    {
    case 201:   strcpy(fmt,R_ERR201);  break;
    case 202:   strcpy(fmt,R_ERR202);  break;
    case 203:   strcpy(fmt,R_ERR203);  break;
    case 204:   strcpy(fmt,R_ERR204);  break;

/***  Updated on 9/7/00  ***/
    case 207:   strcpy(fmt,R_ERR207);  break;

    case 221:   strcpy(fmt,R_ERR221);  break;
    default:    return;
    }
    if (Nrules > 0)
    {
        strcpy(label,t_RULE);
        strcat(label," ");
        strcat(label,Rule[Nrules].label);
    }
    else strcpy(label,t_RULES_SECT);
    sprintf(Msg,fmt);
    strcat(Msg,label);
    strcat(Msg,":");
    writeline(Msg);
    strcpy(fmt,Tok[0]);
    for (i=1; i<Ntokens; i++)
    {
        strcat(fmt," ");
        strcat(fmt,Tok[i]);
    }
    writeline(fmt);
}
#endif
/***************** END OF RULES.C ******************/
