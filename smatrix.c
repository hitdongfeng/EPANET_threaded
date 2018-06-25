/*
*******************************************************************

SMATRIX.C -- Sparse matrix routines for EPANET program.

AUTHOR:     M. Lopez-Ibanez (M. L\'opez-Ib\'a\~nez)
            Napier University, Edinburgh, UK.

$LastChangedBy: manu $ $Revision: 145 $
$LastChangedDate: 2007-11-19 21:03:44 +0100 (Mon, 19 Nov 2007) $

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
AUTHOR:     L. Rossman
            US EPA - NRMRL

This module contains the sparse matrix routines used to solve
a network's hydraulic equations. The entry points into this
module are:
   createsparse() -- called from openhyd() in HYDRAUL.C
   freesparse()   -- called from closehyd() in HYDRAUL.C
   linsolve()     -- called from netsolve() in HYDRAUL.C

Createsparse() does the following:
   1. for each node, builds an adjacency list that identifies
      all links connected to the node (see buildlists())
   2. re-orders the network's nodes to minimize the number
      of non-zero entries in the hydraulic solution matrix
      (see reorder())
   3. converts the adjacency lists into a compact scheme
      for storing the non-zero coeffs. in the lower diagonal
      portion of the solution matrix (see storesparse())
Freesparse() frees the memory used for the sparse matrix.
Linsolve() solves the linearized system of hydraulic equations.

********************************************************************
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
/* #include "smatrix.h" We need to include smatrix.h in types.h for now, so this is not needed */

static smatrix_t allocsparse(void);        /* Allocates matrix memory    */

/* Builds adjacency lists     */
static int  buildlists(Padjlist *adjlist, int *ndx, int paraflag);      

/* Checks for parallel links  */
static int  paralink(Padjlist *adjlist, int *ndx, int i, int j, int k); 

static void    xparalinks(Padjlist *adjlist);           /* Removes parallel links     */
static void  freelists(Padjlist *adjlist); /* Frees adjacency lists      */
static void  countdegree(Padjlist *adjlist, int * Degree);    /* Counts links at each node  */
static int   reordernodes(smatrix_t smatrix, int * Degree);   /* Finds a node re-ordering   */

/* Finds min. degree node     */
static int  mindegree(int *Order, int * Degree, int k, int n);

static int  growlist(smatrix_t smatrix, int *Degree, int knode);     /* Augments adjacency list    */
static int  newlink(smatrix_t smatrix, int * Degree, Padjlist alink); /* Adds fill-ins for a node   */
static int     linked(Padjlist *adjlist, int, int);         /* Checks if 2 nodes linked   */
static int     addlink(Padjlist *adjlist, int, int, int);   /* Creates new fill-in        */
static int  storesparse(smatrix_t smatrix, int n); /* Stores sparse matrix       */
static int  ordersparse (smatrix_t smatrix, int n); /* Orders matrix storage      */
static void    transpose(int,int *,int *,/* Transposes sparse matrix   */
                         int *,int *,int *,int *,int *);

smatrix_t createsparse(void)
/*
**--------------------------------------------------------------
** Input:   none
** Output:  returns error code
** Purpose: creates sparse representation of coeff. matrix
**--------------------------------------------------------------
*/
{
    smatrix_t smatrix;

    int      *Degree;     /* Number of links adjacent to each node  */

    int errcode = 0;

    /* Build node-link adjacency lists with parallel links removed. */
    Degree = calloc (Nnodes+1, sizeof(int));

    /* Allocate data structures */
    smatrix = allocsparse ();
    ERRCODE(buildlists (smatrix->adjlist, smatrix->Ndx, TRUE));
    if (!errcode)
    {
        xparalinks (smatrix->adjlist);    /* Remove parallel links */
        countdegree (smatrix->adjlist, Degree);   /* Find degree of each junction */
    }                   /* (= # of adjacent links)  */

    /* Re-order nodes to minimize number of non-zero coeffs.    */
    /* in factorized solution matrix. At same time, adjacency   */
    /* list is updated with links representing non-zero coeffs. */
    smatrix->Ncoeffs = Nlinks;
    ERRCODE(reordernodes (smatrix, Degree));

    /* Allocate sparse matrix storage */
    smatrix->XLNZ  = calloc(Njuncs + 2, sizeof(int));
    smatrix->NZSUB = calloc(smatrix->Ncoeffs + 2, sizeof(int));
    smatrix->LNZ   = calloc(smatrix->Ncoeffs + 2, sizeof(int));

    /* Allocate memory for sparse storage of positions of non-zero */
    /* coeffs. and store these positions in vector NZSUB. */
    ERRCODE(storesparse (smatrix, Njuncs));

    /* Free memory used for adjacency lists and sort */
    /* row indexes in NZSUB to optimize linsolve().  */
    if (!errcode) freelists (smatrix->adjlist);
    ERRCODE(ordersparse (smatrix, Njuncs));

    /* Re-build adjacency lists without removing parallel */
    /* links for use in future connectivity checking.     */
    ERRCODE(buildlists (smatrix->adjlist, smatrix->Ndx, FALSE));

    /* Free allocated memory */
    free(Degree);
    return smatrix;
}                        /* End of createsparse */


static smatrix_t allocsparse(void)
/*
**--------------------------------------------------------------
** Input:   none
** Output:  returns error code
** Purpose: allocates memory for indexing the solution matrix
**--------------------------------------------------------------
*/
{
    smatrix_t smatrix;
    
    smatrix = malloc (sizeof(struct _smatrix_t));

    smatrix->adjlist = calloc (Nnodes+1,  sizeof(Padjlist));
    smatrix->Order  = calloc (Nnodes+1,  sizeof(int));
    smatrix->Row    = calloc (Nnodes+1,  sizeof(int));
    smatrix->Ndx    = calloc (Nlinks+1,  sizeof(int));
/*    ERRCODE(MEMCHECK(Adjlist));
    ERRCODE(MEMCHECK(Order));
    ERRCODE(MEMCHECK(Row));
    ERRCODE(MEMCHECK(smatrix->Ndx)); */
    return smatrix;
}


void  freesparse(smatrix_t smatrix)
/*
**----------------------------------------------------------------
** Input:   None
** Output:  None
** Purpose: Frees memory used for sparse matrix storage
**----------------------------------------------------------------
*/
{
    freelists (smatrix->adjlist);
    free (smatrix->adjlist);
    free (smatrix->Order);
    free (smatrix->Row);
    free (smatrix->Ndx);
    free (smatrix->XLNZ);
    free (smatrix->NZSUB);
    free (smatrix->LNZ);
    free (smatrix);
}                        /* End of freesparse */


static int  buildlists(Padjlist * adjlist, int *ndx, int paraflag)
/*
**--------------------------------------------------------------
** Input:   paraflag = TRUE if list marks parallel links
** Output:  returns error code
** Purpose: builds linked list of links adjacent to each node
**--------------------------------------------------------------
*/
{
    int    i,j,k;
    int    pmark = 0;
    int    errcode = 0;
    Padjlist  alink;

    /* For each link, update adjacency lists of its end nodes */
    for (k=1; k<=Nlinks; k++)
    {
        i = Link[k].N1;
        j = Link[k].N2;
        if (paraflag) 
            pmark = paralink (adjlist, ndx, i, j, k); /* Parallel link check */

        /* Include link in start node i's list */
        alink = (struct Sadjlist *) malloc(sizeof(struct Sadjlist));
        if (alink == NULL) return(101);
        if (!pmark) alink->node = j;
        else        alink->node = 0;           /* Parallel link marker */
        alink->link = k;
        alink->next = adjlist[i];
        adjlist[i] = alink;

        /* Include link in end node j's list */
        alink = (struct Sadjlist *) malloc(sizeof(struct Sadjlist));
        if (alink == NULL) return(101);
        if (!pmark) alink->node = i;
        else        alink->node = 0;           /* Parallel link marker */
        alink->link = k;
        alink->next = adjlist[j];
        adjlist[j] = alink;
    }
    return(errcode);
}                        /* End of buildlists */


static int  paralink(Padjlist * adjlist, int *ndx, int i, int j, int k)
/*
**--------------------------------------------------------------
** Input:   i = index of start node of link
**          j = index of end node of link
**          k = link index
** Output:  returns 1 if link k parallels another link, else 0
** Purpose: checks for parallel links between nodes i and j
**
**--------------------------------------------------------------
*/
{
    Padjlist alink;
    for (alink = adjlist[i]; alink != NULL; alink = alink->next)
    {
        if (alink->node == j)     /* Link || to k (same end nodes) */
        {
            ndx[k] = alink->link;  /* Assign Ndx entry to this link */
            return(1);
        }
    }
    ndx[k] = k;                  /* Ndx entry if link not parallel */
    return(0);
}                        /* End of paralink */


static void  xparalinks(Padjlist *adjlist)
/*
**--------------------------------------------------------------
** Input:   none
** Output:  none
** Purpose: removes parallel links from nodal adjacency lists
**--------------------------------------------------------------
*/
{
    int    i;
    Padjlist  alink,       /* Current item in adjacency list */
        blink;       /* Previous item in adjacency list */

    /* Scan adjacency list of each node */
    for (i=1; i<=Nnodes; i++)
    {
        alink = adjlist[i];              /* First item in list */
        blink = NULL;
        while (alink != NULL)
        {
            if (alink->node == 0)      /* Parallel link marker found */
            {
                if (blink == NULL)      /* This holds at start of list */
                {
                    adjlist[i] = alink->next;
                    free(alink);             /* Remove item from list */
                    alink = adjlist[i];
                }
                else                    /* This holds for interior of list */
                {
                    blink->next = alink->next;
                    free(alink);             /* Remove item from list */
                    alink = blink->next;
                }
            }
            else
            {
                blink = alink;          /* Move to next item in list */
                alink = alink->next;
            }
        }
    }
}                        /* End of xparalinks */


static void  freelists(Padjlist *adjlist)
/*
**--------------------------------------------------------------
** Input:   none
** Output:  none
** Purpose: frees memory used for nodal adjacency lists
**--------------------------------------------------------------
*/
{
    int   i;
    Padjlist alink;

    for (i=0; i<=Nnodes; i++)
    {
        for (alink = adjlist[i]; alink != NULL; alink = adjlist[i])
        {
            adjlist[i] = alink->next;
            free(alink);
        }
    }
}                        /* End of freelists */


static void  countdegree(Padjlist * adjlist, int * Degree)
/*
**----------------------------------------------------------------
** Input:   none
** Output:  none
** Purpose: counts number of nodes directly connected to each node
**----------------------------------------------------------------
*/
{
    int   i;
    Padjlist alink;
    memset(Degree,0,(Nnodes+1)*sizeof(int));

    /* NOTE: For purposes of node re-ordering, Tanks (nodes with  */
    /*       indexes above Njuncs) have zero degree of adjacency. */

    for (i=1; i<=Njuncs; i++)
        for (alink = adjlist[i]; alink != NULL; alink = alink->next)
            if (alink->node > 0) Degree[i]++;
}


static int   reordernodes(smatrix_t smatrix, int * Degree)
/*
**--------------------------------------------------------------
** Input:   none
** Output:  returns 1 if successful, 0 if not
** Purpose: re-orders nodes to minimize # of non-zeros that
**          will appear in factorized solution matrix
**--------------------------------------------------------------
*/
{
    int *Row = smatrix->Row;
    int *Order = smatrix->Order;

    int k, knode, m, n;

    for (k=1; k<=Nnodes; k++)
    {
        Row[k] = k;
        Order[k] = k;
    }
    n = Njuncs;
    for (k=1; k<=n; k++)                   /* Examine each junction    */
    {
        m = mindegree (Order, Degree, k, n);  /* Node with lowest degree  */
        knode = Order[m];                   /* Node's index             */
        if (!growlist (smatrix, Degree, knode)) /* Augment adjacency list   */
            return(101); 
        Order[m] = Order[k];                /* Switch order of nodes    */
        Order[k] = knode;
        Degree[knode] = 0;                  /* In-activate node         */
    }
    for (k=1; k<=n; k++)                   /* Assign nodes to rows of  */
        Row[Order[k]] = k;                 /*   coeff. matrix          */
    return(0);
}                        /* End of reordernodes */


static int  mindegree(int *Order, int * Degree, int k, int n)
/*
**--------------------------------------------------------------
** Input:   k = first node in list of active nodes
**          n = total number of junction nodes
** Output:  returns node index with fewest direct connections
** Purpose: finds active node with fewest direct connections
**--------------------------------------------------------------
*/
{
    int i, m;
    int min = n,
        imin = n;

    for (i=k; i<=n; i++)
    {
        m = Degree[Order[i]];
        if (m < min)
        {
            min = m;
            imin = i;
        }
    }
    return(imin);
}                        /* End of mindegree */


static int  growlist(smatrix_t smatrix, int *Degree, int knode)
/*
**--------------------------------------------------------------
** Input:   knode = node index
** Output:  returns 1 if successful, 0 if not
** Purpose: creates new entries in knode's adjacency list for
**          all unlinked pairs of active nodes that are
**          adjacent to knode
**--------------------------------------------------------------
*/
{
    int   node;
    Padjlist alink;

    /* Iterate through all nodes connected to knode */
    for (alink = smatrix->adjlist[knode]; alink != NULL; alink = alink -> next)
    {
        node = alink->node;       /* End node of connecting link  */
        if (Degree[node] > 0)     /* End node is active           */
        {
            Degree[node]--;        /* Reduce degree of adjacency   */
            if (!newlink (smatrix, Degree, alink))   /* Add to adjacency list        */
                return(0);
        }
    }
    return(1);
}                        /* End of growlist */


static int  newlink (smatrix_t smatrix, int * Degree, Padjlist alink)
/*
**--------------------------------------------------------------
** Input:   alink = element of node's adjacency list
** Output:  returns 1 if successful, 0 if not
** Purpose: links end of current adjacent link to end nodes of
**          all links that follow it on adjacency list
**--------------------------------------------------------------
*/
{
    int   inode, jnode;
    Padjlist blink;

    /* Scan all entries in adjacency list that follow anode. */
    inode = alink->node;             /* End node of connection to anode */
    for (blink = alink->next; blink != NULL; blink = blink->next)
    {
        jnode = blink->node;          /* End node of next connection */

        /* If jnode still active, and inode not connected to jnode, */
        /* then add a new connection between inode and jnode.       */
        if (Degree[jnode] > 0)        /* jnode still active */
        {
            if (!linked (smatrix->adjlist, inode,jnode))  /* inode not linked to jnode */
            {

                /* Since new connection represents a non-zero coeff. */
                /* in the solution matrix, update the coeff. count.  */
                smatrix->Ncoeffs++;

                /* Update adjacency lists for inode & jnode to */
                /* reflect the new connection.                 */
                if (!addlink (smatrix->adjlist, inode,jnode, smatrix->Ncoeffs)) return 0;
                if (!addlink (smatrix->adjlist, jnode,inode, smatrix->Ncoeffs)) return 0;
                Degree[inode]++;
                Degree[jnode]++;
            }
        }
    }
    return(1);
}                        /* End of newlink */


static int  linked(Padjlist *adjlist, int i, int j)
/*
**--------------------------------------------------------------
** Input:   i = node index
**          j = node index
** Output:  returns 1 if nodes i and j are linked, 0 if not
** Purpose: checks if nodes i and j are already linked.
**--------------------------------------------------------------
*/
{
    Padjlist alink;
    for (alink = adjlist[i]; alink != NULL; alink = alink->next)
        if (alink->node == j) return(1);
    return(0);
}                        /* End of linked */


static int  addlink(Padjlist *adjlist, int i, int j, int n)
/*
**--------------------------------------------------------------
** Input:   i = node index
**          j = node index
**          n = link index
** Output:  returns 1 if successful, 0 if not
** Purpose: augments node i's adjacency list with node j
**--------------------------------------------------------------
*/
{
    Padjlist alink;
    alink = malloc(sizeof(struct Sadjlist));
    if (alink == NULL) return(0);
    alink->node = j;
    alink->link = n;
    alink->next = adjlist[i];
    adjlist[i] = alink;
    return(1);
}                        /* End of addlink */


static int  storesparse(smatrix_t smatrix, int n)
/*
**--------------------------------------------------------------
** Input:   n = number of rows in solution matrix
** Output:  returns error code
** Purpose: stores row indexes of non-zeros of each column of
**          lower triangular portion of factorized matrix
**--------------------------------------------------------------
*/
{
    Padjlist *adjlist = smatrix->adjlist;
    int  *XLNZ = smatrix->XLNZ;   /* Start position of each column in NZSUB  */
    int  *NZSUB = smatrix->NZSUB; /* Row index of each coeff. in each column */
    int  *LNZ = smatrix->LNZ;     /* Position of each coeff. in Aij array    */
    int  *Row = smatrix->Row;
    int  *Order = smatrix->Order;

    Padjlist alink;
    int   i, ii, j, k, l, m;
    int   errcode = 0;

    /* Generate row index pointers for each column of matrix */
    k = 0;
    XLNZ[1] = 1;
    for (i=1; i<=n; i++)             /* column */
    {
        m = 0;
        ii = Order[i];
        for (alink = adjlist[ii]; alink != NULL; alink = alink->next)
        {
            j = Row[alink->node];    /* row */
            l = alink->link;
            if (j > i && j <= n)
            {
                m++;
                k++;
                NZSUB[k] = j;
                LNZ[k] = l;
            }
        }
        XLNZ[i+1] = XLNZ[i] + m;
    }
    return(errcode);
}                        /* End of storesparse */


static int  ordersparse(smatrix_t smatrix, int n)
/*
**--------------------------------------------------------------
** Input:   n = number of rows in solution matrix
** Output:  returns eror code
** Purpose: puts row indexes in ascending order in NZSUB
**--------------------------------------------------------------
*/
{
    int  *XLNZ = smatrix->XLNZ;      /* Start position of each column in NZSUB  */
    int  *NZSUB = smatrix->NZSUB;    /* Row index of each coeff. in each column */
    int  *LNZ = smatrix->LNZ;        /* Position of each coeff. in Aij array    */

    int  i, k;
    int  *xlnzt, *nzsubt, *lnzt, *nzt;
    int  errcode = 0;

    xlnzt  = calloc(n+2, sizeof(int));
    nzsubt = calloc (smatrix->Ncoeffs + 2, sizeof(int));
    lnzt   = calloc (smatrix->Ncoeffs + 2, sizeof(int));
    nzt    = calloc(n+2, sizeof(int));
    ERRCODE(MEMCHECK(xlnzt));
    ERRCODE(MEMCHECK(nzsubt));
    ERRCODE(MEMCHECK(lnzt));
    ERRCODE(MEMCHECK(nzt));
    if (!errcode)
    {

        /* Count # non-zeros in each row */
        for (i=1; i<=n; i++) nzt[i] = 0;
        for (i=1; i<=n; i++)
        {
            for (k=XLNZ[i]; k<XLNZ[i+1]; k++) nzt[NZSUB[k]]++;
        }
        xlnzt[1] = 1;
        for (i=1; i<=n; i++) xlnzt[i+1] = xlnzt[i] + nzt[i];

        /* Transpose matrix twice to order column indexes */
        transpose(n,XLNZ,NZSUB,LNZ,xlnzt,nzsubt,lnzt,nzt);
        transpose(n,xlnzt,nzsubt,lnzt,XLNZ,NZSUB,LNZ,nzt);
    }

    /* Reclaim memory */
    free(xlnzt);
    free(nzsubt);
    free(lnzt);
    free(nzt);
    return(errcode);
}                        /* End of ordersparse */


static void  transpose(int n, int *il, int *jl, int *xl, int *ilt, int *jlt,
                int *xlt, int *nzt)
/*
**---------------------------------------------------------------------
** Input:   n = matrix order
**          il,jl,xl = sparse storage scheme for original matrix
**          nzt = work array
** Output:  ilt,jlt,xlt = sparse storage scheme for transposed matrix
** Purpose: Determines sparse storage scheme for transpose of a matrix
**---------------------------------------------------------------------
*/
{
    int  i, j, k, kk;

    for (i=1; i<=n; i++) nzt[i] = 0;
    for (i=1; i<=n; i++)
    {
        for (k=il[i]; k<il[i+1]; k++)
        {
            j = jl[k];
            kk = ilt[j] + nzt[j];
            jlt[kk] = i;
            xlt[kk] = xl[k];
            nzt[j]++;
        }
    }
}                        /* End of transpose */


int  linsolve(smatrix_t smatrix, int n, REAL *Aii, REAL *Aij, REAL *B)
/*
**--------------------------------------------------------------
** Input:   n    = number of equations
**          Aii  = diagonal entries of solution matrix
**          Aij  = non-zero off-diagonal entries of matrix
**          B    = right hand side coeffs.
** Output:  B    = solution values
**          returns 0 if solution found, or node index 
**          causing system to be ill-conditioned
** Purpose: solves sparse symmetric system of linear
**          equations using Cholesky factorization
**
** NOTE:   This procedure assumes that the solution matrix has
**         been symbolically factorized with the positions of
**         the lower triangular, off-diagonal, non-zero coeffs.
**         stored in the following integer arrays:
**            XLNZ  (start position of each column in NZSUB)
**            NZSUB (row index of each non-zero in each column)
**            LNZ   (position of each NZSUB entry in Aij array)
**
**  This procedure has been adapted from subroutines GSFCT and
**  GSSLV in the book "Computer Solution of Large Sparse
**  Positive Definite Systems" by A. George and J. W-H Liu
**  (Prentice-Hall, 1981).
**--------------------------------------------------------------
*/
{
    int  *XLNZ = smatrix->XLNZ;   /* Start position of each column in NZSUB  */
    int  *NZSUB = smatrix->NZSUB; /* Row index of each coeff. in each column */
    int  *LNZ = smatrix->LNZ;     /* Position of each coeff. in Aij array    */
    int  *Order = smatrix->Order;

    int  *link, *first;
    int  i, istop, istrt, j, k, kfirst;
    
    int  errcode = 0;
    REAL bj, ljk;
    REAL *temp;

    temp = (REAL *) calloc(n+1, sizeof(REAL));
    link = (int *) calloc(n+1,sizeof(int));
    first = (int *) calloc(n+1,sizeof(int));
    ERRCODE(MEMCHECK(temp));
    ERRCODE(MEMCHECK(link));
    ERRCODE(MEMCHECK(first));
    if (errcode)
    {
        errcode = -errcode;
        return errcode;
    }

    /* Begin numerical factorization of matrix A into L */
    /*   Compute column L(*,j) for j = 1,...n */
    for (j=1; j<=n; j++)
    {
        /* For each column L(*,k) that affects L(*,j): */
        REAL diagj = 0.0;
        k = link[j];
        while (k != 0)
        {

            /* Outer product modification of L(*,j) by  */
            /* L(*,k) starting at first[k] of L(*,k).   */
            int newk = link[k];
            kfirst = first[k];
            ljk = Aij[LNZ[kfirst]];
            diagj += ljk*ljk;
            istrt = kfirst + 1;
            istop = XLNZ[k+1];
            if (istop > istrt)
            {

                /* Before modification, update vectors 'first' */
                /* and 'link' for future modification steps.   */
                first[k] = istrt;
                link[k] = link[NZSUB[istrt]];
                link[NZSUB[istrt]] = k;

                /* The actual mod is saved in vector 'temp'. */
                for (i=istrt; i< istop; i++)
                {
                    temp[NZSUB[i]] += Aij[LNZ[i]]*ljk;
                }
            }
            k = newk;
        }

        /* Apply the modifications accumulated */
        /* in 'temp' to column L(*,j).         */
        diagj = Aii[j] - diagj;
        /* FIXME: this is broken. The correct code is:
        if (diagj < 0)
        {
            errcode = j;
            free(temp);
            free(link);
            free(first);
            return errcode;
        }
        diagj = sqrt(diagj);
        if (diagj == 0)
        {
            errcode = j;
            free(temp);
            free(link);
            free(first);
            return errcode;
        }
        */
        if (diagj <= 0.0)        /* Check for ill-conditioning */
        {
            free(temp);
            free(link);
            free(first);
            return Order[j];
        }
        diagj = sqrt(diagj);
        Aii[j] = diagj;
        istrt = XLNZ[j];
        istop = XLNZ[j+1];
        if (istop > istrt)
        {
            first[j] = istrt;
            link[j] = link[NZSUB[istrt]];
            link[NZSUB[istrt]] = j;
            for (i=istrt; i< istop; i++)
            {
                Aij[LNZ[i]] = (Aij[LNZ[i]] - temp[NZSUB[i]])/diagj;
                temp[NZSUB[i]] = 0.0;
            }
        }
    }      /* next j */

    /* Foward substitution */
    for (j=1; j<=n; j++)
    {
        bj = B[j]/Aii[j];
        B[j] = bj;
        istop = XLNZ[j+1];
        for (i=XLNZ[j]; i< istop; i++)
        {
            B[NZSUB[i]] -= Aij[LNZ[i]]*bj;
        }
    }

    /* Backward substitution */
    for (j=n; j>=1; j--)
    {
        bj = B[j];
        istop = XLNZ[j+1];
        for (i=XLNZ[j]; i< istop; i++)
        {
            bj -= Aij[LNZ[i]]*B[NZSUB[i]];
        }
        B[j] = bj/Aii[j];
    }

    free(temp);
    free(link);
    free(first);
    
    return 0;
}                        /* End of linsolve */

int smatrix_non_zero_coeffs (smatrix_t smatrix)
{
    return smatrix->Ncoeffs;
}

Padjlist smatrix_adjlist_item (smatrix_t smatrix, int i)
{
    return smatrix->adjlist[i];
}

/************************ END OF SMATRIX.C ************************/
