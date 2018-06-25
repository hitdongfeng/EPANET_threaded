struct   Sadjlist         /* NODE ADJACENCY LIST ITEM */
{
    int   node;            /* Index of connecting node */
    int   link;            /* Index of connecting link */
    struct Sadjlist *next; /* Next item in list        */
};
/* Pointer to adjacency list item */
typedef struct Sadjlist *Padjlist;

struct _smatrix_t {
/*
** The following arrays store the positions of the non-zero coeffs.
** of the lower triangular portion of A whose values are stored in Aij:
*/
    int      *XLNZ,       /* Start position of each column in NZSUB  */
             *NZSUB,      /* Row index of each coeff. in each column */
             *LNZ;        /* Position of each coeff. in Aij array    */

    int      Ncoeffs,     /* Number of non-0 matrix coeffs       */
             *Order,      /* Row-to-node of A                    */
             *Row,        /* Node-to-row of A                    */
             *Ndx;        /* Index of link's coeff. in Aij       */

    Padjlist *adjlist;              /* Node adjacency lists         */

};

typedef struct _smatrix_t *smatrix_t;

smatrix_t createsparse (void);               /* Creates sparse matrix      */
void  freesparse (smatrix_t smatrix);        /* Frees matrix memory        */
int smatrix_non_zero_coeffs (smatrix_t smatrix); /* Substitutes Ncoeffs */
Padjlist smatrix_adjlist_item (smatrix_t smatrix, int i);
int     linsolve (smatrix_t smatrix, int,    /* Solution of linear eqns.   */ 
                  REAL *, REAL *, REAL *);   /* via Cholesky factorization */

