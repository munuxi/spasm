#include <assert.h>
#include "spasm.h"

/*
 * Permutations matrices are represented by vectors.
 *
 * p[k] = i means that P[k,i] = 1
 */


/*
 * x <-- P.b (or, equivalently, x <-- b.(P^{-1}), for dense vectors x and b;
 * p=NULL denotes identity.
 *
 * This means that x[k] <--- b[ p[k] ]
 */
void spasm_pvec(const int *p, const spasm_GFp * b, spasm_GFp * x, int n) {
    int k;
    assert(x != NULL);
    assert(b != NULL);

    for (k = 0; k < n; k++) {
        x[k] = b[(p != SPASM_IDENTITY_PERMUTATION) ? p[k] : k];
    }
}

/* x <--- (P^{-1}).b (or x <--- b.P), for dense vectors x and b;
 * p=NULL denotes identity.
 *
 * This means that x[ p[k] ] <--- b[ k ]
 *
 * The function is given p, not p^{-1}.
 */
void spasm_ipvec(const int *p, const spasm_GFp * b, spasm_GFp * x, int n) {
    int k;
    assert(x != NULL);
    assert(b != NULL);

    for (k = 0; k < n; k++) {
        x[(p != SPASM_IDENTITY_PERMUTATION) ? p[k] : k] = b[k];
    }
}

/* compute the inverse permutation */
int *spasm_pinv(int const *p, int n) {
    int k, *pinv;
    /* p = NULL denotes identity */
    if (p == NULL) {
      return NULL;
    }
    /* allocate result */
    pinv = spasm_malloc(n * sizeof(int));
    /* invert the permutation */
    for (k = 0; k < n; k++) {
        pinv[ p[k] ] = k;
    }
    return pinv;
}


/*
 * C = P.A.Q where p and q are permutations of 0..m-1 and 0..n-1
 * respectively.
 *
 */
spasm *spasm_permute(const spasm * A, const int *p, const int *q, int values) {
    int t, j, i, nz, m, n, *Ap, *Aj, *Cp, *Cj;
    spasm_GFp *Cx, *Ax;
    spasm *C;

    /* check inputs */
    assert(A != NULL);

    n = A->n;
    m = A->m;
    Ap = A->p;
    Aj = A->j;
    Ax = A->x;

    /* alloc result */
    C = spasm_csr_alloc(n, m, A->nzmax, A->prime, values && (Ax != NULL));
    Cp = C->p;
    Cj = C->j;
    Cx = C->x;
    nz = 0;

    for (i = 0; i < n; i++) {
        /* row i of C is row p[i] of A (denoted by j) */
        Cp[i] = nz;
        j = (p != NULL) ? p[i] : i;
        for(t = Ap[j]; t < Ap[j + 1]; t++) {
            /* col i of A is col pinv[i] of C */
            Cj[nz] = (q != NULL) ? q[Aj[t]] : Aj[t];
            if (Cx != NULL) {
                Cx[nz] = Ax[t];
            }
            nz++;
        }
    }
    /* finalize the last row of C */
    Cp[n] = nz;
    return C;
}
