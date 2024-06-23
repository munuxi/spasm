
#include <assert.h>

#include "flint/perm.h"
#include "flint/nmod.h"
#include "flint/nmod_mat.h"

#include "spasm.h"



// /**
//  * Conversion of a permutation from LAPACK format to Math format
//  */
// inline void LAPACKPerm2MathPerm (slong * MathP, const slong * LapackP,
//                                  const slong N)
// {
//     for (auto i=0; i<N; i++)
//         MathP[i] = i;
//     for (auto i=0; i<N; i++){
//         if (LapackP[i] != i){
//             std::swap(MathP[i],MathP[LapackP[i]]);
//         }
//     }
// }

// /**
//  * Conversion of a permutation from Maths format to LAPACK format
//  */
// inline void MathPerm2LAPACKPerm (slong * LapackP, const slong * MathP,
//                                  const slong N)
// {
//     slong * T = _perm_init(N);
//     slong * Tinv = _perm_init(N);
//     for (slong i=0; i<N; i++){
//         T[i] =i;
//         Tinv[i] = i;
//     }
//     for (slong i=0; i<N; i++){
//         slong j = Tinv [MathP [i]];
//         LapackP [i] = j;
//         slong tmp = T[j];
//         T[j] = T[i];
//         Tinv[T[i]] = j;
//         T[i] = tmp;
//         Tinv[tmp] = i;
//     }
//     _perm_clear(T);
//     _perm_clear(Tinv);
// }

slong
nmod_mat_rref_2(nmod_mat_t A,  slong * pivots_nonpivots)
{
    slong rank, *P;
    if (nmod_mat_is_empty(A))
        return 0;

    if (A->r == 1)
    {
        ulong c, cinv;
        slong i, j;
        slong r = 0;

        for (i = 0; i < A->c; i++)
        {
            c = nmod_mat_entry(A, 0, i);
            if (c != 0)
            {
                r = 1;
                if (c == 1)
                    break;

                cinv = nmod_inv(c, A->mod);
                nmod_mat_set_entry(A, 0, i, 1);
                for (j = i + 1;j < A->c; j++)
                {
                    nmod_mat_set_entry(A, 0, j, nmod_mul(nmod_mat_get_entry(A, 0, j), cinv, A->mod));
                }
                break;
            }
        }
        return r;
    }
    P = _perm_init(A->r);
    // pivots_nonpivots = (slong *)flint_malloc(sizeof(slong) * A->c);
    rank = _nmod_mat_rref(A, pivots_nonpivots, P);
    // flint_free(pivots_nonpivots);
    _perm_clear(P);

    return rank;
}


/* qinv[i]: column of mat with the pivot on row i of S; the pivot is implicitly 1 */
int spasm_flint_rref(nmod_mat_t mat, ulong *S, slong *qinv)
{
	double start = spasm_wtime();
    slong n = mat->r;
    slong m = mat->c;
    ulong prime = mat->mod.n;

	fprintf(stderr, "[flint/rref] Matrix of dimension %lld x %lld mod %" PRId64"... ", n, m, prime);
	fflush(stderr);
	slong rank = nmod_mat_rref_2(mat, qinv);

    // the ordering of entries may be wrong, since the pointers of rows
    // are not continuous in memory, so we need to reset it 
    // and we can perform the permutation of column on the matrix meanwhile

    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            S[i * m + j] = nmod_mat_entry(mat, i, qinv[j]);
        }
    }
    // abort();

	fprintf(stderr, "done in %.1fs. Rank %zd\n", spasm_wtime() - start, rank);
	return rank;
}


// TODO: add LU decomposition
// /*
//  * Q of size m (#cols). On output, Q[0:rank] is significant.
//  * It is implicit that U[i, 0:Q[i]] == 0 and U[i, 0:Q[i]] == 1.
//  * M[i, Q[i]:m] contain the actual data.
//  */
// static int spasm_flint_LU(nmod_mat_t mat, slong *p, slong *qinv)
// {
// 	double start = spasm_wtime();
//     ulong prime = mat->mod.n;
//     slong n = mat->r;
//     slong m = mat->c;
//     fprintf(stderr, "[ffpack/LU] Matrix of dimension %d x %d mod %" PRId64"... ", n, m, prime);
// 	fflush(stderr);
//     auto rank = nmod_mat_lu(p, mat, 0);
// 	// size_t rank = FFPACK::pPLUQ(GFp, FFLAS::FflasUnit, n, m, A, ldA, P, Qt);
// 	fprintf(stderr, "done in %.1fs. Rank %zd\n", spasm_wtime() - start, rank);
// 	start = spasm_wtime();
// 	// FFPACK::getEchelonForm (GFp, FFLAS::FflasUpper, FFLAS::FflasUnit, n, m, rank, Qt, A, ldA, FFPACK::FfpackTileRecursive);
// 	// FFPACK::getReducedEchelonForm(GFp, FFLAS::FflasUpper,  n, m, rank, Qt, A, ldA, FFPACK::FfpackTileRecursive);
// 	/* Qt is in LAPACK representation; convert */
// 	// LAPACKPerm2MathPerm (p, P, n);
//     for (size_t i = 0; i < m; i++)
//         qinv[i] = i;
    
// 	// LAPACKPerm2MathPerm (qinv, Qt, m);
// 	return rank;
// }


// spasm_ZZp spasm_datatype_read(const void *A, size_t i, spasm_datatype datatype)
// {
// 	return ((ulong *) A)[i];
// }

// void spasm_datatype_write(void *A, size_t i, spasm_datatype datatype, spasm_ZZp value)
// {
// 	((ulong *) A)[i] = value; return;
// }

// size_t spasm_datatype_size(spasm_datatype datatype)
// {
// 	return sizeof(ulong);
// }

// spasm_datatype spasm_datatype_choose(ulong prime)
// {
// 	return SPASM_ULONG;
// }

// const char * spasm_datatype_name(spasm_datatype datatype)
// {
// 	return "ulong";
// }