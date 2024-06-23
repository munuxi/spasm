#include "spasm.h"
#include "flint/ulong_extras.h"
/*
 * x = x + beta * A[i], where x is a dense vector
 * 
 * This is where all the heavy lifting should take place.
 */
void spasm_scatter(const struct spasm_csr *A, int i, const spasm_ZZp beta, spasm_ZZp * x)
{
	const i64 *Ap = A->p;
	const int *Aj = A->j;
	const spasm_ZZp *Ax = A->x;
	for (i64 px = Ap[i]; px < Ap[i + 1]; px++) {
		int j = Aj[px];
		ulong betainv = n_mulmod_precomp_shoup(beta, A->field->p);
		x[j] = spasm_ZZp_add(A->field, 
			n_mulmod_shoup(beta, Ax[px], betainv, A->field->p), x[j]);
	}
}
