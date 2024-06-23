#include <string>
#include "spasm.h"
#include "C/WolframLibrary.h"
#include "C/WolframSparseLibrary.h"

// #include <chrono>
// using namespace std::chrono;

EXTERN_C DLLEXPORT int modrref(WolframLibraryData ld, mint Argc, MArgument *Args, MArgument Res) {
    auto mat = MArgument_getMSparseArray(Args[0]);
    auto p = MArgument_getInteger(Args[1]);
    auto nthreads = MArgument_getInteger(Args[2]);

    auto sf = ld->sparseLibraryFunctions;

    auto ranks = sf->MSparseArray_getRank(mat);
    if (ranks != 2 && sf->MSparseArray_getImplicitValue(mat) != 0)
        return LIBRARY_FUNCTION_ERROR;
    
    auto dims = sf->MSparseArray_getDimensions(mat);
    auto nrows = dims[0];
    auto ncols = dims[1];

    auto m_rowptr = sf->MSparseArray_getRowPointers(mat);
    auto m_colptr = sf->MSparseArray_getColumnIndices(mat);
    auto m_valptr = sf->MSparseArray_getExplicitValues(mat);

    // rowptr, valptr, colptr are managed by mathematica
    // do not free them
    mint* rowptr = ld->MTensor_getIntegerData(*m_rowptr);
    mint* valptr = ld->MTensor_getIntegerData(*m_valptr);
    mint* colptr = ld->MTensor_getIntegerData(*m_colptr);

    auto nnz = rowptr[nrows];

    flint_set_num_threads(nthreads);
    spasm_field Fp;
    spasm_field_init(p,Fp);
    struct echelonize_opts dopts;
	auto opts = &dopts;
	// std::cout << "[echelonize] using default settings\n";
	spasm_echelonize_init_opts(opts);

    struct spasm_triplet *T = spasm_triplet_alloc(nrows, ncols, nnz, p, true);
    for (int i = 0; i < nrows; i++) 
        for (int k = rowptr[i]; k < rowptr[i + 1]; k++) 
            spasm_add_entry(T, i, colptr[k] - 1, spasm_ZZp_init(Fp, valptr[k]));
    
    auto mat2 = spasm_compress(T);
	spasm_triplet_free(T);

    auto fact = spasm_echelonize(mat2, opts);
    int *Rqinv = (int*)(spasm_malloc(ncols * sizeof(int)));
	struct spasm_csr *A = spasm_rref(fact, Rqinv);
    spasm_lu_free(fact);

    nnz = spasm_nnz(A);
    MTensor result;
    mint dims_result[2] = {nnz+1, 3};
    ld->MTensor_new(MType_Integer, 2, dims_result, &result);

    // output now
    // mathematica is 1-based
    mint entry[2];
	const int *Aj = A->j;
	const i64 *Ap = A->p;
	const spasm_ZZp *Ax = A->x;
	int n = A->n;
    // int m = nrows - n;
    entry[0] = 0;
    int err;
	for (int i = 0; i < n; i++){
		for (i64 px = Ap[i]; px < Ap[i + 1]; px++) {
            slong x = (Ax != NULL) ? spasm_ZZp_signedrep(A->field,Ax[px]) : 1;
             // (i + 1, Aj[px] + 1, x);
            entry[0]++;
            entry[1] = 1; // row pos
			err = ld->MTensor_setInteger(result, entry, i + 1);
            entry[1] = 2; // col pos
            err = ld->MTensor_setInteger(result, entry, Aj[px] + 1); 
            entry[1] = 3; // value pos
            err = ld->MTensor_setInteger(result, entry, x);
		}
    }
    // last is dimension 
    entry[0]++;
    entry[1] = 1; // row pos
	err = ld->MTensor_setInteger(result, entry, n);
    entry[1] = 2; // col pos
    err = ld->MTensor_setInteger(result, entry, ncols); 
    entry[1] = 3; // value pos
    err = ld->MTensor_setInteger(result, entry, 0);

    // MArgument_setInteger(Res, num);
    MArgument_setMTensor(Res, result);
    spasm_csr_free(A);
    free(Rqinv);

    return LIBRARY_NO_ERROR;
}
