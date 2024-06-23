#ifndef _SPASM_H
#define _SPASM_H

#include <stddef.h>           // size_t
#include <inttypes.h>         // int64_t
#include <stdio.h>            // FILE
#include <stdbool.h>

typedef uint8_t u8;
typedef int64_t i64;
typedef uint64_t u64;
typedef uint32_t u32;
typedef int32_t i32;

#ifdef _OPENMP
#include <omp.h>
#endif

#include "flint/nmod.h"
#include "flint/nmod_mat.h"

#define SPASM_VERSION "1.3"
#define SPASM_BUG_ADDRESS "<charles.bouillaguet@lip6.fr>"

#if defined(_MSC_VER)
//  Microsoft 
#define EXPORT __declspec(dllexport)
#define IMPORT __declspec(dllimport)
#elif defined(__GNUC__)
//  GCC
#define EXPORT __attribute__((visibility("default")))
#define IMPORT
#else
//  do nothing and hope for the best?
#define EXPORT
#define IMPORT
#pragma warning Unknown dynamic link import/export semantics.
#endif

#define errx(code, fmt, ...) do { fprintf(stderr, fmt "\n", ##__VA_ARGS__); exit(code); } while (0)
#define err(code, fmt, ...) do { fprintf(stderr, fmt "\n", ##__VA_ARGS__); exit(code); } while (0)
#define warn(fmt, ...) fprintf(stderr, fmt "\n", ##__VA_ARGS__)
/* --- primary struct spasm_csr routines and data structures --- */

// unfortunately we use "n" for #rows and "m" for #columns whereas the rest of the world (BLAS...)
// does the opposite... 

typedef ulong spasm_ZZp;

struct spasm_field_struct {
	ulong p;
    ulong halfp;
    // ulong mhalfp;
	nmod_t mod;
    // double dinvp;
	ulong invp;
};
typedef struct spasm_field_struct spasm_field[1];

struct spasm_csr {                /* matrix in compressed-sparse row format */
	i64 nzmax;                    /* maximum number of entries */
	int n;                        /* number of rows */
	int m;                        /* number of columns */
	i64 *p;                       /* row pointers (size n+1) */
	int *j;                       /* column indices, size nzmax */
	spasm_ZZp *x;                 /* numerical values, size nzmax (optional) */
	spasm_field field;
	/*
	 * The actual number of entries is p[n]. 
	 * Coefficients of a row need not be sorted by column index.
	 * The numerical values are optional (useful for storing a sparse graph, or the pattern of a matrix).
	 */
};

struct spasm_triplet {             /* matrix in triplet form */
	i64 nzmax;                     /* maximum number of entries */
	i64 nz;                        /* # entries */
	int n;                         /* number of rows */
	int m;                         /* number of columns */
	int *i;                        /* row indices, size nzmax */
	int *j;                        /* column indices (size nzmax) */
	spasm_ZZp *x;                  /* numerical values, size nzmax (optional) */
	spasm_field field;
};

struct spasm_lu {                  /* a PLUQ factorisation */
	int r;                         /* rank of the input matrix */
	bool complete;                 /* if L != NULL, indicates whether A == L*U */
	struct spasm_csr *L;
	struct spasm_csr *U;
	int *qinv;                     /* locate pivots in U (on column j, row qinv[j]) */
	int *p;                        /* locate pivots in L (on column j, row p[j]) */
	struct spasm_triplet *Ltmp;           /* for internal use during the factorization */
};

struct spasm_dm {      /**** a Dulmage-Mendelson decomposition */
				int *p;       /* size n, row permutation */
				int *q;       /* size m, column permutation */
				int *r;       /* size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q) */
				int *c;       /* size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q) */
				int nb;       /* # of blocks in fine decomposition */
				int rr[5];    /* coarse row decomposition */
				int cc[5];    /* coarse column decomposition */
};

struct echelonize_opts {
	/* pivot search sub-algorithms */
	bool enable_greedy_pivot_search;

	/* echelonization sub-algorithms */
	bool enable_tall_and_skinny;
	bool enable_dense;
	bool enable_GPLU;

	/* Parameters of the "root" echelonization procedure itself */
	bool L;                         /* should we compute L / Lp in addition to U / Uqinv ? */
	bool complete;                  /* A == LU / otherwise L is just OK for the pivotal rows */
	double min_pivot_proportion;    /* minimum number of pivots found to keep going; < 0 = keep going */
	int max_round;                  /* maximum number of rounds; < 0 = keep going */
 
 	/* Parameters that determine the choice of a finalization strategy */
 	double sparsity_threshold;      /* denser than this --> dense method; < 0 = keep going */

	/* options of dense methods */
	int dense_block_size;           /* #rows processed in each batch; determine memory consumption */
	double low_rank_ratio;          /* if k rows have rank less than k * low_rank_ratio --> "tall-and-skinny"; <0 = don't */
	double tall_and_skinny_ratio;   /* aspect ratio (#rows / #cols) higher than this --> "tall-and-skinny"; <0 = don't */
	double low_rank_start_weight;   /* compute random linear combinations of this many rows; -1 = auto-select */

};

struct spasm_rank_certificate {
	int r;
	i64 prime;
	u8 hash[32];
	int *i;           /* size r */
	int *j;           /* size r */
	spasm_ZZp *x;     /* size r */
	spasm_ZZp *y;     /* size r */
};

typedef struct {
    u32 h[8];
    u32 Nl, Nh;
    u32 data[16];
    u32 num, md_len;
} spasm_sha256_ctx;

typedef struct {
        u32 block[11];   /* block[0:8] == H(matrix); block[8] = prime; block[9] = ctx, block[10] = seq */
        u32 hash[8];
        u32 prime;
        u32 mask;        /* 2^i - 1 where i is the smallest s.t. 2^i > prime */
        int counter;
        int i;
        spasm_field field;
} spasm_prng_ctx;

// typedef enum {SPASM_DOUBLE, SPASM_FLOAT, SPASM_ULONG} spasm_datatype;

#define SPASM_IDENTITY_PERMUTATION NULL
#define SPASM_IGNORE NULL
#define SPASM_IGNORE_VALUES 0

#ifdef __cplusplus
extern "C" {
#endif

/* spasm_ZZp.c */
EXPORT void spasm_field_init(ulong p, spasm_field F);
EXPORT spasm_ZZp spasm_ZZp_init(const spasm_field F, i64 x);
EXPORT slong spasm_ZZp_signedrep(const spasm_field F, spasm_ZZp x);
EXPORT spasm_ZZp spasm_ZZp_add(const spasm_field F, spasm_ZZp a, spasm_ZZp b);
EXPORT spasm_ZZp spasm_ZZp_neg(const spasm_field F, spasm_ZZp a);
EXPORT spasm_ZZp spasm_ZZp_sub(const spasm_field F, spasm_ZZp a, spasm_ZZp b);
EXPORT spasm_ZZp spasm_ZZp_mul(const spasm_field F, spasm_ZZp a, spasm_ZZp b);
EXPORT spasm_ZZp spasm_ZZp_inverse(const spasm_field F, spasm_ZZp a);
EXPORT spasm_ZZp spasm_ZZp_axpy(const spasm_field F, spasm_ZZp a, spasm_ZZp x, spasm_ZZp y);

/* sha256.c */
EXPORT void spasm_SHA256_init(spasm_sha256_ctx *c);
EXPORT void spasm_SHA256_update(spasm_sha256_ctx *c, const void *data, size_t len);
EXPORT void spasm_SHA256_final(u8 *md, spasm_sha256_ctx *c);

/* spasm_prng.c */
// void spasm_prng_seed(const u8 *seed, i64 prime, u32 seq, spasm_prng_ctx *ctx);
// void spasm_prng_seed_simple(i64 prime, u64 seed, u32 seq, spasm_prng_ctx *ctx);
// u32 spasm_prng_u32(spasm_prng_ctx *ctx);
// spasm_ZZp spasm_prng_ZZp(spasm_prng_ctx *ctx);

/* spasm_util.c */
EXPORT double spasm_wtime();
EXPORT i64 spasm_nnz(const struct spasm_csr * A);
EXPORT void *spasm_malloc(i64 size);
EXPORT void *spasm_calloc(i64 count, i64 size);
EXPORT void *spasm_realloc(void *ptr, i64 size);
EXPORT struct spasm_csr *spasm_csr_alloc(int n, int m, i64 nzmax, i64 prime, bool with_values);
EXPORT void spasm_csr_realloc(struct spasm_csr * A, i64 nzmax);
EXPORT void spasm_csr_resize(struct spasm_csr * A, int n, int m);
EXPORT void spasm_csr_free(struct spasm_csr * A);
EXPORT struct spasm_triplet *spasm_triplet_alloc(int m, int n, i64 nzmax, i64 prime, bool with_values);
EXPORT void spasm_triplet_realloc(struct spasm_triplet * A, i64 nzmax);
EXPORT void spasm_triplet_free(struct spasm_triplet * A);
EXPORT struct spasm_dm *spasm_dm_alloc(int n, int m);
EXPORT void spasm_dm_free(struct spasm_dm * P);
EXPORT void spasm_lu_free(struct spasm_lu *N);
EXPORT void spasm_human_format(int64_t n, char *target);
EXPORT int spasm_get_num_threads();
EXPORT int spasm_get_thread_num();
static inline ulong spasm_get_prime(const struct spasm_csr *A) { return A->field->p; }

/* spasm_triplet.c */
EXPORT void spasm_add_entry(struct spasm_triplet *T, int i, int j, i64 x);
EXPORT void spasm_triplet_transpose(struct spasm_triplet * T);
EXPORT struct spasm_csr *spasm_compress(const struct spasm_triplet * T);

/* spasm_io.c */
EXPORT struct spasm_triplet *spasm_triplet_load(FILE * f, i64 prime, u8 *hash);
EXPORT void spasm_triplet_save(const struct spasm_triplet * A, FILE * f);
EXPORT void spasm_csr_save(const struct spasm_csr * A, FILE * f);
EXPORT void spasm_save_pnm(const struct spasm_csr * A, FILE * f, int x, int y, int mode, struct spasm_dm *DM);

/* spasm_transpose.c */
EXPORT struct spasm_csr *spasm_transpose(const struct spasm_csr * C, int keep_values);

/* spasm_submatrix.c */
EXPORT struct spasm_csr *spasm_submatrix(const struct spasm_csr * A, int r_0, int r_1, int c_0, int c_1, int with_values);

/* spasm_permutation.c */
EXPORT void spasm_pvec(const int *p, const spasm_ZZp * b, spasm_ZZp * x, int n);
EXPORT void spasm_ipvec(const int *p, const spasm_ZZp * b, spasm_ZZp * x, int n);
EXPORT int *spasm_pinv(int const *p, int n);
EXPORT struct spasm_csr *spasm_permute(const struct spasm_csr * A, const int *p, const int *qinv, int with_values);
EXPORT int *spasm_random_permutation(int n);
EXPORT void spasm_range_pvec(int *x, int a, int b, int *p);

/* spasm_scatter.c */
EXPORT void spasm_scatter(const struct spasm_csr *A, int i, const spasm_ZZp beta, spasm_ZZp * x);

/* spasm_reach.c */
EXPORT int spasm_dfs(int i, const struct spasm_csr * G, int top, int *xi, int *pstack, int *marks, const int *pinv);
EXPORT int spasm_reach(const struct spasm_csr * A, const struct spasm_csr * B, int k, int l, int *xj, const int *qinv);

/* spasm_spmv.c */
EXPORT void spasm_xApy(const spasm_ZZp *x, const struct spasm_csr *A, spasm_ZZp *y);
EXPORT void spasm_Axpy(const struct spasm_csr *A, const spasm_ZZp *x, spasm_ZZp *y);

/* spasm_triangular.c */
EXPORT void spasm_dense_back_solve(const struct spasm_csr *L, spasm_ZZp *b, spasm_ZZp *x, const int *p);
EXPORT bool spasm_dense_forward_solve(const struct spasm_csr * U, spasm_ZZp * b, spasm_ZZp * x, const int *q);
EXPORT int spasm_sparse_triangular_solve(const struct spasm_csr *U, const struct spasm_csr *B, int k, int *xj, spasm_ZZp * x, const int *qinv);

/* spasm_schur.c */
EXPORT struct spasm_csr *spasm_schur(const struct spasm_csr *A, const int *p, int n, const struct spasm_lu *fact,
                   double est_density, struct spasm_triplet *L, const int *p_in, int *p_out);
EXPORT double spasm_schur_estimate_density(const struct spasm_csr * A, const int *p, int n, const struct spasm_csr *U, const int *qinv, int R);
EXPORT void spasm_schur_dense(const struct spasm_csr *A, const int *p, int n, const int *p_in,
	struct spasm_lu *fact, ulong *S, int *q, int *p_out);
EXPORT void spasm_schur_dense_randomized(const struct spasm_csr *A, const int *p, int n, const struct spasm_csr *U, const int *qinv,
	ulong *S, int *q, int N, int w);

/* spasm_pivots.c */
EXPORT int spasm_pivots_extract_structural(const struct spasm_csr *A, const int *p_in, struct spasm_lu *fact, int *p, struct echelonize_opts *opts);

/* spasm_matching.c */
EXPORT int spasm_maximum_matching(const struct spasm_csr *A, int *jmatch, int *imatch);
EXPORT int *spasm_permute_row_matching(int n, const int *jmatch, const int *p, const int *qinv);
EXPORT int *spasm_permute_column_matching(int m, const int *imatch, const int *pinv, const int *q);
EXPORT int *spasm_submatching(const int *match, int a, int b, int c, int d);
EXPORT int spasm_structural_rank(const struct spasm_csr *A);

/* spasm_dm.c */
EXPORT struct spasm_dm *spasm_dulmage_mendelsohn(const struct spasm_csr *A);

/* spasm_scc.c */
EXPORT struct spasm_dm *spasm_strongly_connected_components(const struct spasm_csr *A);

/* spasm_flint.c */
EXPORT int spasm_flint_rref(nmod_mat_t mat, ulong* S, slong *qinv);
// int spasm_ffpack_LU(i64 prime, int n, int m, void *A, int ldA, spasm_datatype datatype, size_t *p, size_t *qinv);
// spasm_ZZp spasm_datatype_read(const void *A, size_t i, spasm_datatype datatype);
// void spasm_datatype_write(void *A, size_t i, spasm_datatype datatype, spasm_ZZp value);
// size_t spasm_datatype_size(spasm_datatype datatype);
// spasm_datatype spasm_datatype_choose(ulong prime);
// const char * spasm_datatype_name(spasm_datatype datatype);

/* spasm_echelonize */
EXPORT void spasm_echelonize_init_opts(struct echelonize_opts *opts);
EXPORT struct spasm_lu* spasm_echelonize(const struct spasm_csr *A, struct echelonize_opts *opts);

/* spasm_rref.c */
EXPORT struct spasm_csr * spasm_rref(const struct spasm_lu *fact, int *Rqinv);

/* spasm_kernel.c */
EXPORT struct spasm_csr * spasm_kernel(const struct spasm_lu *fact);
EXPORT struct spasm_csr * spasm_kernel_from_rref(const struct spasm_csr *R, const int *qinv);

/* spasm_solve.c */
EXPORT bool spasm_solve(const struct spasm_lu *fact, const spasm_ZZp *b, spasm_ZZp *x);
EXPORT struct spasm_csr * spasm_gesv(const struct spasm_lu *fact, const struct spasm_csr *B, bool *ok);

/* spasm_certificate.c */
// struct spasm_rank_certificate * spasm_certificate_rank_create(const struct spasm_csr *A, const u8 *hash, const struct spasm_lu *fact);
// bool spasm_certificate_rank_verify(const struct spasm_csr *A, const u8 *hash, const struct spasm_rank_certificate *proof);
// void spasm_rank_certificate_save(const struct spasm_rank_certificate *proof, FILE *f);
// bool spasm_rank_certificate_load(FILE *f, struct spasm_rank_certificate *proof);
// bool spasm_factorization_verify(const struct spasm_csr *A, const struct spasm_lu *fact, u64 seed);


/* utilities */
static inline int spasm_max(int a, int b)
{
	return (a > b) ? a : b;
}

static inline int spasm_min(int a, int b)
{
	return (a < b) ? a : b;
}

static inline int spasm_row_weight(const struct spasm_csr * A, int i)
{
	i64 *Ap = A->p;
	return Ap[i + 1] - Ap[i];
}
#endif

#ifdef __cplusplus
}
#endif