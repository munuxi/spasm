#include <assert.h>

#include "spasm.h"
#include "flint/ulong_extras.h"
#include "flint/nmod.h"

void spasm_field_init(ulong p, spasm_field F)
{
	F->p = p;
	F->halfp = p >> 1;
	// F->dinvp = n_precompute_inverse(p);
	F->invp = n_preinvert_limb(p);
	nmod_init(&(F->mod), p);
}

// static inline spasm_ZZp NORMALISE(const spasm_field F, i64 x)
// {
// 	if (x < F->mhalfp)
// 		x += F->p;
// 	else if (x > F->halfp)
// 		x -= F->p;
// 	return x;
// }


inline slong spasm_ZZp_signedrep(const spasm_field F, spasm_ZZp x)
{
	if (x > F->halfp){
		slong xx = x - (F->p);
		return xx;
	}
	return x;
}


inline spasm_ZZp spasm_ZZp_neg(const spasm_field F, spasm_ZZp a)
{
	return nmod_neg(a, F->mod);
}

inline spasm_ZZp spasm_ZZp_init(const spasm_field F, i64 x)
{
	if (x<0)
		return spasm_ZZp_neg(F,n_mod2_preinv((ulong)(-x), F->p, F->invp));
	else 
		return n_mod2_preinv((ulong)(x), F->p, F->invp);
	// return nmod_set_si(x, F->mod);
}


inline spasm_ZZp spasm_ZZp_add(const spasm_field F, spasm_ZZp a, spasm_ZZp b)
{
	return n_addmod(a, b, F->p);
}

inline spasm_ZZp spasm_ZZp_sub(const spasm_field F, spasm_ZZp a, spasm_ZZp b)
{
	return n_submod(a, b, F->p);
}

inline spasm_ZZp spasm_ZZp_mul(const spasm_field F, spasm_ZZp a, spasm_ZZp b)
{
	// i64 q = ((double) a) * ((double) b) * F->dinvp;
	// return NORMALISE(F, (i64) a * (i64) b - q * F->p);
	return nmod_mul(a, b, F->mod);
}

/* compute bezout relation u*a + v*p == 1; returns u */
// static i64 gcdext(i64 a, i64 p)
// {
// 	assert(a >= 0);
// 	i64 t = 0, u = 1;
// 	i64 r = p, s = a;
// 	while (s != 0) {
// 		i64 q = r / s;
// 		i64 foo = u;
// 		u = t - q * u;
// 		t = foo;

// 		i64 bar = s;
// 		s = r - q * s;
// 		r = bar;
// 	}
// 	return t;
// }

spasm_ZZp spasm_ZZp_inverse(const spasm_field F, spasm_ZZp a)
{
	// i64 aa = a;
	// if (aa < 0)
	// 	aa += F->p;
	// i64 inva = gcdext(aa, F->p);
	// return NORMALISE(F, inva);
	return n_invmod(a, F->p);
}


inline spasm_ZZp spasm_ZZp_axpy(const spasm_field F, spasm_ZZp a, spasm_ZZp x, spasm_ZZp y)
{
	// i64 q = (((((double) a) * ((double) x)) + (double) y) * F->dinvp);
	// i64 aa = a;
	// i64 xx = x;
	// i64 yy = y;
	// return NORMALISE(F, aa * xx + yy - q * F->p);
	return spasm_ZZp_add(F,spasm_ZZp_mul(F,a,x), y);
}
