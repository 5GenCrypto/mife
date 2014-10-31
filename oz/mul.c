#include "mul.h"
#include "oz.h"
#include "util.h"
#include "flint-addons.h"

void fmpz_poly_oz_rem(fmpz_poly_t rem, const fmpz_poly_t f, const long n) {
  fmpz_poly_set(rem, f);
  fmpz_t lead; fmpz_init(lead);
  fmpz_t coef; fmpz_init(coef);
  for(long d=fmpz_poly_degree(rem); d>=n; d--) {
    fmpz_poly_get_coeff_fmpz(lead, rem, d);
    fmpz_poly_get_coeff_fmpz(coef, rem, d-n);
    fmpz_sub(coef, coef, lead);
    fmpz_poly_set_coeff_si(rem, d, 0);
    fmpz_poly_set_coeff_fmpz(rem, d-n, coef);
  }
  fmpz_clear(coef);
  fmpz_clear(lead);
}

void fmpz_mod_poly_oz_rem(fmpz_mod_poly_t rem, const fmpz_mod_poly_t f, const long n) {
  fmpz_mod_poly_set(rem, f);
  fmpz_t lead; fmpz_init(lead);
  fmpz_t coef; fmpz_init(coef);
  for(long d=fmpz_mod_poly_degree(rem); d>=n; d--) {
    fmpz_mod_poly_get_coeff_fmpz(lead, rem, d);
    fmpz_mod_poly_get_coeff_fmpz(coef, rem, d-n);
    fmpz_sub(coef, coef, lead);
    fmpz_mod_poly_set_coeff_ui(rem, d, 0);
    fmpz_mod_poly_set_coeff_fmpz(rem, d-n, coef);
  }
  fmpz_clear(coef);
  fmpz_clear(lead);
}

void fmpq_poly_oz_rem(fmpq_poly_t rem, const fmpq_poly_t f, const long n) {
  const long d = fmpq_poly_degree(f);
  if(d > 2*n-1) {
    // TODO: Don't be so lazy and take care of this case as well
    fmpq_poly_t modulus;
    fmpq_poly_init_oz_modulus(modulus, n);
    fmpq_poly_rem(rem, f, modulus);
    fmpq_poly_clear(modulus);
  } else if (d>=n) {
    assert(d-n+1 <= n);
    fmpq_poly_t r;
    fmpq_poly_init(r);
    fmpq_poly_realloc(r, n);
    mpq_t *z = (mpq_t*)calloc(n, sizeof(mpq_t));
    for (int i = 0; i < n; i++) {
      mpq_init(z[i]);
      fmpq_poly_get_coeff_mpq(z[i], f, i);
    }
    fmpq_poly_set_array_mpq(r, (const mpq_t*)z, n);

    for (int i=0; i<d-n+1; i++)
      fmpq_poly_get_coeff_mpq(z[i], f, n+i);
    fmpq_poly_t t;
    fmpq_poly_init(t);
    fmpq_poly_set_array_mpq(t, (const mpq_t*)z, d-n+1);
    fmpq_poly_sub(r, r, t);
    fmpq_poly_set(rem, r);

    fmpq_poly_clear(t);
    fmpq_poly_clear(r);
    for (int i=0; i<n; i++)
      mpq_clear(z[i]);
    free(z);
  }
  return;
}

int _fmpz_nth_root(fmpz_t rop, const long n, const fmpz_t q) {
  if (fmpz_cmp_si(q, 2) == 0) {
    fmpz_set_ui(rop, 1);
    return 1;
  }

  fmpz_t qm1;
  fmpz_init_set(qm1, q);
  fmpz_sub_ui(qm1, qm1, 1);

  if(!fmpz_divisible_si(qm1, n)) {
    fmpz_clear(qm1);
    return 0;
  }

  fmpz_t a;  fmpz_init_set_ui(a, 2);
  fmpz_t e;  fmpz_init(e);
  fmpz_divexact_ui(e, qm1, 2);

  int found = 0;
  fmpz_t b;  fmpz_init(b);
  while(fmpz_cmp(a, q) < 0) {
    fmpz_powm(b, a, e, q);
    if (fmpz_cmp_ui(b, 1) != 0) {
      found = 1;
      break;
    }
    fmpz_add_ui(a, a, 1);
  }
  fmpz_clear(b);

  if(!found) {
    fmpz_clear(qm1);
    fmpz_clear(a);
    fmpz_clear(e);
    return 0;
  }
  fmpz_divexact_ui(e, qm1, n);
  fmpz_powm(rop, a, e, q);
  fmpz_clear(qm1);
  fmpz_clear(a);
  fmpz_clear(e);
  return 1;
}

void _fmpz_mod_poly_oz_ntt(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, const fmpz_mod_poly_t w, const size_t n) {
  if (n == 1) {
    fmpz_mod_poly_set(rop, op);
    return;
  }
  const fmpz *q = fmpz_mod_poly_modulus(op);
  const size_t k = n_flog(n,2);

  fmpz_mod_poly_t a; fmpz_mod_poly_init2(a, q, n);
  for (size_t i = 0; i < n; i++) {
    size_t ii = i, r=0;
    for (size_t h = 0; h < k; h++) {
      r = (r << 1) | (ii & 1);
      ii >>= 1;
    }
    fmpz_set(a->coeffs+r, op->coeffs+i);
  }
  fmpz_mod_poly_t b; fmpz_mod_poly_init2(b, q, n);

  fmpz *a_hdl = a->coeffs;
  fmpz *b_hdl = b->coeffs;

  fmpz_t tmp;
  fmpz_init(tmp);
  for(size_t i=0; i<k; i++) {
    for(size_t j=0; j<n/2; j++) {
      const size_t tk  = (1UL<<(k-1-i));
      const size_t pij = (j/tk) * tk;
      fmpz_mul(tmp, a_hdl + 2*j+1, w->coeffs + pij);
      if(pij)
        fmpz_mod(tmp, tmp, q);
      fmpz_add(b_hdl + j,       a_hdl + 2*j, tmp);
      if (fmpz_cmpabs(b_hdl + j, q) >= 0)
        fmpz_sub(b_hdl + j, b_hdl + j, q);
      fmpz_sub(b_hdl + j + n/2, a_hdl + 2*j, tmp);
      if (fmpz_sgn(b_hdl + j + n/2) < 0)
        fmpz_add(b_hdl + j + n/2, b_hdl + j + n/2, q);
    }
    if(i!=k-1) {
      fmpz *tmp = a_hdl;
      a_hdl = b_hdl;
      b_hdl = tmp;
    }
  }
  fmpz_mod_poly_realloc(rop, n);
  _fmpz_vec_set(rop->coeffs, b_hdl, n);
  rop->length = n;
  fmpz_mod_poly_clear(b);
  fmpz_mod_poly_clear(a);
}

void fmpz_mod_poly_oz_set_powers(fmpz_mod_poly_t op, const size_t n, const fmpz_t w) {
  fmpz_mod_poly_realloc(op, n);

  fmpz_t acc; fmpz_init_set_ui(acc, 1);
  fmpz_set(op->coeffs + 0, acc);
  for(size_t i=1; i<n; i++) {
    fmpz_mul(acc, acc, w);
    fmpz_mod(acc, acc, fmpz_mod_poly_modulus(op));
    fmpz_set(op->coeffs + i, acc);
  }
  op->length = n;
  _fmpz_mod_poly_normalise(op);
}

void fmpz_mod_poly_oz_ntt(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, const size_t n) {
  const fmpz *q = fmpz_mod_poly_modulus(op);
  fmpz_t w;  fmpz_init(w);
  if (!_fmpz_nth_root(w, n, q)) {
    fmpz_clear(w);
    oz_die("q does not have a n-th root of unity");
  }
  fmpz_mod_poly_t wvec; fmpz_mod_poly_init2(wvec, q, n);
  fmpz_mod_poly_oz_set_powers(wvec, n, w);
  fmpz_clear(w);

  _fmpz_mod_poly_oz_ntt(rop, op, wvec, n);
  fmpz_mod_poly_clear(wvec);
}

void fmpz_mod_poly_oz_intt(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, const size_t n) {
  const fmpz *q = fmpz_mod_poly_modulus(op);
  fmpz_t w;  fmpz_init(w);
  if (!_fmpz_nth_root(w, n, q)) {
    fmpz_clear(w);
    oz_die("q does not have a n-th root of unity");
  }
  fmpz_invmod(w, w, fmpz_mod_poly_modulus(op));
  fmpz_mod_poly_t wvec; fmpz_mod_poly_init2(wvec, q, n);
  fmpz_mod_poly_oz_set_powers(wvec, n, w);
  fmpz_clear(w);

  _fmpz_mod_poly_oz_ntt(rop, op, wvec, n);
  fmpz_mod_poly_clear(wvec);
}

void fmpz_mod_poly_oz_ntt_precomp_init(fmpz_mod_poly_oz_ntt_precomp_t op, const size_t n, const fmpz_t q) {
  fmpz_t w;  fmpz_init(w);
  if (!_fmpz_nth_root(w, n, q)) {
    fmpz_clear(w);
    oz_die("q does not have a n-th root of unity");
  }

  fmpz_t phi;  fmpz_init(phi);
  if(!fmpz_sqrtmod(phi, w, q)) {
    fmpz_clear(phi);
    oz_die("q does not have a 2n-th root of unity");
  }

  op->n = n;

  fmpz_mod_poly_init2(op->w, q, n);
  fmpz_mod_poly_oz_set_powers(op->w, n, w);

  fmpz_invmod(w, w, q);
  fmpz_mod_poly_init2(op->w_inv, q, n);
  fmpz_mod_poly_oz_set_powers(op->w_inv, n, w);
  fmpz_clear(w);

  fmpz_mod_poly_init2(op->phi, q, n);
  fmpz_mod_poly_oz_set_powers(op->phi, n, phi);

  fmpz_invmod(phi, phi, q);
  fmpz_mod_poly_init2(op->phi_inv, q, n);
  fmpz_mod_poly_oz_set_powers(op->phi_inv, n, phi);
  fmpz_clear(phi);

  /** We fold 1/n into φ^-1 **/
  fmpz_t n_inv;  fmpz_init_set_ui(n_inv, n);
  fmpz_invmod(n_inv, n_inv, q);

  for(size_t i=0; i<n; i++) {
    fmpz_mul(op->phi_inv->coeffs+i, n_inv, op->phi_inv->coeffs+i);
    fmpz_mod(op->phi_inv->coeffs+i, op->phi_inv->coeffs+i, q);
  }
  fmpz_clear(n_inv);
}

void fmpz_mod_poly_oz_ntt_precomp_clear(fmpz_mod_poly_oz_ntt_precomp_t op) {
  fmpz_mod_poly_clear(op->w);
  fmpz_mod_poly_clear(op->w_inv);
  fmpz_mod_poly_clear(op->phi);
  fmpz_mod_poly_clear(op->phi_inv);
}

void _fmpz_mod_poly_oz_ntt_mul(fmpz_mod_poly_t h, const fmpz_mod_poly_t f, const fmpz_mod_poly_t g, const size_t n) {
  const fmpz *q = fmpz_mod_poly_modulus(f);
  fmpz_mod_poly_realloc(h, n);

  for(size_t i=0; i<n; i++) {
    fmpz_mul(h->coeffs + i, f->coeffs + i, g->coeffs + i);
    fmpz_mod(h->coeffs + i, h->coeffs+i, q);
  }
  h->length = n;
}


void _fmpz_mod_poly_oz_ntt_set_ui(fmpz_mod_poly_t op, const unsigned long c, const size_t n) {
  fmpz_mod_poly_realloc(op, n);
  for(size_t i=0; i<n; i++)
    fmpz_mod_poly_set_coeff_ui(op, i, c);
}

void _fmpz_mod_poly_oz_ntt_pow_ui(fmpz_mod_poly_t rop, const fmpz_mod_poly_t f, unsigned long e, const size_t n) {
  _fmpz_mod_poly_oz_ntt_set_ui(rop, 1, n);

  fmpz_mod_poly_t tmp;
  fmpz_mod_poly_init2(tmp, fmpz_mod_poly_modulus(f), n);
  fmpz_mod_poly_set(tmp, f);
  while(e>0) {
    if (e&1)
      _fmpz_mod_poly_oz_ntt_mul(rop, rop, f, n);
    e = e>>1;
    _fmpz_mod_poly_oz_ntt_mul(tmp, tmp, tmp, n);
  }
  fmpz_mod_poly_clear(tmp);
}

void _fmpz_mod_poly_oz_ntt_enc(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, const fmpz_mod_poly_oz_ntt_precomp_t precomp) {
  const fmpz *q = fmpz_mod_poly_modulus(op);
  fmpz_mod_poly_realloc(rop, precomp->n);
  for(size_t i=0; i<precomp->n; i++) {
    fmpz_mul(rop->coeffs+i, precomp->phi->coeffs+i, op->coeffs+i);
    fmpz_mod(rop->coeffs+i, rop->coeffs+i, q);
  }
  rop->length = precomp->n;
  _fmpz_mod_poly_oz_ntt(rop, rop, precomp->w, precomp->n);
}

void _fmpz_mod_poly_oz_ntt_dec(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, const fmpz_mod_poly_oz_ntt_precomp_t precomp) {
  const fmpz *q = fmpz_mod_poly_modulus(op);

  _fmpz_mod_poly_oz_ntt(rop, op, precomp->w_inv, precomp->n);

  for(size_t i=0; i<precomp->n; i++) {
    fmpz_mul(rop->coeffs+i, precomp->phi_inv->coeffs+i, rop->coeffs+i);
    fmpz_mod(rop->coeffs+i, rop->coeffs+i, q);
  }
}

void _fmpz_mod_poly_oz_mul_nttnwc(fmpz_mod_poly_t h, const fmpz_mod_poly_t f, const fmpz_mod_poly_t g, const fmpz_mod_poly_oz_ntt_precomp_t precomp) {
  const size_t n = precomp->n;
  const fmpz *q = fmpz_mod_poly_modulus(f);

  /* printf("w: "); fmpz_print(w); printf(", phi: "); fmpz_print(phi); printf("\n"); */

  fmpz_mod_poly_t F;  fmpz_mod_poly_init2(F, q, n);
  fmpz_mod_poly_t G;  fmpz_mod_poly_init2(G, q, n);

  _fmpz_mod_poly_oz_ntt_enc(F, f, precomp);
  _fmpz_mod_poly_oz_ntt_enc(G, g, precomp);

  _fmpz_mod_poly_oz_ntt_mul(h, F, G, n);

  fmpz_mod_poly_clear(F);
  fmpz_mod_poly_clear(G);

  _fmpz_mod_poly_oz_ntt_dec(h, h, precomp);
}

void fmpz_mod_poly_oz_mul_nttnwc(fmpz_mod_poly_t h, const fmpz_mod_poly_t f, const fmpz_mod_poly_t g, const size_t n) {
  assert(fmpz_equal(fmpz_mod_poly_modulus(f), fmpz_mod_poly_modulus(g)));
  assert(fmpz_mod_poly_length(f) == fmpz_mod_poly_length(g));
  assert(fmpz_mod_poly_length(f) == (long)n);

  const fmpz *q = fmpz_mod_poly_modulus(f);
  fmpz_mod_poly_oz_ntt_precomp_t precomp;
  fmpz_mod_poly_oz_ntt_precomp_init(precomp, n, q);

  _fmpz_mod_poly_oz_mul_nttnwc(h, f, g, precomp);

  fmpz_mod_poly_oz_ntt_precomp_clear(precomp);
}
