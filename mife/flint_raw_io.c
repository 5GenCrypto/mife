#include "flint_raw_io.h"

#define xxx_putc(c)        \
do {                       \
    z = fputc((c), file);  \
    if (z <= 0)            \
        return z;          \
} while (0)

#define xxx_flint_printf()                       \
do {                                       \
    z = flint_fprintf(file, "%li %li", r, c);  \
    if (z <= 0)                            \
        return z;                          \
} while (0)

#define xxx_fmpz_print(f)        \
do {                             \
    z = fmpz_out_raw(file, (f));  \
    if (z <= 0)                  \
        return z;                \
} while(0)

int fmpz_mat_fprint_raw(FILE * file, const fmpz_mat_t mat)
{
    int z;
    slong i, j;
    slong r = mat->r;
    slong c = mat->c;

    xxx_flint_printf();
    for (i = 0; (i < r); i++)
    {
        for (j = 0; j < c; j++)
        {
            xxx_fmpz_print(mat->rows[i] + j);
            if (j != c - 1)
                xxx_putc(' ');
        }
        if (i != r - 1)
            xxx_putc(' ');
    }

    return z;
}

#undef xxx_putc
#undef xxx_flint_printf
#undef xxx_fmpz_print

int fmpz_mat_fread_raw(FILE* file, fmpz_mat_t mat) {
    slong r, c, i, j;
    int byte_count;
    mpz_t t;

    /* first number in file should be row dimension */
    mpz_init(t);
    byte_count = mpz_inp_str(t, file, 10);
    if (byte_count == 0)
    {
        mpz_clear(t);
        return 0;
    }
    
    if (!mpz_fits_slong_p(t))
    {
        flint_printf("Exception (fmpz_mat_fread). "
               "Number of rows does not fit into a slong.\n");
        abort();
    }
    r = flint_mpz_get_si(t);

    if(fscanf(file, " ") != 0) {
      abort();
    }

    /* second number in file should be column dimension */
    byte_count = mpz_inp_str(t, file, 10);
    if (byte_count == 0)
    {
        mpz_clear(t);
        return 0;
    }
    
    if (!mpz_fits_slong_p(t))
    {
        flint_printf("Exception (fmpz_mat_fread). "
               "Number of columns does not fit into a slong.\n");
        abort();
    }
    c = flint_mpz_get_si(t);
    mpz_clear(t);
    
    /* if the input is 0 by 0 then set the dimensions to r and c */
    if (mat->r == 0 && mat->c == 0)
    {
        fmpz_mat_clear(mat);
        fmpz_mat_init(mat,r,c);
    }
    else if (mat->r != r || mat->c != c)
    {
        flint_printf("Exception (fmpz_mat_fread). \n"
               "Dimensions are non-zero and do not match input dimensions.\n");
        abort();
    }

    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            if (!fmpz_inp_raw(fmpz_mat_entry(mat, i, j), file))
                return 0;
            if(fscanf(file, " ") != 0) {
              abort();
    }


        }
    }

    /* a return value of 0 means a problem with 
       the file stream a value of 1 means success*/
    return 1;
}



int fmpz_poly_fread_raw(FILE * file, fmpz_poly_t poly)
{
    int r;
    slong i, len;
    mpz_t t;

    mpz_init(t);
    r = mpz_inp_str(t, file, 10);
    if (r == 0)
    {
        mpz_clear(t);
        return 0;
    }
    if (!mpz_fits_slong_p(t))
    {
        flint_printf("Exception (fmpz_poly_fread). Length does not fit into a slong.\n");
        abort();
    }
    len = flint_mpz_get_si(t);
    mpz_clear(t);

    fmpz_poly_fit_length(poly, len);

    for (i = 0; i < len; i++)
    {
        r = fmpz_inp_raw(poly->coeffs + i, file);
        if (r <= 0)
            return r;
    }

    _fmpz_poly_set_length(poly, len);
    _fmpz_poly_normalise(poly);

    return 1;
}


int _fmpz_poly_fprint_raw(FILE * file, const fmpz * vec, slong len)
{
    int r;
    slong i;

    r = flint_fprintf(file, "%li", len);
    if ((len > 0) && (r > 0))
    {
        for (i = 0; (i < len) && (r > 0); i++)
        {
            r = fmpz_out_raw(file, vec + i);
        }
    }

    return r;
}

int fmpz_poly_fprint_raw(FILE * file, const fmpz_poly_t poly)
{
    return _fmpz_poly_fprint_raw(file, poly->coeffs, poly->length);
}



int _fmpz_mod_poly_fprint_raw(FILE * file, const fmpz *poly, slong len, 
                          const fmpz_t p)
{
    int r;
    slong i;

    r = flint_fprintf(file, "%wd ", len);
    if (r <= 0)
        return r;

    r = fmpz_out_raw(file, p);
    if (r <= 0)
        return r;

    if (len == 0)
        return r;

    r = flint_fprintf(file, " ");
    if (r <= 0)
        return r;

    for (i = 0; (r > 0) && (i < len); i++)
    {
        r = flint_fprintf(file, " ");
        if (r <= 0)
            return r;
        r = fmpz_out_raw(file, poly + i);
        if (r <= 0)
            return r;
    }

    return r;
}

int fmpz_mod_poly_fprint_raw(FILE * file, const fmpz_mod_poly_t poly)
{
    return _fmpz_mod_poly_fprint_raw(file, poly->coeffs, poly->length,
        &(poly->p));
}

/* copied from fmpz_mod_poly_fread, mostly (but fixed) */
int fmpz_mod_poly_fread_raw(FILE * f, fmpz_mod_poly_t poly)
{
    slong i, length;
    fmpz_t coeff;
    ulong res;

    fmpz_init(coeff);
    if (flint_fscanf(f, "%wd ", &length) != 1) {
        fmpz_clear(coeff);
        return 0;
    }

    fmpz_inp_raw(coeff,f);
    fmpz_mod_poly_clear(poly);
    fmpz_mod_poly_init(poly, coeff);
    fmpz_mod_poly_fit_length(poly, length);

    poly->length = length;
    flint_fscanf(f, " ");
  
    for (i = 0; i < length; i++)
    {
        flint_fscanf(f, " ");
        res = fmpz_inp_raw(coeff, f);

        fmpz_mod_poly_set_coeff_fmpz(poly,i,coeff);

        if (!res)
        {
            poly->length = i;
            fmpz_clear(coeff);
            return 0;
        }
    }

    fmpz_clear(coeff);
    _fmpz_mod_poly_normalise(poly);

    return 1;
}

