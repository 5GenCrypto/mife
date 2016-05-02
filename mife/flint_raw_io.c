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
