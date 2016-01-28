#ifndef _ORE_H_
#define _ORE_H_

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <gmp.h>
#include <gghlite/gghlite.h>
#include <gghlite/gghlite-defs.h>
#include <mife/mife.h>

/* these are used for get_matrix_bit() */
#define X_TYPE 0
#define Y_TYPE 1
#define NONZERO_VAL 1


#define DEFAULT_SHA_SEED "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"

int CT_SIZE;

typedef enum {
   //!< the MBP produced is a series of 2n (d+3) x (d+3) matrices
  ORE_MBP_NORMAL = 0x08,

  //!< the MBP produced is degree-compressed by reading the input as: x0 (y0 y1) 
  //(x1 x2) (y2 y3) ...
  ORE_MBP_DC = 0x10,

  //!< the MBP produced is matrix-compressed, with dimensions 3 x d. This is not 
  //compatible with ORE_MBP_DC, though.
  ORE_MBP_MC = 0x20,

  ORE_DEFAULT = ORE_MBP_DC,
} ore_flag_t;

ore_flag_t ORE_GLOBAL_FLAGS = ORE_DEFAULT;

typedef struct _ore_params_struct {
  int d; // the base
  int bitstr_len; // length of the plaintexts in d-ary
} ore_params_t;

typedef struct ore_cmdline {
  int lambda;
  int d1;
  int d2;
  int challenge_index;
  int num_messages;
  char *m_file;
  char *seed;

  char *pp_file;
  char *sk_file;
  int num_files;
  char **files;

  int is_cgen;
  int is_plaintexts_gen;
  int is_pp_gen;
  int is_tests_only;
  int is_evaluate_only;

} ore_cmdline_t[1];

void ore_print_help_and_exit();
void ore_parse_cmdline(int argc, char *argv[], ore_cmdline_t params);

/* functions dealing with ORE challenge generation */
void generate_plaintexts(int num_messages, int d, int n, char *seed);
int test_ciphertexts(char *pp_file, char *ct1_file, char *ct2_file);
void test_all_ciphertexts(char *pp_file, int n);
void ore_pp_sk_gen(char *pp_file, char *sk_file, int lambda, int d, int n,
    char *seed);
void ore_challenge_gen(char *m_file, int challenge_index, int lambda,
    int d, int n, char *seed);
void ore_set_best_params(mife_pp_t pp, int lambda, fmpz_t message_space_size,
    ore_params_t *params);
int ore_get_matrix_bit_normal_mbp(int input, int i, int j, int type);

int ore_mbp_param(mife_pp_t pp, int index);
void ore_mbp_kilian(mife_pp_t pp, int *dims);
void ore_mbp_set_matrices(mife_pp_t pp, mife_mat_clr_t met, void *message_untyped);
void ore_mbp_ordering(mife_pp_t pp, int index, int *ip, int *jp);
int ore_mbp_parse(mife_pp_t pp, f2_matrix m);

/* test functions */
void run_tests();
int test_ore(int lambda, int mspace_size, int num_messages, int d,
    int bitstr_len, ore_flag_t flags, int verbose);
void print_random_matrices_with_inverse(int n, char *p_file,
    char *a_file, char *inv_file);


/* benchmarking info for choosing best params */

int MAX_KAPPA_BENCH = 28;
long KAPPA_BENCH[] = {
-1,
-1,
13743796,
39428602,
50011646,
127578831,
149971402,
172363974,
194756546,
217149117,
504072750,
551305371,
598537993,
645770614,
693003236,
740235857,
787468479,
834701100,
881933722,
929166343,
2051532015,
2150882456,
2250232897,
2349583338,
2448933779,
2548284220,
2647634661,
2746985102,
2846335543,
2945685984,
3045036425,
3144386867,
3243737308,
3343087749,
3442438190,
3541788631,
3641139072,
3740489513,
3839839954,
};


#endif /* _ORE_H_ */
