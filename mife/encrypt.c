#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <unistd.h>
#include <string.h>

#include <mife/mife.h>
#include <mmap/mmap_gghlite.h>
#include <mmap/mmap_clt.h>

#include "cmdline.h"
#include "mbp_types.h"
#include "mbp_glue.h"
#include "mife.h"
#include "parse.h"
#include "util.h"

#define UID_TRIES 100
#define INT_STR_LEN (3 * sizeof(int))

typedef struct {
    mbp_plaintext pt;
    location record_location;
    fmpz_t partition;
    mife_sk_t sk;
    mife_pp_t pp;
    aes_randstate_t seed;
} encrypt_inputs;


void mife_encrypt_parse_cmdline(int argc, char **argv, encrypt_inputs *const ins, bool *use_clt);
bool mife_encrypt_print_output(const_mmap_vtable mmap, mife_pp_t pp, int global_index, mmap_enc_mat_t ct, location record_location);
void mife_encrypt_cleanup(const_mmap_vtable mmap, encrypt_inputs *const ins);

int main(int argc, char **argv) {
    encrypt_inputs ins;
    mife_mat_clr_t clr;
    int ***partitions;
    bool success = true;
    bool use_clt = false;

    PRINT_TIMERS = 1; /* prints timing/progress info */

    mife_encrypt_parse_cmdline(argc, argv, &ins, &use_clt);

    const mmap_vtable *mmap;
    if (use_clt) {
        mmap = &clt_vtable;
    } else {
        mmap = &gghlite_vtable;
    }

    mife_encrypt_setup(ins.pp, ins.partition, &ins.pt, clr, &partitions);

    /**
     * determine the total # of encodings we need to make for this ciphertext, for
     * benchmarking and progress bar purposes
     */
    int encoding_count = 0;
    const mbp_template *const template = ((mbp_template_stats *)ins.pp->mbp_params)->template;
    for(unsigned int i = 0; i < template->steps_len; i++) {
        const f2_matrix *const m = template->steps[i].matrix;
        encoding_count += m->num_rows * m->num_cols;
    }

    set_NUM_ENC(encoding_count);

    // now, perform the actual encryption

    reset_T();
    for(unsigned int i = 0; i < template->steps_len; i++) {
        mmap_enc_mat_t ct;
        mife_encrypt_single(mmap, ins.pp, ins.sk, ins.seed, i, clr, partitions, ct);
        success &= mife_encrypt_print_output(mmap, ins.pp, i, ct, ins.record_location);
        mmap_enc_mat_clear(mmap, ct);
    }
    timer_printf("\n");
    mife_encrypt_clear(ins.pp, clr, partitions);

    mife_encrypt_cleanup(mmap, &ins);

    {
        struct rusage usage;
        (void) getrusage(RUSAGE_SELF, &usage);
        (void) printf("Max memory usage: %ld\n", usage.ru_maxrss);
    }

    return success ? 0 : -1;
}

static void mife_encrypt_usage(const int code) {
    /* separate the diagnostic information from the usage information a little bit */
    if(0 != code) printf("\n\n");
    printf(
        "USAGE: encrypt [OPTIONS] PLAINTEXT\n"
        "The encryption operation hides the information in a single plaintext. The\n"
        "plaintext is should be represented as a JSON array containing strings naming\n"
        "symbols from the template available in the public parameters directory.\n"
        "Brackets indicate default values for each argument.\n"
        "\n"
        "Common options:\n"
        "  -h, --help               Display this usage information\n"
        "  -r, --private            A directory for private parameters [private]\n"
        "  -u, --public             A directory for public parameters [public]\n"
        "  -d, --db, --database     A directory to store encrypted values in [database]\n"
        "  -C, --clt13              Use CLT13 as the underlying multilinear map\n"
        "  -s, --sequential         Disable parallelism\n"
        "\n"
        "Encryption-specific options:\n"
        "  -i, --uid                A string that uniquely identifies this record and\n"
        "                           will be publically visible [insecurely random,\n"
        "                           reported on stdout]\n"
        "  -a, --partition          A number that uniquely identifies this record and\n"
        "                           is at most the database size specified during\n"
        "                           the keygen phase; not published by this tool (but\n"
        "                           you may choose to publish it without threatening\n"
        "                           the security of the data) [securely random]\n"
        "\n"
        "Files used:\n"
        "  <database>/<uid>/*/*.bin   W binary  the encrypted record\n"
        "  <public>/template.json    R  JSON    a description of the function being\n"
        "                                       encrypted\n"
        "  <public>/mife.pub         R  custom  public parameters for evaluating\n"
        "  <private>/mife.priv       R  custom  private parameters for encrypting\n"
        "  <private>/seed.bin        R  binary  %d-byte seed for PRNG\n"
        "  /dev/urandom              R  binary  used in case above file is missing\n"
        , AES_SEED_BYTE_SIZE
        );
    exit(code);
}

/* reads (num_bits/8 + 1)*8 bytes into n, mod 2^num_bits */
static bool fmpz_read_bits(fmpz_t n, FILE *file, int num_bits) {
    fmpz_zero(n);
    while(num_bits >= 8) {
        int c = fgetc(file);
        if(EOF == c) return false;
        fmpz_mul_ui(n, n, 256);
        fmpz_add_ui(n, n, c);
        num_bits -= 8;
    }
    if(num_bits >= 0) {
        int c = fgetc(file);
        if(EOF == c) return false;
        c = c % (1 << num_bits);
        fmpz_mul_ui(n, n, 1 << num_bits);
        fmpz_add_ui(n, n, c);
    }
    return true;
}

void mife_encrypt_parse_cmdline(int argc, char **argv, encrypt_inputs *const ins, bool *use_clt) {
    unsigned int i, j;
    bool done = false;
    char *uid = NULL;
    bool have_partition = false;

    /* set defaults */
    location public_location = {   "public", true },
            private_location = {  "private", true },
           database_location = { "database", true };
    fmpz_init(ins->partition);

    struct option long_opts[] =
        { {"db"       , required_argument, NULL, 'd'}
        , {"database" , required_argument, NULL, 'd'}
        , {"help"     ,       no_argument, NULL, 'h'}
        , {"uid"      , required_argument, NULL, 'i'}
        , {"partition", required_argument, NULL, 'a'}
        , {"private"  , required_argument, NULL, 'r'}
        , {"public"   , required_argument, NULL, 'u'}
        , {"clt"      ,       no_argument, NULL, 'C'}
        , {"sequential",      no_argument, NULL, 's'}
        };

    g_parallel = 1;

    while(!done) {
        int c = getopt_long(argc, argv, "a:d:hi:Cr:su:", long_opts, NULL);
        switch(c) {
            case  -1: done = true; break;
            case   0: break; /* a long option with non-NULL flag; should never happen */
            case '?': mife_encrypt_usage(1); break; /* braking is good defensive driving */
            case 'a':
                if(have_partition) fmpz_clear(ins->partition);
                if(0 != fmpz_set_str(ins->partition, optarg, 10)) {
                    fprintf(stderr, "%s: could not read %s as a base-10 number\n", *argv, optarg);
                    mife_encrypt_usage(2);
                }
                have_partition = true;
                break;
            case 'd':
                location_free(database_location);
                database_location = (location) { .path = optarg, .stack_allocated = true };
                break;
            case 'h': mife_encrypt_usage(0); break;
            case 'i':
                uid = optarg;
                break;
            case 'C':
                *use_clt = true;
                break;
            case 'r':
                location_free(private_location);
                private_location = (location) { optarg, true };
                break;
            case 's':
                g_parallel = 0;
                break;
            case 'u':
                location_free(public_location);
                public_location = (location) { optarg, true };
                break;
            default:
                fprintf(stderr, "The impossible happened! getopt returned %d (%c)\n", c, c);
                exit(-1);
                break;
        }
    }

    /* read the plaintext */
    if(optind != argc-1) {
        fprintf(stderr, "%s: specify exactly one plaintext (found %d)\n", *argv, argc-optind);
        mife_encrypt_usage(2);
    }
    if(!jsmn_parse_mbp_plaintext_string(argv[optind], &ins->pt)) {
        fprintf(stderr, "%s: could not parse plaintext as JSON array of strings\n", *argv);
        mife_encrypt_usage(3);
    }

    /* read the template */
    location template_location = location_append(public_location, "template.json");
    mbp_template *template = NULL;
    mbp_template_stats *stats = NULL;
    if(template_location.path == NULL || ALLOC_FAILS(template, 1) || ALLOC_FAILS(stats, 1)) {
        fprintf(stderr, "%s: out of memory while loading template\n", *argv);
        exit(-1);
    }
    if(!jsmn_parse_mbp_template_location(template_location, template)) {
        fprintf(stderr, "%s: could not parse '%s' as a\nJSON representation of a matrix branching program template over the field F_2\n", *argv, template_location.path);
        mife_encrypt_usage(4);
    }
    if(!mbp_template_to_mife_pp(ins->pp, template, stats)) {
        fprintf(stderr, "%s: internal error while computing statistics for template\n", *argv);
        exit(-1);
    }
    location_free(template_location);

    /* read the public parameters */
    location pp_location = location_append(public_location, "mife.pub");
    if(pp_location.path == NULL) {
        fprintf(stderr, "%s: out of memory while loading public parameters\n", *argv);
        exit(-1);
    }

    // If we're going to use the mmap in this function, we should know which
    // one to use.
    const mmap_vtable *mmap;
    if (*use_clt) {
        mmap = &clt_vtable;
    } else {
        mmap = &gghlite_vtable;
    }

    /* TODO: some error-checking would be nice here */
    fread_mife_pp(mmap, ins->pp, pp_location.path);
    location_free(pp_location);

    /* read the secret key */
    location sk_location = location_append(private_location, "mife.priv");
    if(sk_location.path == NULL) {
        fprintf(stderr, "%s: out of memory while loading private key\n", *argv);
        exit(-1);
    }
    /* TODO: error-checking */
    fread_mife_sk(mmap, ins->sk, sk_location.path);
    location_free(sk_location);

    /* initialize record_path, ensuring uid is initialized as a side effect */
    if(NULL != uid) {
        ins->record_location = location_append(database_location, uid);
        if(NULL == ins->record_location.path) {
            fprintf(stderr, "out of memory while building path %s/%s\n", database_location.path, uid);
            exit(-1);
        }
    } else {
        const size_t offset = strlen(database_location.path) + 1;
        /* ought to be able to fit 5 bits per base-62 character */
        const size_t record_location_size = offset + ins->pp->L / 5 + 1;
        /* it is not important that the uid be cryptographically random, so just
         * use /dev/urandom */
        FILE *urandom = fopen("/dev/urandom", "rb");
        struct stat file_stat;
        fmpz_t uid_num;
        if(NULL == urandom) {
            fprintf(stderr, "no uid specified and could not open /dev/urandom to choose one\n");
            exit(-1);
        }
        if(ALLOC_FAILS(ins->record_location.path, record_location_size)) {
            fclose(urandom);
            fprintf(stderr, "out of memory while finding unused uid\n");
            exit(-1);
        }

        memcpy(ins->record_location.path, database_location.path, offset-1);
        ins->record_location.path[offset-1] = '/';
        ins->record_location.stack_allocated = false;
        uid = ins->record_location.path+offset;

        for(i = 0; i < UID_TRIES; i++) {
            fmpz_init(uid_num);
            /* from the birthday bound: if we draw 2^L elements from a pool of
             * 2^(2L+8) possible elements, we have a roughly 0.2% chance of
             * failing; that should be sufficiently low for our purposes */
            fmpz_read_bits(uid_num, urandom, 2*ins->pp->L + 8);
            fmpz_get_str(uid, 62, uid_num);
            /* probably not totally foolproof, but is a decent quick check that
             * a file doesn't exist
             */
            if(stat(ins->record_location.path, &file_stat)) break;
            fmpz_clear(uid_num);
        }

        fclose(urandom);

        if(UID_TRIES == i) {
            fprintf(stderr, "tried %d random uids, but they were all in use\n", UID_TRIES);
            exit(-1);
        }

        printf("%s\n", uid);
    }

    /* initialize the random seed */
    const char function_name[] = "encrypt";
    const size_t context_size = strlen(uid) + sizeof(function_name), context_len = context_size-1;
    char *context;
    if(ALLOC_FAILS(context, context_size)) {
        fprintf(stderr, "%s: out of memory when generating context for RNG seed\n", *argv);
        exit(-1);
    }
    const unsigned int tmp = snprintf(context, context_size, "%s%s", function_name, uid);
    if(context_len != tmp) {
        fprintf(stderr, "The impossible happened: expecting '%s%s' to have length %lu, but snprintf reports it has length %d.\n",
            function_name, uid, context_len, tmp);
        exit(-1);
    }
    check_parse_result(load_seed(private_location, context, ins->seed), mife_encrypt_usage, 5);

    /* initialize partition if it wasn't specified on the command line */
    if(!have_partition)
        /* TODO: this cast -- from int to mp_bitcnt_t -- is probably fine...
         * right??? the FLINT docs are surprisingly quiet about mp_bitcnt_t */
        fmpz_randbits_aes(ins->partition, ins->seed, ins->pp->L);
    /* probably unnecessary, but might as well keep this up to date */
    have_partition = true;

    /* from here on out, it's all sanity-checks */

    /* check that the partition is in range */
    if(fmpz_sizeinbase(ins->partition, 2) > (size_t)ins->pp->L) {
        fprintf(stderr, "partition uses %zu bits, but the current key supports only up to %d bits\n", fmpz_sizeinbase(ins->partition, 2), ins->pp->L);
        mife_encrypt_usage(6);
    }

    /* check that the template and plaintext have the same length */
    if(template->steps_len != ins->pt.symbols_len) {
        fprintf(stderr, "the number of symbols in the plaintext (%d)\ndoes not match the number of steps in the template (%d)\n",
                ins->pt.symbols_len, template->steps_len);
        mife_encrypt_usage(7);
    }

    /* check that the template and plaintext match up appropriately */
    bool match_everywhere = true;
    for(i = 0; i < template->steps_len; i++) {
        bool match_here = false;
        for(j = 0; j < template->steps[i].symbols_len; j++)
            match_here |= !strcmp(template->steps[i].symbols[j], ins->pt.symbols[i]);
        if(!match_here) {
            fprintf(stderr, "the plaintext symbol %s at index %d is unknown\n", ins->pt.symbols[i], i);
            fprintf(stderr, "\t(known symbols: ");
            for(j = 0; j < template->steps[i].symbols_len-1; j++)
                fprintf(stderr, "%s, ", template->steps[i].symbols[j]);
            fprintf(stderr, "%s)\n", template->steps[i].symbols[j]);
        }
        match_everywhere &= match_here;
    }
    if(!match_everywhere) mife_encrypt_usage(8);

    /* TODO: check that `ins->pp` and `stats` match up */
    /* TODO: check that the secret key is appropriately dimensioned */

    /* these do nothing for now, and are only here as a defensive measure
     * against future refactorings */
    location_free(public_location);
    location_free(private_location);
}

static FILE *mife_encrypt_fopen_bin(const location record_location, const char *const position, const int local_index) {
    int tmp;
    char *dir, *path;
    const size_t  dir_len = strlen(record_location.path) + 1 + strlen(position), dir_size = dir_len + 1;
    const size_t path_len = dir_len + 1 + INT_STR_LEN + 4, path_size = path_len + 1;
    FILE *result = NULL;

    if(ALLOC_FAILS( dir,  dir_size)) goto done;
    if(ALLOC_FAILS(path, path_size)) goto free_dir;
    tmp = snprintf( dir,  dir_size, "%s/%s"    , record_location.path, position);
    assert(tmp < dir_size);
    tmp = snprintf(path, path_size, "%s/%d.bin", dir, local_index);
    assert(tmp < path_size);
    if(!create_directory_if_missing(dir)) goto free_path;
    result = fopen(path, "wb");

free_path:
    free(path);
free_dir:
    free(dir);
done:
    return result;
}

bool mife_encrypt_print_output(const_mmap_vtable mmap, mife_pp_t pp, int global_index, mmap_enc_mat_t ct, location record_location) {
    const mbp_template_stats *const stats    = pp->mbp_params;
    const mbp_template       *const template = stats->template;

    const int local_index = stats->local_index[global_index];
    const char *const position = template->steps[global_index].position;
    FILE *dest = mife_encrypt_fopen_bin(record_location, position, local_index);
    if(NULL == dest) {
        fprintf(stderr, "could not write %s/%d.bin\n", position, local_index);
        return false;
    }
    fwrite_mmap_enc_mat(mmap, ct, dest);
    fclose(dest);
    return true;
}

void mife_encrypt_cleanup(const_mmap_vtable mmap, encrypt_inputs *const ins) {
    mbp_plaintext_free(ins->pt);
    location_free(ins->record_location);
    fmpz_clear(ins->partition);
    mife_clear_sk(mmap, ins->sk);
    mbp_template_stats *stats = ins->pp->mbp_params;
    mbp_template *templ = (mbp_template *)stats->template;
    mbp_template_stats_free(*stats); free(stats);
    mbp_template_free(*templ); free(templ);
    mife_clear_pp_read(mmap, ins->pp);
    aes_randclear(ins->seed);
}
