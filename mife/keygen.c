#include <errno.h>
#include <getopt.h>
#include <sys/stat.h>
#include <sys/resource.h>

#include <mife/mife.h>
#include <mmap/mmap_dummy.h>
#include <mmap/mmap_clt.h>

#include "cmdline.h"
#include "mbp_glue.h"
#include "mife.h"
#include "parse.h"
#include "util.h"

typedef struct {
  int sec_param, log_db_size;
  aes_randstate_t seed;
  mbp_template template;
} keygen_inputs;

typedef struct {
  location public, private;
} keygen_locations;

void mife_keygen_parse_cmdline(int argc, char **argv, keygen_inputs *const ins, keygen_locations *const outs, bool *use_clt, int *ncores);
bool mife_keygen_print_outputs(const_mmap_vtable mmap, keygen_locations outs, mife_pp_t pp, mife_sk_t sk);
void mife_keygen_cleanup(const_mmap_vtable mmap, keygen_inputs *const ins, keygen_locations *const outs, mife_pp_t pp, mife_sk_t sk);

int main(int argc, char **argv) {
  keygen_inputs ins;
  keygen_locations outs;
  mife_pp_t pp;
  mife_sk_t sk;
  mbp_template_stats stats;
  bool success;
  bool use_clt = false;
  int ncores;

  PRINT_TIMERS = 1; /* prints timing/progress info */

  mife_keygen_parse_cmdline(argc, argv, &ins, &outs, &use_clt, &ncores);

  const mmap_vtable *mmap;
  if (use_clt)
    mmap = &clt_vtable;
  else
    mmap = &dummy_vtable;

  if (!mbp_template_to_mife_pp(pp, &ins.template, &stats))
      return -1;

  mife_setup(mmap, pp, sk, ins.log_db_size, ins.sec_param, ncores, ins.seed);
  timer_printf("Finished calling mife_setup. Starting to write outputs...\n");
  start_timer();

  success = mife_keygen_print_outputs(mmap, outs, pp, sk);
  timer_printf("Finished writing outputs");
  print_timer();
  timer_printf("\n");

  timer_printf("Starting cleanup...\n");
  start_timer();
  mife_keygen_cleanup(mmap, &ins, &outs, pp, sk);
  timer_printf("Finished cleanup");
  print_timer();
  timer_printf("\n");

  {
      struct rusage usage;
      (void) getrusage(RUSAGE_SELF, &usage);
      (void) printf("Max memory usage: %ld\n", usage.ru_maxrss);
  }

  return success ? 0 : -1;
}

static void mife_keygen_usage(const int code) {
  /* separate the diagnostic information from the usage information a little bit */
  if(0 != code) printf("\n\n");
  printf(
    "The key generation phase initializes the parameters that are shared across\n"
    "an entire database. Brackets indicate default values for each argument.\n"
    "\n"
    "Common options:\n"
    "  -h, --help         Display this usage information\n"
    "  -r, --private      A directory for private parameters [private]\n"
    "  -u, --public       A directory for public parameters [public]\n"
    "  -C, --clt13        Use CLT13 as the underlying multilinear map\n"
    "  -c, --ncores       Number of cores to use [0]\n"
    "\n"
    "Keygen-specific options:\n"
    "  -s, --secparam     Security parameter [80]\n"
    "  -n, --dbsize       Allow up to 2^n records [80]\n"
    "\n"
    "Files used:\n"
    "  <public>/template.json  R  JSON    a description of the function being\n"
    "                                     encrypted\n"
    "  <public>/mife.pub        W custom  public parameters for evaluating\n"
    "  <private>/mife.priv      W custom  private parameters for encrypting\n"
    "  <private>/seed.bin      R  binary  %d-byte seed for PRNG\n"
    "  /dev/urandom            R  binary  used in case above file is missing\n"
    , AES_SEED_BYTE_SIZE
    );
  exit(code);
}

void mife_keygen_parse_cmdline(int argc, char **argv, keygen_inputs *const ins,
                               keygen_locations *const outs, bool *use_clt, int *ncores)
{
  bool done = false;

  /* set defaults */
  ins->sec_param = 80;
  ins->log_db_size = 80;
  *ncores = 0;
  *outs = (keygen_locations) { { "public", true }, { "private", true } };

  struct option long_opts[] =
    { {"help"     ,       no_argument, NULL, 'h'}
    , {"dbsize"   , required_argument, NULL, 'n'}
    , {"private"  , required_argument, NULL, 'r'}
    , {"secparam" , required_argument, NULL, 's'}
    , {"public"   , required_argument, NULL, 'u'}
    , {"clt"      ,       no_argument, NULL, 'C'}
    , {"ncores"   , required_argument, NULL, 'c'}
    , {NULL, 0, NULL, 0}
    };

  while(!done) {
    int c = getopt_long(argc, argv, "hn:c:Cr:s:u:", long_opts, NULL);
    switch(c) {
      case  -1: done = true; break;
      case   0: break; /* a long option with non-NULL flag; should never happen */
      case '?': mife_keygen_usage(1); break; /* braking is good defensive driving */
      case 'h': mife_keygen_usage(0); break;
    case 'c':
        *ncores = atoi(optarg);
        break;
      case 'n':
        if((ins->log_db_size = atoi(optarg)) < 1) {
          fprintf(stderr, "%s: unparseable database size '%s', should be positive number\n", *argv, optarg);
          mife_keygen_usage(2);
        }
        break;
      case 'r':
        outs->private.path = optarg;
        break;
      case 's':
        if((ins->sec_param = atoi(optarg)) < 1) {
          fprintf(stderr, "%s: unparseable security parameter '%s', should be positive number\n", *argv, optarg);
          mife_keygen_usage(3);
        }
        break;
      case 'u':
        outs->public.path = optarg;
        break;
      case 'C':
        *use_clt = true;
        break;
      default:
        fprintf(stderr, "The impossible happened! getopt returned %d (%c)\n", c, c);
        exit(-1);
        break;
    }
  }

  if(optind < argc) {
    fprintf(stderr, "%s: unexpected non-option argument %s\n", *argv, argv[optind]);
    exit(4);
  }

  /* read template */
  location template_location = location_append(outs->public, "template.json");
  if(NULL == template_location.path) {
    fprintf(stderr, "%s: out of memory when trying to create path to template\n", *argv);
    exit(-1);
  }
  if(!jsmn_parse_mbp_template_location(template_location, &ins->template)) {
    fprintf(stderr, "%s: could not parse '%s' as a\nJSON representation of a matrix branching program template over the field F_2\n", *argv, template_location.path);
    mife_keygen_usage(7);
  }
  location_free(template_location);

  /* read seed */
  check_parse_result(load_seed(outs->private, "keygen", ins->seed), mife_keygen_usage, 8);
}

/* for now, use the mife library's custom format; would be good to upgrade this
 * to something a bit more redundant and human-readable to improve error
 * detection, error reporting, versioning and just generally make this tool a
 * bit more robust
 */
bool mife_keygen_print_outputs(const_mmap_vtable mmap, keygen_locations outs, mife_pp_t pp, mife_sk_t sk) {
  /* the public directory has to exist -- we read template.json out of it! --
   * but the private one might not yet */
  if(!create_directory_if_missing(outs.private.path)) {
    fprintf(stderr, "could not create output directory %s\n", outs.private.path);
    return false;
  }

  location  public_location = location_append(outs.public , "mife.pub" );
  location private_location = location_append(outs.private, "mife.priv");
  if( public_location.path == NULL ||
     private_location.path == NULL) {
    location_free(public_location);
    fprintf(stderr, "out of memory when generating output paths\n");
    return false;
  }

  fwrite_mife_pp(mmap, pp,  public_location.path);
  fwrite_mife_sk(mmap, sk, private_location.path);

  location_free( public_location);
  location_free(private_location);
  return true;
}

void mife_keygen_cleanup(const_mmap_vtable mmap, keygen_inputs *const ins, keygen_locations *const outs, mife_pp_t pp, mife_sk_t sk) {
  aes_randclear(ins->seed);
  mbp_template_stats_free(*(mbp_template_stats *)pp->mbp_params);
  mbp_template_free(ins->template);
  location_free(outs-> public);
  location_free(outs->private);
  mife_clear_pp(pp);
  mife_clear_sk(mmap, sk);
}
