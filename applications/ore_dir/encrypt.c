#include "matrix.h"
#include "mife.h"
#include "parse.h"
#include "util.h"

typedef struct {
	plaintext pt;
	fmpz_t uid;
	mife_sk_t sk;
	mife_pp_t pp;
	aes_randstate_t seed;
} encrypt_inputs;

void parse_cmdline(int argc, char **argv, encrypt_inputs *const ins, location *const outs);
void print_outputs(location outs, mife_pp_t pp, mife_ciphertext_t ct);
void cleanup(encrypt_inputs *const ins, location *const outs);

int main(int argc, char **argv) {
	encrypt_inputs ins;
	location outs;
	mife_ciphertext_t ct;

	parse_cmdline(argc, argv, &ins, &outs);
	/* TODO: keeping the whole ciphertext in memory is probably infeasible */
	mife_encrypt(ct, &ins.pt, ins.pp, ins.sk, ins.seed);
	print_outputs(outs, ins.pp, ct);

	cleanup(&ins, &outs);
}

void parse_cmdline(int argc, char **argv, encrypt_inputs *const ins, location *const outs) {
	/* TODO */
}

void print_outputs(location outs, mife_pp_t pp, mife_ciphertext_t ct) {
	/* TODO */
}

void cleanup(encrypt_inputs *const ins, location *const outs) {
	/* TODO */
}
