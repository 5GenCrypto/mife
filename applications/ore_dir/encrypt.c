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

void parse_cmdline(int argc, char **argv, encrypt_inputs *const ins, location *const database_location);
void print_outputs(location database, mife_pp_t pp, mife_ciphertext_t ct);
void cleanup(encrypt_inputs *const ins, location *const database_location);

int main(int argc, char **argv) {
	encrypt_inputs ins;
	location database_location;
	mife_ciphertext_t ct;

	parse_cmdline(argc, argv, &ins, &database_location);
	/* TODO: keeping the whole ciphertext in memory is probably infeasible */
	mife_encrypt(ct, &ins.pt, ins.pp, ins.sk, ins.seed);
	print_outputs(database_location, ins.pp, ct);

	cleanup(&ins, &database_location);
}

void parse_cmdline(int argc, char **argv, encrypt_inputs *const ins, location *const database_location) {
	/* TODO */
}

void print_outputs(location database_location, mife_pp_t pp, mife_ciphertext_t ct) {
	/* TODO */
}

void cleanup(encrypt_inputs *const ins, location *const database_location) {
	/* TODO */
}
