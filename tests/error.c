#include "polycap-error.h"

#include <assert.h>
#include <string.h>

static void test_literal(void) {
	polycap_error *error = NULL;

	polycap_set_error_literal(&error, POLYCAP_ERROR_MEMORY, "%s %d %x");

	assert(polycap_error_matches(error, POLYCAP_ERROR_MEMORY) == true);
	assert(strcmp(error->message, "%s %d %x") == 0);
	polycap_error_free(error);
}

static void test_copy(void) {
	polycap_error *error = NULL, *copy = NULL;

	polycap_set_error_literal(&error, POLYCAP_ERROR_MEMORY, "%s %d %x");
	copy = polycap_error_copy(error);

	assert(polycap_error_matches(copy, POLYCAP_ERROR_MEMORY) == true);
	assert(strcmp(copy->message, "%s %d %x") == 0);
	polycap_error_free(error);
	polycap_error_free(copy);
}

int main(int argc, char *argv[]) {

	test_literal();
	test_copy();

	return 0;
}
