#include <polycap-private.h>
#include <assert.h>

void test_polycap_capil_segment() {
	polycap_error *error = NULL;
	int test = 0;

	//won't work
	test = polycap_capil_segment(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &error);
	assert(test == -1);
	assert(polycap_error_matches(error, POLYCAP_ERROR_INVALID_ARGUMENT));

	//Should work
	polycap_clear_error(&error);
	test = polycap_capil_segment(cap_coord0, cap_coord1, cap_rad0, cap_rad1, photon_coord, photon_dir, surface_norm, alfa, &error);
	assert(test == 0);
	
}

void test_polycap_refl() {


}

void test_polycap_capil_reflect() {


}

void test_polycap_capil_trace() {


}

int main(int argc, char *argv[]) {

	test_polycap_capil_segment();
	test_polycap_refl();
	test_polycap_capil_reflect();
	test_polycap_capil_trace();

	return 0;
}
