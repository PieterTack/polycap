/*
 * Copyright (C) 2018 Pieter Tack, Tom Schoonjans and Laszlo Vincze
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * */

#include "polycap.h"

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
