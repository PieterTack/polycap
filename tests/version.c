#include "polycap-private.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv){

	char *header_version = NULL;

#ifdef _WIN32
	int bytes_needed = _scprintf("%d.%d", POLYCAP_VERSION_MAJOR, POLYCAP_VERSION_MINOR);
	if (bytes_needed < 0)
		return 1;
	header_version = malloc((bytes_needed + 1) * sizeof(char));
	if (_snprintf(header_version, bytes_needed + 1, "%d.%d", POLYCAP_VERSION_MAJOR, POLYCAP_VERSION_MINOR) < 0) {
		return 1;
	}
#else
	if (asprintf(&header_version, "%d.%d", POLYCAP_VERSION_MAJOR, POLYCAP_VERSION_MINOR) < 0) {
		fprintf(stderr, "vasprintf error\n");
		return 1;
	}
#endif

	if(strcmp(header_version, PACKAGE_VERSION) == 0){
		return 0;
	}
	return 1;
}
