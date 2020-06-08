#include "polycap-private.h"
#include <string.h>
#include <stdio.h>

int main(int argc, char **argv){

	char *header_version = NULL;
	asprintf(&header_version, "%d.%d", POLYCAP_VERSION_MAJOR, POLYCAP_VERSION_MINOR);
	if(strcmp(header_version, PACKAGE_VERSION) == 0){
		return 0;
	}
	return 1;
}
