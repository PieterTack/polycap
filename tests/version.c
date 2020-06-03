#include "polycap-private.h"
#include <string.h>
#include <stdio.h>

int main(int argc, char **argv){

	//gchar **header_version = g_strdup_printf("%d.%d", POLYCAP_VERSION_MAJOR, POLYCAP_VERSION_MINOR);
	char **header_version = vasprintf("%d.%d", POLYCAP_VERSION_MAJOR, POLYCAP_VERSION_MINOR);
	if(strcmp(header_version, PACKAGE_VERSION) == 0){
		return 0;
	}
	return 1;
}
