#include "polycap-private.h"

//===========================================
// get a new polycap_source by providing all its properties 
polycap_source* polycap_source_new(double d_source, double src_x, double src_y, double src_sigx, double src_sigy, double src_shiftx, double src_shifty)
{
	polycap_source *source;

	source = malloc(sizeof(polycap_source));
	if(source == NULL){
		printf("Could not allocate source memory.\n");
		exit(1);
	}

	source->d_source = d_source;
	source->src_x = src_x;
	source->src_y = src_y;
	source->src_sigx = src_sigx;
	source->src_sigy = src_sigy;
	source->src_shiftx = src_shiftx;
	source->src_shifty = src_shifty;

	return source;
}
//===========================================
// free a polycap_source struct
void polycap_source_free(polycap_source *source)
{
	free(source);

	return;
}

