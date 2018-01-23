polycap_profile* polycap_profile_new_from_file(
	const char *single_cap_profile_file,
	const char *central_axis_file,
	const char *external_shape_file)
{
	struct cap_prof_arrays *shape_arr;
	struct inp_file *cap;
	struct cap_profile *profile;

	//duplicate file names to structure as accepted by read_cap_profile
	cap->prf = strdup(single_cap_profile_file);
	cap->axs = strdup(central_axis_file);
	cap->ext = strdup(external_shape_file);

	//create profile struct for output as made by read_cap_profile
	profile = malloc(sizeof(struct cap_profile));

	//read files
	read_cap_profile(cap, profile);

	//transfer shape information to structure with identical output as polycap_profile_new()
	shape_arr = malloc(sizeof(struct cap_prof_arrays)*(profile->nmax+1));
	for(i=0;i<=profile->nmax;i++){
		shape_arr[i].zarr = profile->arr[i].zarr;
		shape_arr[i].profil = profile->arr[i].profil;
		shape_arr[i].d_arr = profile->arr[i].d_arr;
	}

	free(cap);
	free(profile->arr);
	free(profile);

	return shape_arr;

}
