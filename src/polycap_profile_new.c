polycap_profile* polycap_profile_new(
	enum polycap_profile_type type,
	double length,
	double rad_ext[2],
	double rad_int[2],
	double focal_dist[2])
{
	struct cap_prof_arrays *shape_arr;

	shape_arr = def_cap_profile(enum polycap_profile_type type, double length, double rad_ext[2], double rad_int[2], double focal_dist[2]);

	return shape_arr;

}

