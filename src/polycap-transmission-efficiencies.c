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

#include "polycap-private.h"
#include "polycap-aux.h"
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <hdf5.h>
#include <math.h>

/* Ideally, we would need to do something similar on Windows, instead of just disabling the use of a mutex...
 * Still, seems unlikely to cause problems anytime soon...
 */
#ifndef _WIN32
#include <pthread.h>

static pthread_mutex_t tables_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

/* error handling borrowed from h5py */

struct _minor_table_entry {
	H5E_minor_t min_num;
	enum polycap_error_code code;
};

/* "Fudge" table to accommodate annoying inconsistencies in HDF5's use
 of the minor error codes.  If a (major, minor) entry appears here,
 it will override any entry in the minor error table.*/
struct _exact_table_entry {
	H5E_major_t maj_num;
	H5E_minor_t min_num;
	enum polycap_error_code code;
};

/* If the minor error code matches an entry
 in this dict, the generated exception will be used.
*/
static struct _minor_table_entry *_minor_table = NULL;
static unsigned int _minor_table_len = 0;
static unsigned int _minor_table_cap = 0;
static struct _exact_table_entry *_exact_table = NULL;
static unsigned int _exact_table_len = 0;
static unsigned int _exact_table_cap = 0;

static void append_minor_table_entry(H5E_minor_t h5_min_num, enum polycap_error_code pc_code) {
	if (_minor_table_len == _minor_table_cap) {
		_minor_table_cap *= 2;
		_minor_table = realloc(_minor_table, sizeof(struct _minor_table_entry) * _minor_table_cap);
	}
	_minor_table[_minor_table_len].min_num = h5_min_num;
	_minor_table[_minor_table_len].code = pc_code;

	_minor_table_len++;
}

static void append_exact_table_entry(H5E_major_t h5_maj_num, H5E_minor_t h5_min_num, enum polycap_error_code pc_code) {
	if (_exact_table_len == _exact_table_cap) {
		_exact_table_cap *= 2;
		_exact_table = realloc(_exact_table, sizeof(struct _exact_table_entry) * _exact_table_cap);
	}
	_exact_table[_exact_table_len].min_num = h5_min_num;
	_exact_table[_exact_table_len].maj_num = h5_maj_num;
	_exact_table[_exact_table_len].code = pc_code;

	_exact_table_len++;
}

static void tables_init() {
#ifndef _WIN32
	pthread_mutex_lock(&tables_mutex);
#endif

	// silence HDF5
	H5Eset_auto(H5E_DEFAULT, NULL, NULL);

	if (_minor_table_len == 0) {
		_minor_table_len = 0;
		_minor_table_cap = 1;

		_minor_table = malloc(sizeof(struct _minor_table_entry));

		// populate table
		append_minor_table_entry(H5E_SEEKERROR,      POLYCAP_ERROR_IO);    // Seek failed
		append_minor_table_entry(H5E_READERROR,      POLYCAP_ERROR_IO);    // Read failed
		append_minor_table_entry(H5E_WRITEERROR,     POLYCAP_ERROR_IO);    // Write failed
		append_minor_table_entry(H5E_CLOSEERROR,     POLYCAP_ERROR_IO);    // Close failed
		append_minor_table_entry(H5E_OVERFLOW,       POLYCAP_ERROR_IO);    // Address overflowed
		append_minor_table_entry(H5E_FCNTL,          POLYCAP_ERROR_IO);    // File control (fcntl) failed
		append_minor_table_entry(H5E_FILEEXISTS,     POLYCAP_ERROR_IO);    // File already exists
		append_minor_table_entry(H5E_FILEOPEN,       POLYCAP_ERROR_IO);    // File already open
		append_minor_table_entry(H5E_CANTCREATE,     POLYCAP_ERROR_IO);    // Unable to create file
		append_minor_table_entry(H5E_CANTOPENFILE,   POLYCAP_ERROR_IO);    // Unable to open file
		append_minor_table_entry(H5E_CANTCLOSEFILE,  POLYCAP_ERROR_IO);    // Unable to close file
		append_minor_table_entry(H5E_NOTHDF5,        POLYCAP_ERROR_IO);    // Not an HDF5 file
		append_minor_table_entry(H5E_BADFILE,        POLYCAP_ERROR_INVALID_ARGUMENT); // Bad file ID accessed
		append_minor_table_entry(H5E_TRUNCATED,      POLYCAP_ERROR_IO);    // File has been truncated
		append_minor_table_entry(H5E_MOUNT,          POLYCAP_ERROR_IO);    // File mount error
		append_minor_table_entry(H5E_NOFILTER,       POLYCAP_ERROR_IO);    // Requested filter is not available
		append_minor_table_entry(H5E_CALLBACK,       POLYCAP_ERROR_IO);    // Callback failed
		append_minor_table_entry(H5E_CANAPPLY,       POLYCAP_ERROR_IO);    // Error from filter 'can apply' callback
		append_minor_table_entry(H5E_SETLOCAL,       POLYCAP_ERROR_IO);    // Error from filter 'set local' callback
		append_minor_table_entry(H5E_NOENCODER,      POLYCAP_ERROR_IO);    // Filter present but encoding disabled
		append_minor_table_entry(H5E_BADGROUP,       POLYCAP_ERROR_INVALID_ARGUMENT);  // Unable to find ID group information
		append_minor_table_entry(H5E_BADSELECT,      POLYCAP_ERROR_INVALID_ARGUMENT);  // Invalid selection (hyperslabs)
		append_minor_table_entry(H5E_UNINITIALIZED,  POLYCAP_ERROR_INVALID_ARGUMENT);  // Information is uninitialized
		append_minor_table_entry(H5E_UNSUPPORTED,    POLYCAP_ERROR_UNSUPPORTED);    // Feature is unsupported
		append_minor_table_entry(H5E_NOTFOUND,       POLYCAP_ERROR_INVALID_ARGUMENT);    // Object not found
		append_minor_table_entry(H5E_CANTINSERT,     POLYCAP_ERROR_INVALID_ARGUMENT);   // Unable to insert object
		append_minor_table_entry(H5E_BADTYPE,        POLYCAP_ERROR_TYPE);   // Inappropriate type
		append_minor_table_entry(H5E_BADRANGE,       POLYCAP_ERROR_INVALID_ARGUMENT);  // Out of range
		append_minor_table_entry(H5E_BADVALUE,       POLYCAP_ERROR_INVALID_ARGUMENT);  // Bad value
		append_minor_table_entry(H5E_EXISTS,         POLYCAP_ERROR_INVALID_ARGUMENT);  // Object already exists
		append_minor_table_entry(H5E_ALREADYEXISTS,  POLYCAP_ERROR_INVALID_ARGUMENT);  // Object already exists, part II
		append_minor_table_entry(H5E_CANTCONVERT,    POLYCAP_ERROR_TYPE);   // Can't convert datatypes
		append_minor_table_entry(H5E_CANTDELETE,     POLYCAP_ERROR_INVALID_ARGUMENT);    // Can't delete message
		append_minor_table_entry(H5E_CANTOPENOBJ,    POLYCAP_ERROR_INVALID_ARGUMENT);
		append_minor_table_entry(H5E_CANTMOVE,       POLYCAP_ERROR_INVALID_ARGUMENT);  // Can't move a link
	}

	if (_exact_table_len == 0) {
		_exact_table_len = 0;
		_exact_table_cap = 1;

		_exact_table = malloc(sizeof(struct _exact_table_entry));

		// populate table
		append_exact_table_entry(H5E_CACHE, H5E_BADVALUE,      POLYCAP_ERROR_IO);  // obj create w/o write intent 1.8
		append_exact_table_entry(H5E_RESOURCE, H5E_CANTINIT,   POLYCAP_ERROR_IO);  // obj create w/o write intent 1.6
		append_exact_table_entry(H5E_INTERNAL, H5E_SYSERRSTR,  POLYCAP_ERROR_IO);  // e.g. wrong file permissions
		append_exact_table_entry(H5E_DATATYPE, H5E_CANTINIT,   POLYCAP_ERROR_TYPE);  // No conversion path
		append_exact_table_entry(H5E_DATASET, H5E_CANTINIT,    POLYCAP_ERROR_INVALID_ARGUMENT);  // bad param for dataset setup
		append_exact_table_entry(H5E_ARGS, H5E_CANTINIT,       POLYCAP_ERROR_TYPE);  // Illegal operation on object
		append_exact_table_entry(H5E_ARGS, H5E_BADTYPE,        POLYCAP_ERROR_INVALID_ARGUMENT); // Invalid location in file
		append_exact_table_entry(H5E_REFERENCE, H5E_CANTINIT,  POLYCAP_ERROR_INVALID_ARGUMENT); // Dereferencing invalid ref

#if H5_VERSION_GE(1, 12, 0)
		append_exact_table_entry(H5E_DATASET, H5E_CANTCREATE,  POLYCAP_ERROR_INVALID_ARGUMENT);  // bad param for dataset setup
#endif

#if H5_VERSION_LE(1, 13, 0)
		append_minor_table_entry(H5E_BADATOM,                  POLYCAP_ERROR_INVALID_ARGUMENT);  // Unable to find atom information (already closed?)
		append_exact_table_entry(H5E_SYM, H5E_CANTINIT,        POLYCAP_ERROR_INVALID_ARGUMENT);  // # Object already exists/1.8
#else
		append_minor_table_entry(H5E_BADID,                    POLYCAP_ERROR_INVALID_ARGUMENT);  // Unable to find ID information
		append_exact_table_entry(H5E_SYM, H5E_CANTCREATE,      POLYCAP_ERROR_INVALID_ARGUMENT);  // # Object already exists
#endif
	}

#ifndef _WIN32
	pthread_mutex_unlock(&tables_mutex);
#endif
}

struct err_data_t {
    H5E_error_t err;
    int n;
};

static herr_t walk_cb(unsigned n, const H5E_error2_t *desc, void *e) {

    struct err_data_t *ee = e;

    ee->err.maj_num = desc->maj_num;
    ee->err.min_num = desc->min_num;
    ee->err.desc = desc->desc;
    ee->n = n;
    return 0;
}

static void set_exception(polycap_error **error) {
	if (error == NULL)
		return;

	struct err_data_t err;
	char *desc = NULL;
	const char *desc_bottom = NULL;

	if (H5Ewalk(H5E_DEFAULT, H5E_WALK_UPWARD, walk_cb, &err) < 0) {
		polycap_set_error_literal(error, POLYCAP_ERROR_RUNTIME, "Failed to walk HDF5 error stack");
		return;
	}

	if (err.n < 0) {
		polycap_set_error_literal(error, POLYCAP_ERROR_RUNTIME, "No HDF5 exception information found");
        	return;
	}

    	enum polycap_error_code code = POLYCAP_ERROR_RUNTIME;
	int i;
	for (i = 0 ; i < _minor_table_len ; i++) {
		if (_minor_table[i].min_num == err.err.min_num) {
			code = _minor_table[i].code;
			break;
		}
	} 

	for (i = 0 ; i < _exact_table_len ; i++) {
		if (_exact_table[i].min_num == err.err.min_num && _exact_table[i].maj_num == err.err.maj_num) {
			code = _exact_table[i].code;
			break;
		}
	} 

	if (err.err.desc == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_RUNTIME, "Failed to extract top-level HDF5 error description");
		return;
	}

	desc = polycap_strdup(err.err.desc);

	// Second, retrieve the bottom-most error description for additional info

	err.n = -1;

	if (H5Ewalk(H5E_DEFAULT, H5E_WALK_DOWNWARD, walk_cb, &err) < 0) {
		polycap_set_error_literal(error, POLYCAP_ERROR_RUNTIME, "Failed to walk HDF5 error stack");
		return;
	}

	desc_bottom = err.err.desc;
	if (desc_bottom == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_RUNTIME, "Failed to extract top-level HDF5 error description");
		return;
	}

	polycap_set_error(error, code, "%s (%s)", desc, desc_bottom);
	free(desc);
}

//===========================================
// Write data set in HDF5 file
static bool polycap_h5_write_dataset(hid_t file, int rank, hsize_t *dim, char *dataset_name, double *data, char *unitname, polycap_error **error) {
	herr_t status;
	hid_t dataset;
	hid_t dataspace, attr_id, attr_type, attr_dataspace_id; //handles

	//Describe size of the array and make fixed data space
	dataspace = H5Screate_simple(rank, dim, NULL);
	if (dataspace < 0) {
		set_exception(error);
		return false;
	}

	//Create new dataset within the HDF5 file with default creation properties
	dataset = H5Dcreate(file, dataset_name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (dataset < 0) {
		set_exception(error);
		return false;
	}

	//Write data to the dataset with default transfer properties
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	if(status < 0){
		set_exception(error);
		return false;
	}

	//Write unit attributes
	attr_dataspace_id = H5Screate(H5S_SCALAR);
	attr_type = H5Tcopy(H5T_C_S1);	
	if (H5Tset_size(attr_type,(hsize_t)strlen(unitname)) < 0) {
		set_exception(error);
		return false;
	}
	attr_id = H5Acreate(dataset, "Units", attr_type, attr_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
	if (attr_id < 0) {
		set_exception(error);
		return false;
	}
	if (H5Awrite(attr_id, attr_type, unitname) < 0) {
		set_exception(error);
		return false;
	}

	//Close release sources
	status = H5Sclose(attr_dataspace_id);
	if(status < 0){
		set_exception(error);
		return false;
	}
	status = H5Tclose(attr_type);
	if(status < 0){
		set_exception(error);
		return false;
	}
	status = H5Aclose(attr_id);
	if(status < 0){
		set_exception(error);
		return false;
	}
	status = H5Dclose(dataset);
	if(status < 0){
		set_exception(error);
		return false;
	}
	status = H5Sclose(dataspace);
	if(status < 0){
		set_exception(error);
		return false;
	}
	return true;
}
//===========================================
// Write efficiencies output in a hdf5 file
bool polycap_transmission_efficiencies_write_hdf5(polycap_transmission_efficiencies *efficiencies, const char *filename, polycap_error **error) {
	hid_t file, PC_Exit_id, PC_Start_id, Leaks_id, Recap_id, Input_id;
	hsize_t n_energies_temp, dim[2];
	double *data_temp;
	int j,k;

	tables_init();

	//argument sanity check
	if (filename == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_transmission_efficiencies_write_hdf5: filename cannot be NULL");
		return false;
	}
	if (efficiencies == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_transmission_efficiencies_write_hdf5: efficiencies cannot be NULL");
		return false;
	}
	//Create new HDF5 file using H5F_ACC_TRUNC and default creation and access properties
	file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); 
	if (file < 0) {
		set_exception(error);
		return false;
	}

	//Write energies array
	//Define temporary dataset dimension
	n_energies_temp = efficiencies->n_energies;
	if (!polycap_h5_write_dataset(file, 1, &n_energies_temp, "/Energies", efficiencies->energies,"keV", error))
		return false;

	//Write efficiencies array
	if (!polycap_h5_write_dataset(file, 1, &n_energies_temp, "/Transmission_Efficiencies", efficiencies->efficiencies,"a.u.", error))
		return false;

	//Write simulated polycap start coordinates
	//Create PC_Start group
	PC_Start_id = H5Gcreate2(file, "/PC_Start", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//Copy coordiante data to temporary array for straightforward HDF5 writing
	data_temp = malloc(sizeof(double)*efficiencies->images->i_exit*2);
	if(data_temp == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_MEMORY, strerror(errno));
		return false;
	}

	for(j=0;j<efficiencies->images->i_exit;j++){
		data_temp[j] = efficiencies->images->pc_start_coords[0][j];
		data_temp[j+efficiencies->images->i_exit] = efficiencies->images->pc_start_coords[1][j];
	}
	//Define temporary dataset dimension
	dim[0] = 2;
	dim[1] = efficiencies->images->i_exit;
	if (!polycap_h5_write_dataset(file, 2, dim, "/PC_Start/Coordinates", data_temp,"[cm,cm]", error))
		return false;
	//Free data_temp
	free(data_temp);
	
	//Write simulated polycap start direction
	//Copy direction data to temporary array for straightforward HDF5 writing
	data_temp = malloc(sizeof(double)*efficiencies->images->i_exit*2);
	if(data_temp == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_MEMORY, strerror(errno));
		return false;
	}
	for(j=0;j<efficiencies->images->i_exit;j++){
		data_temp[j] = efficiencies->images->pc_start_dir[0][j];
		data_temp[j+efficiencies->images->i_exit] = efficiencies->images->pc_start_dir[1][j];
	}
	//Define temporary dataset dimension
	dim[0] = 2;
	dim[1] = efficiencies->images->i_exit;
	if (!polycap_h5_write_dataset(file, 2, dim, "/PC_Start/Direction", data_temp,"[cm,cm]", error))
		return false;
	//Free data_temp
	free(data_temp);

	//Write simulated polycap start electric vectors
	//Copy direction data to temporary array for straightforward HDF5 writing
	data_temp = malloc(sizeof(double)*efficiencies->images->i_exit*2);
	if(data_temp == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_MEMORY, strerror(errno));
		return false;
	}
	for(j=0;j<efficiencies->images->i_exit;j++){
		data_temp[j] = efficiencies->images->pc_start_elecv[0][j];
		data_temp[j+efficiencies->images->i_exit] = efficiencies->images->pc_start_elecv[1][j];
	}
	//Define temporary dataset dimension
	dim[0] = 2;
	dim[1] = efficiencies->images->i_exit;
	if (!polycap_h5_write_dataset(file, 2, dim, "/PC_Start/Electric_Vector", data_temp,"[cm,cm]", error))
		return false;
	//Free data_temp
	free(data_temp);

	//Write simulated polycap exit coordinates
	//Create PC_Exit group
	PC_Exit_id = H5Gcreate2(file, "/PC_Exit", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//Copy coordinate data to temporary array for straightforward HDF5 writing
	data_temp = malloc(sizeof(double)*efficiencies->images->i_exit*3);
	if(data_temp == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_MEMORY, strerror(errno));
		return false;
	}
	for(j=0;j<efficiencies->images->i_exit;j++){
		data_temp[j] = efficiencies->images->pc_exit_coords[0][j];
		data_temp[j+efficiencies->images->i_exit] = efficiencies->images->pc_exit_coords[1][j];
		data_temp[j+efficiencies->images->i_exit*2] = efficiencies->images->pc_exit_coords[2][j];
	}
	//Define temporary dataset dimension
	dim[0] = 3;
	dim[1] = efficiencies->images->i_exit;
	if (!polycap_h5_write_dataset(file, 2, dim, "/PC_Exit/Coordinates", data_temp,"[cm,cm,cm]", error))
		return false;
	//Free data_temp
	free(data_temp);

	//Write n_reflections for each exited photon
	n_energies_temp = efficiencies->images->i_exit;
	data_temp = malloc(sizeof(double)*n_energies_temp);
	if(data_temp == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_MEMORY, strerror(errno));
		return false;
	}
	for(j=0; j<n_energies_temp; j++)
		data_temp[j] = (double)efficiencies->images->pc_exit_nrefl[j];
	if (!polycap_h5_write_dataset(file, 1, &n_energies_temp, "/PC_Exit/N_Reflections", data_temp,"a.u.", error))
		return false;
	//Free data_temp
	free(data_temp);

	//Write simulated polycap exit direction
	//Copy direction data to temporary array for straightforward HDF5 writing
	data_temp = malloc(sizeof(double)*efficiencies->images->i_exit*2);
	if(data_temp == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_MEMORY, strerror(errno));
		return false;
	}
	for(j=0;j<efficiencies->images->i_exit;j++){
		data_temp[j] = efficiencies->images->pc_exit_dir[0][j];
		data_temp[j+efficiencies->images->i_exit] = efficiencies->images->pc_exit_dir[1][j];
	}
	//Define temporary dataset dimension
	dim[0] = 2;
	dim[1] = efficiencies->images->i_exit;
	if (!polycap_h5_write_dataset(file, 2, dim, "/PC_Exit/Direction", data_temp,"[cm,cm]", error))
		return false;
	//Free data_temp
	free(data_temp);

	//Write simulated source start coordinates
	//Copy coordiante data to temporary array for straightforward HDF5 writing
	data_temp = malloc(sizeof(double)*efficiencies->images->i_exit*2);
	if(data_temp == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_MEMORY, strerror(errno));
		return false;
	}
	for(j=0;j<efficiencies->images->i_exit;j++){
		data_temp[j] = efficiencies->images->src_start_coords[0][j];
		data_temp[j+efficiencies->images->i_exit] = efficiencies->images->src_start_coords[1][j];
	}
	//Define temporary dataset dimension
	dim[0] = 2;
	dim[1] = efficiencies->images->i_exit;
	if (!polycap_h5_write_dataset(file, 2, dim, "/Source_Start_Coordinates", data_temp,"[cm,cm]", error))
		return false;

	//Free data_temp
	free(data_temp);

	//Write simulated polycap exit electric vectors
	//Copy direction data to temporary array for straightforward HDF5 writing
	data_temp = malloc(sizeof(double)*efficiencies->images->i_exit*2);
	if(data_temp == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_MEMORY, strerror(errno));
		return false;
	}
	for(j=0;j<efficiencies->images->i_exit;j++){
		data_temp[j] = efficiencies->images->pc_exit_elecv[0][j];
		data_temp[j+efficiencies->images->i_exit] = efficiencies->images->pc_exit_elecv[1][j];
	}
	//Define temporary dataset dimension
	dim[0] = 2;
	dim[1] = efficiencies->images->i_exit;
	if (!polycap_h5_write_dataset(file, 2, dim, "/PC_Exit/Electric_Vector", data_temp,"[cm,cm]", error))
		return false;
	//Free data_temp
	free(data_temp);

	//Write transmitted photon weights
	//Define temporary dataset dimension
	dim[1] = efficiencies->n_energies;
	dim[0] = efficiencies->images->i_exit;
	if (!polycap_h5_write_dataset(file, 2, dim, "/PC_Exit/Weights", efficiencies->images->exit_coord_weights,"[keV,a.u.]", error))
		return false;

	//Write transmitted photon traveled distance
	n_energies_temp = efficiencies->images->i_exit;
	if (!polycap_h5_write_dataset(file, 1, &n_energies_temp, "/PC_Exit/D_Travel", efficiencies->images->pc_exit_dtravel,"[cm]", error))
		return false;

	if(efficiencies->images->i_extleak > 0){
		//Write leak photons data
		//Make Leaks group
		Leaks_id = H5Gcreate2(file, "/ExternalLeaks", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			//write coordinates
		//Copy coordinate data to temporary array for straightforward HDF5 writing
		data_temp = malloc(sizeof(double)*efficiencies->images->i_extleak*3);
		if(data_temp == NULL){
			polycap_set_error_literal(error, POLYCAP_ERROR_MEMORY, strerror(errno));
			return false;
		}
		for(j=0;j<efficiencies->images->i_extleak;j++){
			data_temp[j] = efficiencies->images->extleak_coords[0][j];
			data_temp[j+efficiencies->images->i_extleak] = efficiencies->images->extleak_coords[1][j];
			data_temp[j+efficiencies->images->i_extleak*2] = efficiencies->images->extleak_coords[2][j];
		}
		//Define temporary dataset dimension
		dim[0] = 3;
		dim[1] = efficiencies->images->i_extleak;
		if (!polycap_h5_write_dataset(file, 2, dim, "/ExternalLeaks/Coordinates", data_temp,"[cm,cm,cm]", error))
			return false;
		//Free data_temp
		free(data_temp);
			//write direction
		//Copy direction data to temporary array for straightforward HDF5 writing
		data_temp = malloc(sizeof(double)*efficiencies->images->i_extleak*2);
		if(data_temp == NULL){
			polycap_set_error_literal(error, POLYCAP_ERROR_MEMORY, strerror(errno));
			return false;
		}
		for(j=0;j<efficiencies->images->i_extleak;j++){
			data_temp[j] = efficiencies->images->extleak_dir[0][j];
			data_temp[j+efficiencies->images->i_extleak] = efficiencies->images->extleak_dir[1][j];
		}
		//Define temporary dataset dimension
		dim[0] = 2;
		dim[1] = efficiencies->images->i_extleak;
		if (!polycap_h5_write_dataset(file, 2, dim, "/ExternalLeaks/Direction", data_temp,"[cm,cm]", error))
			return false;
		//Free data_temp
		free(data_temp);
			//Write leaked photon weights
		//Define temporary dataset dimension
		dim[1] = efficiencies->n_energies;
		dim[0] = efficiencies->images->i_extleak;
		if (!polycap_h5_write_dataset(file, 2, dim, "/ExternalLeaks/Weights", efficiencies->images->extleak_coord_weights,"[keV,a.u.]", error))
			return false;
		//Save weight average as function of energy
		data_temp = malloc(sizeof(double)*efficiencies->n_energies);
		for(j=0; j<efficiencies->n_energies; j++){
			data_temp[j] = 0.0;
			for(k=0; k<efficiencies->images->i_extleak; k++){
				data_temp[j] += efficiencies->images->extleak_coord_weights[k*efficiencies->n_energies+j];
			}
			data_temp[j] = data_temp[j] / (double)efficiencies->images->i_start;
		}
		n_energies_temp = efficiencies->n_energies;
		if (!polycap_h5_write_dataset(file, 1, &n_energies_temp, "/ExternalLeaks/Weight_Total", data_temp,"a.u.", error))
			return false;
		//Free data_temp
		free(data_temp);
		//Write n_reflections for each leaked photon
		n_energies_temp = efficiencies->images->i_extleak;
		data_temp = malloc(sizeof(double)*n_energies_temp);
		if(data_temp == NULL){
			polycap_set_error_literal(error, POLYCAP_ERROR_MEMORY, strerror(errno));
			return false;
		}
		for(j=0; j<n_energies_temp; j++)
			data_temp[j] = (double)efficiencies->images->extleak_n_refl[j];
		if (!polycap_h5_write_dataset(file, 1, &n_energies_temp, "/ExternalLeaks/N_Reflections", data_temp,"a.u.", error))
			return false;
		//Free data_temp
		free(data_temp);

		if (H5Gclose(Leaks_id) < 0)
			set_exception(error);
	}

	//Write intleak photons data
	//Make Recap group
	if(efficiencies->images->i_intleak > 0){
		Recap_id = H5Gcreate2(file, "/InternalLeaks", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			//write coordinates
		//Copy coordinate data to temporary array for straightforward HDF5 writing
		data_temp = malloc(sizeof(double)*efficiencies->images->i_intleak*3);
		if(data_temp == NULL){
			polycap_set_error_literal(error, POLYCAP_ERROR_MEMORY, strerror(errno));
			return false;
		}
		for(j=0;j<efficiencies->images->i_intleak;j++){
			data_temp[j] = efficiencies->images->intleak_coords[0][j];
			data_temp[j+efficiencies->images->i_intleak] = efficiencies->images->intleak_coords[1][j];
			data_temp[j+efficiencies->images->i_intleak*2] = efficiencies->images->intleak_coords[2][j];
		}
		//Define temporary dataset dimension
		dim[0] = 3;
		dim[1] = efficiencies->images->i_intleak;
		if (!polycap_h5_write_dataset(file, 2, dim, "/InternalLeaks/Coordinates", data_temp,"[cm,cm,cm]", error))
			return false;
		//Free data_temp
		free(data_temp);
			//write direction
		//Copy direction data to temporary array for straightforward HDF5 writing
		data_temp = malloc(sizeof(double)*efficiencies->images->i_intleak*2);
		if(data_temp == NULL){
			polycap_set_error_literal(error, POLYCAP_ERROR_MEMORY, strerror(errno));
			return false;
		}
		for(j=0;j<efficiencies->images->i_intleak;j++){
			data_temp[j] = efficiencies->images->intleak_dir[0][j];
			data_temp[j+efficiencies->images->i_intleak] = efficiencies->images->intleak_dir[1][j];
		}
		//Define temporary dataset dimension
		dim[0] = 2;
		dim[1] = efficiencies->images->i_intleak;
		if (!polycap_h5_write_dataset(file, 2, dim, "/InternalLeaks/Direction", data_temp,"[cm,cm]", error))
			return false;
		//Free data_temp
		free(data_temp);
		//Write intleak electric vectors
		data_temp = malloc(sizeof(double)*efficiencies->images->i_intleak*2);
		if(data_temp == NULL){
			polycap_set_error_literal(error, POLYCAP_ERROR_MEMORY, strerror(errno));
			return false;
		}
		for(j=0;j<efficiencies->images->i_intleak;j++){
			data_temp[j] = efficiencies->images->intleak_elecv[0][j];
			data_temp[j+efficiencies->images->i_intleak] = efficiencies->images->intleak_elecv[1][j];
		}
		//Define temporary dataset dimension
		dim[0] = 2;
		dim[1] = efficiencies->images->i_intleak;
		if (!polycap_h5_write_dataset(file, 2, dim, "/InternalLeaks/Electric_Vector", data_temp,"[cm,cm]", error))
			return false;
		//Free data_temp
		free(data_temp);
		//Write intleak photon weights
		//Define temporary dataset dimension
		dim[1] = efficiencies->n_energies;
		dim[0] = efficiencies->images->i_intleak;
		if (!polycap_h5_write_dataset(file, 2, dim, "/InternalLeaks/Weights", efficiencies->images->intleak_coord_weights,"[keV,a.u.]", error))
			return false;
		//Save weight average as function of energy
		data_temp = malloc(sizeof(double)*efficiencies->n_energies);
		for(j=0; j<efficiencies->n_energies; j++){
			data_temp[j] = 0.0;
			for(k=0; k<efficiencies->images->i_intleak; k++){
				data_temp[j] += efficiencies->images->intleak_coord_weights[k*efficiencies->n_energies+j];
			}
			data_temp[j] = data_temp[j] / (double)efficiencies->images->i_start;
		}
		n_energies_temp = efficiencies->n_energies;
		if (!polycap_h5_write_dataset(file, 1, &n_energies_temp, "/InternalLeaks/Weight_Total", data_temp,"a.u.", error))
			return false;
		//Free data_temp
		free(data_temp);
		//Write	n_reflections for each intleak photon
		n_energies_temp = efficiencies->images->i_intleak;
		data_temp = malloc(sizeof(double)*n_energies_temp);
		if(data_temp == NULL){
			polycap_set_error_literal(error, POLYCAP_ERROR_MEMORY, strerror(errno));
			return false;
		}
		for(j=0; j<n_energies_temp; j++)
			data_temp[j] = (double)efficiencies->images->intleak_n_refl[j];
		if (!polycap_h5_write_dataset(file, 1, &n_energies_temp, "/InternalLeaks/N_Reflections", data_temp,"a.u.", error))
			return false;
		//Free data_temp
		free(data_temp);

		if (H5Gclose(Recap_id) < 0)
			set_exception(error);
	}

	//Write Input parameters
	//Make Input group
	Input_id = H5Gcreate2(file, "/Input", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//Copy direction data to temporary array for straightforward HDF5 writing
	data_temp = malloc(sizeof(double)*efficiencies->source->description->profile->nmax*2);
	if(data_temp == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_MEMORY, strerror(errno));
		return false;
	}
	for(j=0;j<efficiencies->source->description->profile->nmax;j++){
		data_temp[j] = efficiencies->source->description->profile->z[j];
		data_temp[j+efficiencies->source->description->profile->nmax] = efficiencies->source->description->profile->ext[j];
	}
	//Define temporary dataset dimension
	dim[0] = 2;
	dim[1] = efficiencies->source->description->profile->nmax;
	if (!polycap_h5_write_dataset(file, 2, dim, "/Input/PC_Shape", data_temp,"[cm,cm]", error))
		return false;
	//Free data_temp
	free(data_temp);

	//Copy direction data to temporary array for straightforward HDF5 writing
	data_temp = malloc(sizeof(double)*efficiencies->source->description->profile->nmax*2);
	if(data_temp == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_MEMORY, strerror(errno));
		return false;
	}
	for(j=0 ; j < efficiencies->source->description->profile->nmax ; j++){
		data_temp[j] = efficiencies->source->description->profile->z[j];
		data_temp[j+ efficiencies->source->description->profile->nmax] = efficiencies->source->description->profile->cap[j];
	}
	//Define temporary dataset dimension
	dim[0] = 2;
	dim[1] = efficiencies->source->description->profile->nmax;
	if (!polycap_h5_write_dataset(file, 2, dim, "/Input/Cap_Shape", data_temp,"[cm,cm]", error))
		return false;
	//Free data_temp
	free(data_temp);
	
	//Write ncap and other input parameters
	n_energies_temp = 1;
	data_temp = malloc(sizeof(double));
	if(data_temp == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_MEMORY, strerror(errno));
		return false;
	}
	*data_temp = (double)efficiencies->source->description->n_cap;
	if (!polycap_h5_write_dataset(file, 1, &n_energies_temp, "/Input/N_Capillaries", data_temp,"a.u.", error))
		return false;
	free(data_temp);
	if (!polycap_h5_write_dataset(file, 1, &n_energies_temp, "/Input/Surface_Roughness", &efficiencies->source->description->sig_rough,"Angstrom", error))
		return false;
	if (!polycap_h5_write_dataset(file, 1, &n_energies_temp, "/Input/Open_Area", &efficiencies->source->description->open_area,"a.u.", error))
		return false;
	dim[0] = 2;
	dim[1] = efficiencies->source->description->nelem;
	data_temp = malloc(sizeof(double)*efficiencies->source->description->nelem*2);
	if(data_temp == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_MEMORY, strerror(errno));
		return false;
	}
	for(j=0 ; j < efficiencies->source->description->nelem ; j++){
		data_temp[j] = efficiencies->source->description->iz[j];
		data_temp[j+ efficiencies->source->description->nelem] = efficiencies->source->description->wi[j];
	}
	if (!polycap_h5_write_dataset(file, 2, dim, "/Input/PC_Composition", data_temp,"[Z,w%]", error))
		return false;
	free(data_temp);
	if (!polycap_h5_write_dataset(file, 1, &n_energies_temp, "/Input/PC_Density", &efficiencies->source->description->density,"g/cm3", error))
		return false;
	if (!polycap_h5_write_dataset(file, 1, &n_energies_temp, "/Input/Src_PC_Dist", &efficiencies->source->d_source,"cm", error))
		return false;


	//Close Group access
	if (H5Gclose(PC_Exit_id) < 0)
		set_exception(error);
	if (H5Gclose(PC_Start_id) < 0)
		set_exception(error);
	if (H5Gclose(Input_id) < 0)
		set_exception(error);

	//Close file
	if (H5Fclose(file) < 0)
		set_exception(error);
	
	return true;
}
//===========================================
bool polycap_transmission_efficiencies_get_start_data(polycap_transmission_efficiencies *efficiencies, int64_t *n_start, int64_t *n_exit, polycap_vector3 **start_coords, polycap_vector3 **start_direction, polycap_vector3 **start_elecv, polycap_vector3 **src_start_coords, polycap_error **error)
{
	int i;

	if (efficiencies == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_start_data: efficiencies cannot be NULL");
		return false;
	}
	if (efficiencies->images == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_start_data: source->images cannot be NULL");
		return false;
	}

	*n_start = efficiencies->images->i_start;
	if (efficiencies->images->i_start == 0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_start_data: no photon start events in efficiencies");
		return false;
	}
	*n_exit = efficiencies->images->i_exit;
	if (efficiencies->images->i_exit == 0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_start_data: no photon exit events in efficiencies");
		return false;
	}

	*start_coords = malloc(sizeof(polycap_vector3) * efficiencies->images->i_exit);
	if (*start_coords == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_start_data: could not allocate memory for start_coords -> %s", strerror(errno));
		return false;
	}
	*start_direction = malloc(sizeof(polycap_vector3) * efficiencies->images->i_exit);
	if (*start_direction == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_start_data: could not allocate memory for start_direction -> %s", strerror(errno));
		return false;
	}
	*start_elecv = malloc(sizeof(polycap_vector3) * efficiencies->images->i_exit);
	if (*start_elecv == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_start_data: could not allocate memory for start_elecv -> %s", strerror(errno));
		return false;
	}
	*src_start_coords = malloc(sizeof(polycap_vector3) * efficiencies->images->i_exit);
	if (*src_start_coords == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_start_data: could not allocate memory for src_start_coords -> %s", strerror(errno));
		return false;
	}

	for(i = 0; i < efficiencies->images->i_exit; i++) {
		(*start_coords)[i].x = efficiencies->images->pc_start_coords[0][i];
		(*start_coords)[i].y = efficiencies->images->pc_start_coords[1][i];
		(*start_coords)[i].z = 0.;

		(*start_direction)[i].x = efficiencies->images->pc_start_dir[0][i];
		(*start_direction)[i].y = efficiencies->images->pc_start_dir[1][i];
		(*start_direction)[i].z = sqrt(1.- (*start_direction)[i].x * (*start_direction)[i].x - (*start_direction)[i].y * (*start_direction)[i].y);

		(*start_elecv)[i].x = efficiencies->images->pc_start_elecv[0][i];
		(*start_elecv)[i].y = efficiencies->images->pc_start_elecv[1][i];
		(*start_elecv)[i].z = sqrt(1.- (*start_elecv)[i].x * (*start_elecv)[i].x - (*start_elecv)[i].y * (*start_elecv)[i].y);

		(*src_start_coords)[i].x = efficiencies->images->src_start_coords[0][i];
		(*src_start_coords)[i].y = efficiencies->images->src_start_coords[1][i];
		(*src_start_coords)[i].z = 0.;
	}

	return true;
}
//===========================================
bool polycap_transmission_efficiencies_get_exit_data(polycap_transmission_efficiencies *efficiencies, int64_t *n_exit, polycap_vector3 **exit_coords, polycap_vector3 **exit_direction, polycap_vector3 **exit_elecv, int64_t **n_refl, double **d_travel, size_t *n_energies, double ***exit_weights, polycap_error **error)
{
	int i,j;

	if (efficiencies == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_exit_data: efficiencies cannot be NULL");
		return false;
	}
	if (efficiencies->images == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_exit_data: source->images cannot be NULL");
		return false;
	}

	*n_exit = efficiencies->images->i_exit;
	*n_energies = efficiencies->n_energies;
	if (efficiencies->images->i_start == 0){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_exit_data: no photon start events in efficiencies");
		return false;
	}

	*exit_coords = malloc(sizeof(polycap_vector3) * efficiencies->images->i_exit);
	if (*exit_coords == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_exit_data: could not allocate memory for exit_coords -> %s", strerror(errno));
		return false;
	}
	*exit_direction = malloc(sizeof(polycap_vector3) * efficiencies->images->i_exit);
	if (*exit_direction == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_exit_data: could not allocate memory for exit_direction -> %s", strerror(errno));
		return false;
	}
	*exit_elecv = malloc(sizeof(polycap_vector3) * efficiencies->images->i_exit);
	if (*exit_elecv == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_exit_data: could not allocate memory for exit_elecv -> %s", strerror(errno));
		return false;
	}
	*n_refl = malloc(sizeof(int64_t) * efficiencies->images->i_exit);
	if (*n_refl == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_exit_data: could not allocate memory for n_refl -> %s", strerror(errno));
		return false;
	}
	*d_travel = malloc(sizeof(double) * efficiencies->images->i_exit);
	if (*d_travel == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_exit_data: could not allocate memory for d_travel -> %s", strerror(errno));
		return false;
	}
	*exit_weights = malloc(sizeof(double*) * efficiencies->images->i_exit);
	if (*exit_weights == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_exit_data: could not allocate memory for exit_weights -> %s", strerror(errno));
		return false;
	}


	for(i = 0; i < efficiencies->images->i_exit; i++) {
		(*n_refl)[i] = efficiencies->images->pc_exit_nrefl[i];
		(*d_travel)[i] = efficiencies->images->pc_exit_dtravel[i];

		(*exit_coords)[i].x = efficiencies->images->pc_exit_coords[0][i];
		(*exit_coords)[i].y = efficiencies->images->pc_exit_coords[1][i];
		(*exit_coords)[i].z = efficiencies->images->pc_exit_coords[2][i];

		(*exit_direction)[i].x = efficiencies->images->pc_exit_dir[0][i];
		(*exit_direction)[i].y = efficiencies->images->pc_exit_dir[1][i];
		(*exit_direction)[i].z = sqrt(1.-(*exit_direction)[i].x*(*exit_direction)[i].x - (*exit_direction)[i].y*(*exit_direction)[i].y);

		(*exit_elecv)[i].x = efficiencies->images->pc_exit_elecv[0][i];
		(*exit_elecv)[i].y = efficiencies->images->pc_exit_elecv[1][i];
		(*exit_elecv)[i].z = sqrt(1.-(*exit_elecv)[i].x*(*exit_elecv)[i].x - (*exit_elecv)[i].y*(*exit_elecv)[i].y);

		(*exit_weights)[i] = malloc(sizeof(double) * efficiencies->n_energies);
		if ( (*exit_weights)[i] == NULL){
			polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_exit_data: could not allocate memory for (*exit_weights)[i] -> %s", strerror(errno));
			return false;
		}
		for(j=0; j < efficiencies->n_energies; j++){
			(*exit_weights)[i][j] = efficiencies->images->exit_coord_weights[i*efficiencies->n_energies+j];
		}
	}

	return true;
}
//===========================================
bool polycap_transmission_efficiencies_get_extleak_data(polycap_transmission_efficiencies *efficiencies, polycap_leak ***leaks, int64_t *n_leaks, polycap_error **error)
{
	int i,j;
	polycap_vector3 temp;

	if (efficiencies == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_extleak_data: efficiencies cannot be NULL");
		return false;
	}
	if (efficiencies->images == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_extleak_data: source->images cannot be NULL");
		return false;
	}

	*n_leaks = efficiencies->images->i_extleak;
	//fprintf(stderr, "C: n_extleak: %lld\n", efficiencies->images->i_extleak);
	if (efficiencies->images->i_extleak == 0){
		*leaks = NULL;
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_extleak_data: no extleak events in efficiencies");
		return false;
	}

	*leaks = malloc(sizeof(polycap_leak*) * efficiencies->images->i_extleak);
	if (*leaks == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_extleak_data: could not allocate memory for leaks -> %s", strerror(errno));
		return false;
	}

	for(i = 0; i < efficiencies->images->i_extleak; i++) {
		(*leaks)[i] = malloc(sizeof(polycap_leak));
		if ((*leaks)[i] == NULL){
			polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_extleak_data: could not allocate memory for (*leaks)[i] -> %s", strerror(errno));
			return false;
		}
		temp.x = efficiencies->images->extleak_coords[0][i];
		temp.y = efficiencies->images->extleak_coords[1][i];
		temp.z = efficiencies->images->extleak_coords[2][i];
		memcpy(&(*leaks)[i]->coords, &temp, sizeof(polycap_vector3));
		temp.x = efficiencies->images->extleak_dir[0][i];
		temp.y = efficiencies->images->extleak_dir[1][i];
		temp.z = sqrt(1.-temp.x*temp.x - temp.y*temp.y);
		memcpy(&(*leaks)[i]->direction, &temp, sizeof(polycap_vector3));
		//temp.x = efficiencies->images->extleak_elecv[0][i];
		//temp.y = efficiencies->images->extleak_elecv[1][i];
		//temp.z = sqrt(1.-temp.x*temp.x - temp.y*temp.y);
		//memcpy(&(*leaks)[i]->elecv, &temp, sizeof(polycap_vector3));
		
		(*leaks)[i]->n_energies = efficiencies->n_energies;
		(*leaks)[i]->n_refl = efficiencies->images->extleak_n_refl[i];
		
		(*leaks)[i]->weight = malloc(sizeof(double) * efficiencies->n_energies);
		if ( (*leaks)[i]->weight == NULL){
			polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_extleak_data: could not allocate memory for (*leaks[i])->weight -> %s", strerror(errno));
			return false;
		}

		for(j=0; j < efficiencies->n_energies; j++){
			(*leaks)[i]->weight[j] = efficiencies->images->extleak_coord_weights[i*efficiencies->n_energies+j];
		}

	}

	return true;
}
//===========================================
bool polycap_transmission_efficiencies_get_intleak_data(polycap_transmission_efficiencies *efficiencies, polycap_leak ***leaks, int64_t *n_leaks, polycap_error **error)
{
	int i,j;
	polycap_vector3 temp;

	if (efficiencies == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_intleak_data: efficiencies cannot be NULL");
		return false;
	}
	if (efficiencies->images == NULL){
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_get_intleak_data: source->images cannot be NULL");
		return false;
	}

	*n_leaks = efficiencies->images->i_intleak;
	if (efficiencies->images->i_intleak == 0){
		*leaks = NULL;
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_source_getintleak_data: no intleak events in efficiencies");
		return false;
	}

	*leaks = malloc(sizeof(polycap_leak*) * efficiencies->images->i_intleak);
	if (*leaks == NULL){
		polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_intleak_data: could not allocate memory for leaks -> %s", strerror(errno));
		return false;
	}

	for(i = 0; i < efficiencies->images->i_intleak; i++) {
		(*leaks)[i] = malloc(sizeof(polycap_leak));
		if ((*leaks)[i] == NULL){
			polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_intleak_data: could not allocate memory for (*leaks)[i] -> %s", strerror(errno));
			return false;
		}
		temp.x = efficiencies->images->intleak_coords[0][i];
		temp.y = efficiencies->images->intleak_coords[1][i];
		temp.z = efficiencies->images->intleak_coords[2][i];
		memcpy(&(*leaks)[i]->coords, &temp, sizeof(polycap_vector3));
		temp.x = efficiencies->images->intleak_dir[0][i];
		temp.y = efficiencies->images->intleak_dir[1][i];
		temp.z = sqrt(1.-temp.x*temp.x - temp.y*temp.y);
		memcpy(&(*leaks)[i]->direction, &temp, sizeof(polycap_vector3));
		temp.x = efficiencies->images->intleak_elecv[0][i];
		temp.y = efficiencies->images->intleak_elecv[1][i];
		temp.z = sqrt(1.-temp.x*temp.x - temp.y*temp.y);
		memcpy(&(*leaks)[i]->elecv, &temp, sizeof(polycap_vector3));
		
		(*leaks)[i]->n_energies = efficiencies->n_energies;
		(*leaks)[i]->n_refl = efficiencies->images->intleak_n_refl[i];

		(*leaks)[i]->weight = malloc(sizeof(double) * efficiencies->n_energies);
		if ( (*leaks)[i]->weight == NULL){
			polycap_set_error(error, POLYCAP_ERROR_MEMORY, "polycap_source_get_intleak_data: could not allocate memory for (*leaks[i])->weight -> %s", strerror(errno));
			return false;
		}
		for(j=0; j < efficiencies->n_energies; j++){
			(*leaks)[i]->weight[j] = efficiencies->images->intleak_coord_weights[i*efficiencies->n_energies+j];
		}
	}

	return true;
}
//===========================================
void polycap_images_free(struct _polycap_images *images)
{
	if (images == NULL)
		return;

	if (images->src_start_coords[0])
		free(images->src_start_coords[0]);
	if (images->src_start_coords[1])
		free(images->src_start_coords[1]);
	if (images->pc_start_coords[0])
		free(images->pc_start_coords[0]);
	if (images->pc_start_coords[1])
		free(images->pc_start_coords[1]);
	if (images->pc_start_dir[0])
		free(images->pc_start_dir[0]);
	if (images->pc_start_dir[1])
		free(images->pc_start_dir[1]);
	if (images->pc_start_elecv[0])
		free(images->pc_start_elecv[0]);
	if (images->pc_start_elecv[1])
		free(images->pc_start_elecv[1]);
	if (images->pc_exit_coords[0])
		free(images->pc_exit_coords[0]);
	if (images->pc_exit_coords[1])
		free(images->pc_exit_coords[1]);
	if (images->pc_exit_coords[2])
		free(images->pc_exit_coords[2]);
	if (images->pc_exit_dir[0])
		free(images->pc_exit_dir[0]);
	if (images->pc_exit_dir[1])
		free(images->pc_exit_dir[1]);
	if (images->pc_exit_elecv[0])
		free(images->pc_exit_elecv[0]);
	if (images->pc_exit_elecv[1])
		free(images->pc_exit_elecv[1]);
	if (images->pc_exit_nrefl)
		free(images->pc_exit_nrefl);
	if (images->pc_exit_dtravel)
		free(images->pc_exit_dtravel);
	if (images->exit_coord_weights)
		free(images->exit_coord_weights);
	if (images->extleak_coords[0])
		free(images->extleak_coords[0]);
	if (images->extleak_coords[1])
		free(images->extleak_coords[1]);
	if (images->extleak_coords[2])
		free(images->extleak_coords[2]);
	if (images->extleak_dir[0])
		free(images->extleak_dir[0]);
	if (images->extleak_dir[1])
		free(images->extleak_dir[1]);
	if (images->extleak_coord_weights)
		free(images->extleak_coord_weights);
	if (images->extleak_n_refl)
		free(images->extleak_n_refl);
	if (images->intleak_coords[0])
		free(images->intleak_coords[0]);
	if (images->intleak_coords[1])
		free(images->intleak_coords[1]);
	if (images->intleak_coords[2])
		free(images->intleak_coords[2]);
	if (images->intleak_dir[0])
		free(images->intleak_dir[0]);
	if (images->intleak_dir[1])
		free(images->intleak_dir[1]);
	if (images->intleak_elecv[0])
		free(images->intleak_elecv[0]);
	if (images->intleak_elecv[1])
		free(images->intleak_elecv[1]);
	if (images->intleak_coord_weights)
		free(images->intleak_coord_weights);
	if (images->intleak_n_refl)
		free(images->intleak_n_refl);
	free(images);

}
//===========================================
void polycap_transmission_efficiencies_free(polycap_transmission_efficiencies *efficiencies)
{
	if (efficiencies == NULL)
		return;
	if (efficiencies->energies)
		free(efficiencies->energies);
	if (efficiencies->efficiencies)
		free(efficiencies->efficiencies);
	if (efficiencies->images) {
		polycap_images_free(efficiencies->images);
	}
	free(efficiencies);
}

bool polycap_transmission_efficiencies_get_data(polycap_transmission_efficiencies *efficiencies, size_t *n_energies, double **energies_arr, double **efficiencies_arr, polycap_error **error) {
	if (efficiencies == NULL) {
		polycap_set_error_literal(error, POLYCAP_ERROR_INVALID_ARGUMENT, "polycap_transmission_efficiencies_get_data: efficiencies cannot be NULL");
		return false;
	}

	if (n_energies)
		*n_energies = efficiencies->n_energies;

	if (energies_arr) {
		*energies_arr = malloc(sizeof(double) * efficiencies->n_energies);
		memcpy(*energies_arr, efficiencies->energies, sizeof(double) * efficiencies->n_energies);
	}

	if (efficiencies_arr) {
		*efficiencies_arr = malloc(sizeof(double) * efficiencies->n_energies);
		memcpy(*efficiencies_arr, efficiencies->efficiencies, sizeof(double) * efficiencies->n_energies);
	}

	return true;
}

void polycap_free(void *data) {
	if (data)
		free(data);
}
