#include "polycap-private.h"
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include "hdf5.h"

//===========================================
// Write data set in HDF5 file
void polycap_h5_write_dataset(hid_t file, int rank, hsize_t *dim, char *dataset_name, double *data, char *unitname)
{
	herr_t status;
	hid_t dataset;
	hid_t dataspace, attr_id, attr_type, attr_dataspace_id; //handles

	//Describe size of the array and make fixed data space
	dataspace = H5Screate_simple(rank, dim, NULL);

	//Create new dataset within the HDF5 file with default creation properties
	dataset = H5Dcreate2(file, dataset_name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	//Write data to the dataset with default transfer properties
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	if(status < 0){
		printf("Error: H5D write finished on error.\n");
		exit(4);
	}

	//Write unit attributes
	attr_dataspace_id = H5Screate(H5S_SCALAR);
	attr_type = H5Tcopy(H5T_C_S1);	
	H5Tset_size(attr_type,(hsize_t)strlen(unitname));
	attr_id = H5Acreate2(dataset, "Units", attr_type, attr_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr_id, attr_type, unitname);

	//Close release sources
	status = H5Sclose(attr_dataspace_id);
	if(status < 0){
		printf("Error: H5S close finished on error.\n");
		exit(5);
	}
	status = H5Tclose(attr_type);
	if(status < 0){
		printf("Error: H5T close finished on error.\n");
		exit(7);
	}
	status = H5Aclose(attr_id);
	if(status < 0){
		printf("Error: H5A close finished on error.\n");
		exit(8);
	}
	status = H5Dclose(dataset);
	if(status < 0){
		printf("Error: H5D close finished on error.\n");
		exit(6);
	}
	status = H5Sclose(dataspace);
	if(status < 0){
		printf("Error: H5S close finished on error.\n");
		exit(5);
	}

}
//===========================================
// Write efficiencies output in a hdf5 file
void polycap_transmission_efficiencies_write_hdf5(polycap_transmission_efficiencies *efficiencies, const char *filename)
{
	hid_t file;
	hsize_t n_energies_temp, dim[2];
	double *data_temp;
	int j;

	//Create new HDF5 file using H5F_ACC_TRUNC and default creation and access properties
	file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); 

	//Write energies array
	//Define temporary dataset dimension
	n_energies_temp = efficiencies->n_energies;
	polycap_h5_write_dataset(file, 1, &n_energies_temp, "/Energies", efficiencies->energies,"keV");
	//Write efficiencies array
	polycap_h5_write_dataset(file, 1, &n_energies_temp, "/Transmission_Efficiencies", efficiencies->efficiencies,"a.u.");

	//Write simulated polycap start coordinates
	//Copy coordiante data to temporary array for straightforward HDF5 writing
	data_temp = malloc(sizeof(double)*efficiencies->images->i_start*2);
	if(data_temp == NULL){
		printf("Could not allocate data_temp memory.\n");
		exit(1);
	}
	for(j=0;j<efficiencies->images->i_start;j++){
		data_temp[j] = efficiencies->images->pc_start_coords[0][j];
		data_temp[j+efficiencies->images->i_start] = efficiencies->images->pc_start_coords[1][j];
	}
	//Define temporary dataset dimension
	dim[0] = 2;
	dim[1] = efficiencies->images->i_start;
	polycap_h5_write_dataset(file, 2, dim, "/PC_Start_Coordinates", data_temp,"cm");
	//Free data_temp
	free(data_temp);
	
	//Write simulated polycap start direction
	//Copy direction data to temporary array for straightforward HDF5 writing
	data_temp = malloc(sizeof(double)*efficiencies->images->i_start*2);
	if(data_temp == NULL){
		printf("Could not allocate data_temp memory.\n");
		exit(1);
	}
	for(j=0;j<efficiencies->images->i_start;j++){
		data_temp[j] = efficiencies->images->pc_start_dir[0][j];
		data_temp[j+efficiencies->images->i_start] = efficiencies->images->pc_start_dir[1][j];
	}
	//Define temporary dataset dimension
	dim[0] = 2;
	dim[1] = efficiencies->images->i_start;
	polycap_h5_write_dataset(file, 2, dim, "/PC_Start_Direction", data_temp,"cm");
	//Free data_temp
	free(data_temp);

	//Write simulated polycap exit coordinates
	//Copy coordiante data to temporary array for straightforward HDF5 writing
	data_temp = malloc(sizeof(double)*efficiencies->images->i_exit*2);
	if(data_temp == NULL){
		printf("Could not allocate data_temp memory.\n");
		exit(1);
	}
	for(j=0;j<efficiencies->images->i_exit;j++){
		data_temp[j] = efficiencies->images->pc_exit_coords[0][j];
		data_temp[j+efficiencies->images->i_exit] = efficiencies->images->pc_exit_coords[1][j];
	}
	//Define temporary dataset dimension
	dim[0] = 2;
	dim[1] = efficiencies->images->i_exit;
	polycap_h5_write_dataset(file, 2, dim, "/PC_Exit_Coordinates", data_temp,"cm");
	//Free data_temp
	free(data_temp);
	
	//Write simulated polycap exit direction
	//Copy direction data to temporary array for straightforward HDF5 writing
	data_temp = malloc(sizeof(double)*efficiencies->images->i_exit*2);
	if(data_temp == NULL){
		printf("Could not allocate data_temp memory.\n");
		exit(1);
	}
	for(j=0;j<efficiencies->images->i_exit;j++){
		data_temp[j] = efficiencies->images->pc_exit_dir[0][j];
		data_temp[j+efficiencies->images->i_exit] = efficiencies->images->pc_exit_dir[1][j];
	}
	//Define temporary dataset dimension
	dim[0] = 2;
	dim[1] = efficiencies->images->i_exit;
	polycap_h5_write_dataset(file, 2, dim, "/PC_Exit_Direction", data_temp,"cm");
	//Free data_temp
	free(data_temp);

	//Write simulated source start coordinates
	//Copy coordiante data to temporary array for straightforward HDF5 writing
	data_temp = malloc(sizeof(double)*efficiencies->images->i_start*2);
	if(data_temp == NULL){
		printf("Could not allocate data_temp memory.\n");
		exit(1);
	}
	for(j=0;j<efficiencies->images->i_start;j++){
		data_temp[j] = efficiencies->images->src_start_coords[0][j];
		data_temp[j+efficiencies->images->i_start] = efficiencies->images->src_start_coords[1][j];
	}
	//Define temporary dataset dimension
	dim[0] = 2;
	dim[1] = efficiencies->images->i_start;
	polycap_h5_write_dataset(file, 2, dim, "/Source_Start_Coordinates", data_temp,"cm");
	//Free data_temp
	free(data_temp);

	//Write simulated source start coordinates
	//Define temporary dataset dimension
	dim[1] = efficiencies->n_energies;
	dim[0] = efficiencies->images->i_exit;
	polycap_h5_write_dataset(file, 2, dim, "/Exit_Coordinate_Weights", efficiencies->images->exit_coord_weights,"a.u.");

	//Close file
	H5Fclose(file);

	printf("%s was written.\n",filename);
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
		if (efficiencies->images->src_start_coords[0])
			free(efficiencies->images->src_start_coords[0]);
		if (efficiencies->images->src_start_coords[1])
			free(efficiencies->images->src_start_coords[1]);
		if (efficiencies->images->pc_start_coords[0])
			free(efficiencies->images->pc_start_coords[0]);
		if (efficiencies->images->pc_start_coords[1])
			free(efficiencies->images->pc_start_coords[1]);
		if (efficiencies->images->pc_start_dir[0])
			free(efficiencies->images->pc_start_dir[0]);
		if (efficiencies->images->pc_start_dir[1])
			free(efficiencies->images->pc_start_dir[1]);
		if (efficiencies->images->pc_exit_coords[0])
			free(efficiencies->images->pc_exit_coords[0]);
		if (efficiencies->images->pc_exit_coords[1])
			free(efficiencies->images->pc_exit_coords[1]);
		if (efficiencies->images->pc_exit_dir[0])
			free(efficiencies->images->pc_exit_dir[0]);
		if (efficiencies->images->pc_exit_dir[1])
			free(efficiencies->images->pc_exit_dir[1]);
		if (efficiencies->images->exit_coord_weights)
			free(efficiencies->images->exit_coord_weights);
		free(efficiencies->images);
	}
	free(efficiencies);
}
