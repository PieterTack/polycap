#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include "hdf5.h"

//===========================================
void polycap_out(size_t n_energies, double *energies, double *efficiencies)
{
	hid_t file, dataset;
	hid_t datatype, dataspace; //handles
	hsize_t n_energies_temp;
	herr_t status;
	char h5filename[] = "polycap_out.h5";


	//Create new HDF5 file using H5F_ACC_TRUNC and default creation and access properties
	file = H5Fcreate(h5filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); 

	//Write energies array
	//Define temporary dataset dimension
	n_energies_temp = n_energies;
	//Describe size of the array and make fixed data space
	dataspace = H5Screate_simple(1, &n_energies_temp, NULL); //RANK = 1 due to 1D dataset
	//Define datatype for the dataspace
	datatype = H5Tcopy(H5T_NATIVE_DOUBLE); //dataset consists of doubles
	status = H5Tset_order(datatype, H5T_ORDER_LE); //Little Endian storage
	//Create new dataset within the HDF5 file with default creation properties
	dataset = H5Dcreate2(file, "Energies", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//Write data to the dataset with default transfer properties
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, energies);
	//Close release sources
	status = H5Sclose(dataspace);
	status = H5Tclose(datatype);
	status = H5Dclose(dataset);

	//Write efficiencies array
	//Define temporary dataset dimension
	n_energies_temp = n_energies;
	//Describe size of the array and make fixed data space
	dataspace = H5Screate_simple(1, &n_energies_temp, NULL); //RANK = 1 due to 1D dataset
	//Define datatype for the dataspace
	datatype = H5Tcopy(H5T_NATIVE_DOUBLE); //dataset consists of doubles
	status = H5Tset_order(datatype, H5T_ORDER_LE); //Little Endian storage
	//Create new dataset within the HDF5 file with default creation properties
	dataset = H5Dcreate2(file, "Transmission Efficiencies", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//Write data to the dataset with default transfer properties
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, efficiencies);
	//Close release sources
	status = H5Sclose(dataspace);
	status = H5Tclose(datatype);
	status = H5Dclose(dataset);

	//Close file
	H5Fclose(file);

	printf("%s was written.\n",h5filename);
}
