#include "io_hdf5.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

//******************************************************************************
//Creates an empty trajectory and opens it
Trajectory traj_create(const char* fn_traj)
{
    Trajectory  traj;
    hid_t       file_id;
    hid_t       group_id;

    //Create trajectory file
    file_id = H5Fcreate(fn_traj, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    //Create "/common" group
    group_id = H5Gcreate(file_id, "/common", H5P_DEFAULT, H5P_DEFAULT,
            H5P_DEFAULT);
    H5Gclose(group_id);

    traj.fn = fn_traj;
    traj.file_id = file_id;
    traj.isopen = true;
    traj.length = 0;
    traj.common_id = -1;
    traj.tvd_id = -1;

    return traj;
}

//******************************************************************************

//Opens a trajectory 
Trajectory traj_open(const char* fn_traj, const char* mode)
{
    Trajectory  traj;
    hid_t       file_id;
    hsize_t     num_obj;
    herr_t      status;
    bool        new_traj;

    new_traj = false;
    
    //Truncate existing file and open in readwrite mode
    if (strcmp(mode, "rw+")==0){
        traj = traj_create(fn_traj);
        new_traj = true;
    }
    //Open in readwrite mode if file exists
    else if (strcmp(mode, "rw")==0){
        //Create new file and open in readwrite mode if file does not exist
        //Function access() defined in unistd.h
        if (access(fn_traj, F_OK) == -1){
            traj = traj_create(fn_traj);
            new_traj = true;
        }
        else {
            file_id = H5Fopen(fn_traj, H5F_ACC_RDWR, H5P_DEFAULT);
        }
    }
    //Open in readonly mode if file exists
    else if (strcmp(mode, "r")==0){
        file_id = H5Fopen(fn_traj, H5F_ACC_RDONLY, H5P_DEFAULT);
    }

    //If an existing file is opened, set values for the struct fields
    if (!new_traj){
        status = H5Gget_num_objs(file_id, &num_obj);

        traj.fn = fn_traj;
        traj.file_id = file_id;
        traj.isopen = true;
        traj.length = num_obj-1;
        traj.common_id = -1;
        traj.tvd_id = -1;
    }

    return traj;
}

//******************************************************************************

//Closes an open trajectory file. Before calling this ensure that no sets are
//open.
void traj_close(Trajectory* traj)
{
    herr_t  status;

    status = H5Fclose(traj->file_id);

    traj->isopen = false;
}

//******************************************************************************

//Create a tvd set
bool traj_create_tvd_set(Trajectory* traj, const double time)
{
    char    name[16];
    int     n;
    hid_t   group_id;
    hid_t   attr_id;
    hid_t   dspace_id;
    herr_t  status;

    n = sprintf(name, "/ts-%li", traj->length);
    group_id = H5Gcreate(traj->file_id, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    //Add time as an attribute
    dspace_id = H5Screate(H5S_SCALAR);
    attr_id = H5Acreate(group_id, "time", H5T_NATIVE_DOUBLE, dspace_id,
            H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &time);
    status = H5Sclose(dspace_id);
    status = H5Aclose(attr_id);

    //Close group
    status = H5Gclose(group_id);

    //Increment length
    traj->length += 1;

    return (status >= 0);
}

//******************************************************************************

//Open a set
bool traj_open_set(Trajectory* traj, const char* set_type, const long int index)
{
    char    name[16];
    int     n;
    bool    ret;

    if (strcmp(set_type, "common")==0){
        traj->common_id = H5Gopen(traj->file_id, "/common", H5P_DEFAULT);
        if (traj->common_id < 0) ret = false;
    }
    else if (strcmp(set_type, "tvd")==0){
        n = sprintf(name, "/ts-%li", index);
        traj->tvd_id = H5Gopen(traj->file_id, name, H5P_DEFAULT);
        if (traj->tvd_id < 0) ret = false;
    }
    else {
        ret = false;
    }

    return ret;
}

//******************************************************************************

//Close a set
void traj_close_set(Trajectory* traj, const char* set_type)
{
    herr_t  status;

    if (strcmp(set_type, "common")==0){
        status = H5Gclose(traj->common_id);
        traj->common_id = -1;
    }
    else if (strcmp(set_type, "tvd")==0){
        status = H5Gclose(traj->tvd_id);
        traj->tvd_id = -1;
    }
}

//******************************************************************************

//Write data to an already open set
void traj_write_to_set(Trajectory* traj, const char* set_type, const char* data_name,
        const char* data_type, const int data_rank, const int* data_extent,
        const void* data_buffer)
{
    herr_t      status;
    hid_t       group_id;
    hid_t       dtype_id;
    hid_t       dspace_id;
    hid_t       dset_id;
    hsize_t     extent[data_rank]; //Type must be hsize_t, won't work if declared as int

    for (int i=0; i<data_rank; ++i){
        extent[i] = data_extent[i];
    }

    if (strcmp(set_type, "common")==0){
        group_id = traj->common_id;
    }
    else if (strcmp(set_type, "tvd")==0){
        group_id = traj->tvd_id;
    }

    //If dataset already exists
    if (H5Lexists(group_id, data_name, H5P_DEFAULT) > 0){
        printf("dataset exits\n");
        return ;
    }

    //Create dataspace
    if (data_rank==0){
        dspace_id = H5Screate(H5S_SCALAR); 
    }
    else {
        dspace_id = H5Screate_simple(data_rank, extent, NULL);
    }

    //Separate calls based on data type. Note that the predefined HDF5 types 
    //are macro definitions.
    if (strcmp(data_type, "int")==0){
        //Create dataset in group with identifier `group_id`
        dset_id = H5Dcreate(group_id, data_name, H5T_NATIVE_INT, dspace_id,
                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        //Write dataset
        status = H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                data_buffer);
    }
    else if (strcmp(data_type, "double")==0){
        //Create dataset in group with identifier `group_id`
        dset_id = H5Dcreate(group_id, data_name, H5T_NATIVE_DOUBLE, dspace_id,
                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        //Write dataset
        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                H5P_DEFAULT, data_buffer);
    }

    //Close dataspace
    status = H5Sclose(dspace_id);

    //Close datasset
    status = H5Dclose(dset_id);
}

//******************************************************************************

//Read data from an already open set
void traj_read_from_set(Trajectory* traj, const char* set_type, 
        const char* data_name, const char* data_type, void* data_buffer)
{
    herr_t  status;
    hid_t   group_id;
    hid_t   dtype_id;
    hid_t   dset_id;

    if (strcmp(set_type, "common")==0){
        group_id = traj->common_id;
    }
    else if (strcmp(set_type, "tvd")==0){
        group_id = traj->tvd_id;
    }

    //Open dataset in group with identifier `group_id`
    dset_id = H5Dopen(group_id, data_name, H5P_DEFAULT);

    //Separate calls to read dataset based on data type
    if (strcmp(data_type, "int")==0){
        status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                data_buffer);
    }
    else if (strcmp(data_type, "double")==0){
        status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                data_buffer);
    }

    //Close dataset
    status = H5Dclose(dset_id);
}

//******************************************************************************

//Read time from a tvd_set specified by index `index`

double traj_get_time(Trajectory* traj, const long int index)
{
    hid_t   attr_id;
    hid_t   group_id;
    herr_t  status;
    double  time;
    char    name[16];

    if (index >= traj->length){
        time = nan("");
    }
    else{
        //Open tvd_set
        sprintf(name, "/ts-%li", index);
        group_id = H5Gopen(traj->file_id, name, H5P_DEFAULT);

        //Open attribute
        attr_id = H5Aopen(group_id, "time", H5P_DEFAULT);

        //Read attribute
        status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &time);

        //Close attribute
        status = H5Aclose(attr_id);

        //Close tvd set
        status = H5Gclose(group_id);
    }

    return time;
}
//******************************************************************************
