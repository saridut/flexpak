#ifndef IO_HDF5_H
#define IO_HDF5_H

#ifdef __cplusplus
extern "C" {
#endif

#include <hdf5.h>
#include <stdbool.h>

struct trajectory {   
    const char* fn;
    hid_t       file_id;
    bool        isopen;
    long int    length;
    hid_t       common_id;
    hid_t       tvd_id;
} ;

typedef struct trajectory Trajectory;

//Create a new trajectory
Trajectory traj_create(const char* fn_traj);

//Open an existing trajectory
Trajectory traj_open(const char* fn_traj, const char* mode);

//Close an open trajectory
void traj_close(Trajectory* traj);

//Create a tvd set in an open trajectory
bool traj_create_tvd_set(Trajectory* traj, const double time);

//Open a set in an open trajectory
bool traj_open_set(Trajectory* traj, const char* set_type, const long int index);

//Close a set in an open trajectory
void traj_close_set(Trajectory* traj, const char* set_type);

//Write to an open set in an open trajectory
void traj_write_to_set(Trajectory* traj, const char* set_type,
                        const char* data_name, const char* data_type,
                        const int data_rank, const int* data_extent,
                        const void* data_buffer);

//Read from an open set in an open trajectory
void traj_read_from_set(Trajectory* traj, const char* set_type, 
                        const char* data_name, const char* data_type,
                        void* data_buffer);

//Get the value of time from an open tvd set in an open trajectory
double traj_get_time(Trajectory* traj, const long int index);

#ifdef __cplusplus
}
#endif

#endif // IO_HDF5_H
