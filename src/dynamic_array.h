#ifndef DYNAMIC_ARRAY_H
#define DYNAMIC_ARRAY_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>

/* A struct for dynamic arrays of integers */
struct iDynamicArray{
    int      init_size;
    int      size;
    int*     data;
    int      max_size;
};

typedef struct iDynamicArray iDynamicArray;

/* Creates a DynamicArray */
iDynamicArray idyar_create(const int init_size);


/* Deletes a DynamicArray */
void idyar_delete(iDynamicArray* this);


/* Clears a DynamicArray */
void idyar_clear(iDynamicArray* this);


/* Returns length of array */
int idyar_size(const iDynamicArray* this);


/* Returns the ith element */
int idyar_retrieve(const iDynamicArray* this, const int i);


/* Sets the value of element i */
void idyar_set_val(iDynamicArray* this, const int i, const int val);


/* Adds an element to the end */
void idyar_append(iDynamicArray* this, const int val);


/* Returns a pointer to the underlying data. WARNING: The pointer may change after */
/* every reallocation                                                              */
int* idyar_data(const iDynamicArray* this);


/* Releases additional memory to fit underlying data */
void idyar_shrink_to_fit(iDynamicArray* this);


/* Sorts an array. Order is 'a' for ascending and 'd' for descending */
void idyar_sort(iDynamicArray* this, const char order);


/* Removes all duplicate entries */
void idyar_unique(iDynamicArray* this);


#ifdef __cplusplus
}
#endif

#endif /* DYNAMIC_ARRAY_H */
