#include "dynamic_array.h"
#include <string.h>

/******************************************************************************/

/* Creates a DynamicArray initialized to zero */
iDynamicArray idyar_create(const int init_size)
{
    iDynamicArray ida;

    if (init_size <= 0){
        ida.init_size = 8;
    }
    else {
        ida.init_size = init_size;
    }
    ida.max_size = ida.init_size;
    ida.size = 0;
    ida.data = calloc(ida.max_size, sizeof(int));

    return ida;
}

/******************************************************************************/

/* Deletes a DynamicArray*/
void idyar_delete(iDynamicArray* this)
{
    free(this->data);
    this->size = 0;
}

/******************************************************************************/

/* Clears a DynamicArray -- returns to a state same as that after a call to
 * idyar_create*/
void idyar_clear(iDynamicArray* this)
{
    this->max_size = this->init_size;
    this->size = 0;
    this->data = realloc(this->data, this->max_size*sizeof(int));
    for(int i=0; i < this->max_size; ++i){
        this->data[i] = 0;
    }
}

/******************************************************************************/

/* Returns length of array */
int idyar_size(const iDynamicArray* this)
{
    return this->size;
}

/******************************************************************************/

/* Returns the ith element */
int idyar_retrieve(const iDynamicArray* this, const int i)
{
    return this->data[i];
}

/******************************************************************************/

/* Sets the value of element i */
void idyar_set_val(iDynamicArray* this, const int i, const int val)
{
    this->data[i] = val;
}

/******************************************************************************/

/* Adds an element to the end */
void idyar_append(iDynamicArray* this, const int val)
{
    /* Double max_size if current size is equal to max_size */
    if (this->size == this->max_size){
        this->max_size *= 2;
        this->data = realloc(this->data, this->max_size*sizeof(int));
    }
    this->data[this->size] = val;
    this->size += 1;
}

/******************************************************************************/

/* Returns a pointer to the underlying data.                 */
/* WARNING: The pointer may change after every reallocation. */
int* idyar_data(const iDynamicArray* this)
{
    return this->data;
}

/******************************************************************************/

/* Releases additional memory to fit underlying data */
void idyar_shrink_to_fit(iDynamicArray* this)
{
    if (this->size % 2 == 0){
        this->max_size = this->size;
    }
    else {
        this->max_size = this->size + 1;
    }
    this->data = realloc(this->data, this->max_size*sizeof(int));
}

/******************************************************************************/

/* Sorts an array. Order is 'a' for ascending and 'd' for descending */

/* Comparison function for ascending order */
int cmp_is_lesser(const void* a, const void* b)
    /* If return value is negative, *a goes before *b */
{
    return *(int*)a - *(int*)b;
}

/* Comparison function for descending order */
int cmp_is_greater(const void* a, const void* b)
    /* If return value is negative, *a goes before *b */
{
    return *(int*)b - *(int*)a;
}

/* Qsort */
void idyar_sort(iDynamicArray* this, const char order)
{
    /* Sort in ascending order */
    if (order == 'a'){
        qsort(this->data, this->size, sizeof(int), cmp_is_lesser);
    }

    /* Sort in descending order */
    if (order == 'd'){
        qsort(this->data, this->size, sizeof(int), cmp_is_greater);
    }
}

/******************************************************************************/

/* Removes all duplicate entries */
void idyar_unique(iDynamicArray* this)
{
    int* unique_data = malloc(this->size*sizeof(int));

    /* The first element is always unique */
    unique_data[0] = this->data[0];
    int n = 0; /* Index for unique elements */

    for (int i=1; i < this->size; ++i){
        if (this->data[i] != this->data[i-1]){
            n += 1;
            unique_data[n] = this->data[i];
        }
    }
    memcpy(this->data, unique_data, (n+1)*sizeof(int));

    this->size = n + 1;

    free(unique_data);

    idyar_shrink_to_fit(this);
}

/******************************************************************************/
