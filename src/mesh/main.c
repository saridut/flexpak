#include <stdio.h>
#include "dynamic_array.h"

int main()
{
    iDynamicArray A = idyar_create(2);
    printf("Size = %d\n", idyar_size(&A));

    for (int i=0; i <10; ++i){
        idyar_append(&A, 2*i);
    }

    printf("Size = %d\n", idyar_size(&A));

    for (int i=0; i <10; ++i){
        printf("%d\n", idyar_retrieve(&A, i));
    }

    idyar_clear(&A);
    printf("Size = %d\n", idyar_size(&A));

    idyar_append(&A, 2);
    idyar_append(&A, 2);
    idyar_append(&A, -1);
    idyar_append(&A, 10);
    idyar_append(&A, 10);
    idyar_append(&A, 10);
    idyar_append(&A, 0);
    idyar_append(&A, -4);
    idyar_append(&A, -4);

    idyar_sort(&A, 'a');

    printf("\n");

    for (int i=0; i < idyar_size(&A); ++i){
        printf("%d  %d\n", i, idyar_retrieve(&A, i));
    }

    idyar_unique(&A);

    printf("\n");
    for (int i=0; i < idyar_size(&A); ++i){
        printf("%d  %d\n", i, idyar_retrieve(&A, i));
    }

    idyar_delete(&A);
}
