#ifndef RAN_NUM_H
#define RAN_NUM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <mkl_vsl.h>

typedef VSLStreamStatePtr RandomStream;

unsigned int get_random_seed();

RandomStream init_stream(const unsigned int seed);

void delete_stream(RandomStream* stream);

RandomStream load_stream(const char* fn_stream );

void save_stream(const RandomStream stream, const char* fn_stream);

void iuni(RandomStream stream, const int lb, const int ub, const int n, int* u);

void duni(RandomStream stream, const double lb, const double ub, const int n, double* u);

void dnor(RandomStream stream, const double mean, const double std_dev,
        const int n, double* u);

#ifdef __cplusplus
}
#endif

#endif // RAN_NUM_H
