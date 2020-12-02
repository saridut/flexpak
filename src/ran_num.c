#include "ran_num.h"
#include <stdio.h>
#include <string.h>

//******************************************************************************

//Returns an unsigned int for use as seed from a file named fn.
//If fn="", read from /dev/urandom.
unsigned int get_random_seed(const char* fn){

    unsigned int seed;
    FILE*       fp;

    if (strcmp(fn, "")==0){
        fp = fopen("/dev/urandom", "r");
        fread(&seed, sizeof(seed), 1, fp);
        fclose(fp);
    }
    else {
        fp = fopen(fn, "r");
        fscanf(fp, "%u", &seed);
        fclose(fp);
    }
    return seed;
}

//******************************************************************************

RandomStream init_stream(const unsigned int seed){

    RandomStream    stream;
    
    //Using Mersenne Twister 19937
    vslNewStream(&stream, VSL_BRNG_MT19937, seed);
    return stream;
}

//******************************************************************************

//Deletes a stream
void delete_stream(RandomStream* stream){
    vslDeleteStream(stream);
}

//******************************************************************************

//Loads a stream from a binary file
RandomStream load_stream(const char* fn_stream){

    RandomStream    stream;

    vslLoadStreamF(&stream, fn_stream);
    return stream;
}

//******************************************************************************

//Saves a stream to a binary file
void save_stream(const RandomStream stream, const char* fn_stream){

    vslSaveStreamF(stream, fn_stream);
}
//******************************************************************************

//Generates n uniformly distributed integers over [lb, ub)
void iuni(RandomStream stream, const int lb, const int ub,
        const int n, int* u){

    viRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, n, u, lb, ub);
}

//******************************************************************************

//Generates n uniformly distributed doubles over [lb, ub)
void duni(RandomStream stream, const double lb, const double ub,
        const int n, double* u){

    vdRngUniform( VSL_RNG_METHOD_UNIFORM_STD_ACCURATE, stream, n, u, lb, ub );
}

//******************************************************************************
//Generates n noramlly distributed doubles 
void dnor(RandomStream stream, const double mean, const double std_dev,
        const int n, double* u){

    vdRngGaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, n, u, mean, std_dev);
}

//******************************************************************************
