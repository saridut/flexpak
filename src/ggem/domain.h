#ifndef DOMAIN_H
#define DOMAIN_H

#ifdef __cplusplus
extern "C" {
#endif

struct Domain{
    double lx;
    double ly;
    double lz;
    int [3] pbc;
};

typedef struct Domain Domain;

#ifdef __cplusplus
}
#endif

#endif /* DOMAIN_H */
