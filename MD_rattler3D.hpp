#ifndef _MDRAT3_HPP_
#define _MDRAT3_HPP_

#include <vector>

void MD_rattler3D_cleanByCn(const int N, const double * Dn, const double * x, bool * is_rattler, int * Cn, const double gam);

std::vector<double> MD_rattler3D_getNeighbors(const int N, const double * Dn, const double * x, bool * is_rattler, const double gam, const int nn);

double euclid3_norm(double * n);

double get_acos3(double * v1, double * v2);

bool MD_rattler3D_check(std::vector<double> neighbors);

void MD_getRattler3D(const int N, const double * Dn, const double * x, bool * is_rattler, const double gam);

#endif