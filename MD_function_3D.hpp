#ifndef _MDFUNC3_HPP_
#define _MDFUNC3_HPP_

#include <iostream>
#include <vector>
#include <tuple>

void logspace(int low, int high, int intervals, double * list);

void MD_getRattler(int N, double * Dn, double * x, bool * is_rattler, double gam);

void MD_getStressTensor(const int &N, double * Dn, double * x, double * stress, const double &alpha, const double &gam);

double MD_getP(int N, double * Dn, double * x, double alpha, double gam);

void MD_getCn(const int &N, const double * Dn, const double * x, const bool * is_rattler, int * Cn, const double &gam);

void MD_getCtcNwk(const int &N, double * Dn, double * x, bool * contact_network, const double &gam);

void unjamToJustTouch(int &N, double * Dn, double * m, double * x, double gam);

double MD_getP_DS(const int &N, double * Dn, double * x, const bool * contact_network, const double &alpha, const double &gam);

double MD_getFtol(const int &N, double * Dn, double * x, double * Fx, const double &alpha, const double &gam);

bool ifUpdateNL(const int &N, double * Dn, double * x, double * Dn_skin, double * x_skin, const double &gam);

double buildNeighborList(const int &N, double * Dn, double * x, double * Dn_skin, double * x_skin, double * Fx, const double &alpha, const double &gam, const double &skin_scale, bool * NL);

double getFtolFromNL(const int &N, double * Dn, double * x, double * Fx, const double &alpha, const double &gam, const bool * NL);

double FIRE(int N, double * vx, double * Fx, double a);

double FIRE_getFv(const int &N, double * vx, double * Fx);

void FIRE_mixing(const int &N, double * vx, double * Fx, const double &a);

void MD_CMzeroing(int N, double * vx, double * m);

double MD_cmpjam_minimization(int N, double * Dn, double * m, double * x, double dt, double alpha, double dphi, double gam, const double& tol);

double MD_cmpjam_minimization_VV_FIRE(int N, double * Dn, double * m, double * x, double dt, double alpha, double dphi, double gam, const double& tol);

double MD_shrjam_minimization(const int &N, double * Dn, double * m, double * x, double dt, const double &alpha, const double &dG, double &gam, const double& tol);

void MD_cmpjam_main(int N, double * Dn, double * m, double * x, double dt, double alpha, double Ptol, double dphi, double gam);

void MD_cmpjam_main_DS(int N, double * Dn, double * m, double * x, bool * contact_network, double dt, double alpha, double Ptol, double dphi, double gam);

void scale(int N, double * Dn, double * m, double rescale);

bool anytouch(int N, double * pos, double * sc);

bool fileExist(const char * name);

bool IsPathExist(const char * s);

void writeResult(const char * path, int N, double * Dn, double * m, double * x);

void writeResultToFile(FILE * xyz, const int &N, double * Dn, double * m, double * x);

void loadResult(const char * path, int N, double * Dn, double * m, double * x);

void MD_shearModulus_main(const int &N, const double &alpha, const char * loadpath, const char * savepath, const char * saveCPpath, const double &Ptol, const bool &getCPstate, const bool &positive_shear, const bool &HQ_flag);

void MD_mapping_shear_func(const int &N, const double &alpha, const char * loadpath, double * Plist, const double &dG);

void MD_mapping_shearCP_func(const int &N, const double &alpha, const char * loadpath, const char * savepath, const double &Ptol, const double &dG);

int getNumberofLines(const char * path);

void loadResultwLineNumber(const char * path, int N, double * Dn, double * m, double * x, const int &target_lines);

void MD_mapping_CPfamily_func(const int &N, const double &alpha, const char * loadpath, const char * savepath, const double &Ptol, const double &dG, const bool &early_termination);

#endif