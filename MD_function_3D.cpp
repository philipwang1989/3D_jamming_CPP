#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <math.h>
#include <limits>
#include <algorithm>
#include <vector>
#include <tuple>
#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <ctime>
#include <sys/stat.h>
#include "MD_rattler3D.hpp"

using namespace std;

void logspace(int low, int high, int intervals, double * list)
{
    double delta = double(high - low) / double(intervals - 1);
    for (int i=0; i<intervals; i++)
    {
        double curr = low + i*delta;
        list[i] = pow(10, curr);
    }
}

void MD_getCn(const int &N, const double * Dn, const double * x, const bool * is_rattler, int * Cn, const double &gam)
{
    double dx, dy, dz, im, dnm, Dnm;
    for (int i=0; i<N; i++) Cn[i] = 0;
    for (int nn=0; nn<N; nn++)
    {
        for (int mm=nn+1; mm<N; mm++)
        {
            if (!is_rattler[nn] && !is_rattler[mm])
            {
                dy = x[3*mm+1]-x[3*nn+1];
                im = round(dy);
                dy -= im;
                dx = x[3*mm]-x[3*nn];
                // dx -= round(dx);
                dx -= (round(dx-im*gam)+im*gam);
                dz = x[3*mm+2]-x[3*nn+2];
                dz -= round(dz);
                Dnm = 0.5 * (Dn[nn]+Dn[mm]);
                if (fabs(dy)<Dnm)
                {
                    dnm = sqrt(dx*dx+dy*dy+dz*dz);
                    if (dnm < Dnm)
                    {
                        Cn[nn] += 1;
                        Cn[mm] += 1;
                    }
                }
            }
        }
    }
}

void MD_getStressTensor(const int &N, double * Dn, double * x, double * stress, const double &alpha, const double &gam)
{
    // x runs x0,y0,x1,y1...xN,yN
    int nn, mm;
    double P, dx, dy, dz, Dnm, dnm, F, im, D;
    // double stress[4] = {0}; // Sxx,Sxy,Syx,Syy
    for (int i=0; i<9; i++) stress[i] = 0.0;
    bool * is_rattler = new bool[N];
    for (int i=0; i<N; i++) is_rattler[i] = false;
    // MD_getRattler(N, Dn, x, is_rattler, gam);
    D = 1.0;
    for (nn=0; nn<N; nn++)
    {
        D = min(D, Dn[nn]);
        for (mm=nn+1; mm<N; mm++)
        {
            if (!is_rattler[nn] && !is_rattler[mm])
            {
                dy = x[3*mm+1]-x[3*nn+1];
                im = round(dy);
                dy -= im;
                dx = x[3*mm]-x[3*nn];
                // dx -= round(dx);
                dx -= (round(dx-im*gam)+im*gam);
                dz = x[3*mm+2]-x[3*nn+2];
                dz -= round(dz);
                Dnm = 0.5 * (Dn[nn]+Dn[mm]);
                if (fabs(dy)<Dnm)
                {
                    dnm = sqrt(dx*dx+dy*dy+dz*dz);
                    if (dnm < Dnm)
                    {
                        F = -pow((1-dnm/Dnm),alpha-1.0)/Dnm/dnm;
                        stress[0] -= F*dx*dx;             // Sxx
                        stress[1] -= 0.5*F*(dx*dy+dy*dx); // Sxy
                        stress[2] -= 0.5*F*(dz*dx+dx*dz); // Sxz
                        stress[3] -= 0.5*F*(dy*dx+dx*dy); // Syx
                        stress[4] -= F*dy*dy;             // Syy
                        stress[5] -= 0.5*F*(dy*dz+dz*dy); // Syz
                        stress[6] -= 0.5*F*(dz*dx+dx*dz); // Szx
                        stress[7] -= 0.5*F*(dz*dy+dy*dz); // Szy
                        stress[8] -= F*dz*dz;             // Szz
                    }
                }
            }
            
        }
    }
    for (int i=0; i<9; i++) stress[i] *= (D*D*D/8.0);
}

double MD_getP(int N, double * Dn, double * x, double alpha, double gam)
{
    // x runs x0,y0,x1,y1...xN,yN
    int nn, mm;
    double P, dx, dy, dz, Dnm, dnm, F, im, D;
    double stress[9] = {0.0}; // Sxx,Sxy,Syx,Syy
    bool * is_rattler = new bool[N];
    for (int i=0; i<N; i++) is_rattler[i] = false;
    MD_getRattler3D(N, Dn, x, is_rattler, gam);
    D = 1.0;
    for (nn=0; nn<N; nn++)
    {
        D = min(D, Dn[nn]);
        for (mm=nn+1; mm<N; mm++)
        {
            if (!is_rattler[nn] && !is_rattler[mm])
            {
                dy = x[3*mm+1]-x[3*nn+1];
                im = round(dy);
                dy -= im;
                dx = x[3*mm]-x[3*nn];
                // dx -= round(dx);
                dx -= (round(dx-im*gam)+im*gam);
                dz = x[3*mm+2]-x[3*nn+2];
                dz -= round(dz);
                Dnm = 0.5 * (Dn[nn]+Dn[mm]);
                // if (fabs(dy)<Dnm)
                {
                    dnm = sqrt(dx*dx+dy*dy+dz*dz);
                    if (dnm < Dnm)
                    {
                        F = -pow((1-dnm/Dnm),alpha-1.0)/Dnm/dnm;
                        stress[0] -= F*dx*dx;             // Sxx
                        stress[1] -= 0.5*F*(dx*dy+dy*dx); // Sxy
                        stress[2] -= 0.5*F*(dz*dx+dx*dz); // Sxz
                        stress[3] -= 0.5*F*(dy*dx+dx*dy); // Syx
                        stress[4] -= F*dy*dy;             // Syy
                        stress[5] -= 0.5*F*(dy*dz+dz*dy); // Syz
                        stress[6] -= 0.5*F*(dz*dx+dx*dz); // Szx
                        stress[7] -= 0.5*F*(dz*dy+dy*dz); // Szy
                        stress[8] -= F*dz*dz;             // Szz
                    }
                }
            }
            
        }
    }
    for (int i=0; i<9; i++) stress[i] *= (D*D*D/8.0);
    P = (stress[0] + stress[4] + stress[8])/3.0;
    return P;
}

// STILL 2D
void MD_getCtcNwk(const int &N, double * Dn, double * x, bool * contact_network, const double &gam)
{
    int nn, mm;
    double dx, dy, Dnm, dnm, F, im;
    bool * is_rattler = new bool[N];
    for (int i=0; i<N; i++) is_rattler[i] = false;
    MD_getRattler3D(N, Dn, x, is_rattler, gam);
    for (int i=0; i<N*N; i++) contact_network[i] = false;
    for (nn=0; nn<N; nn++)
    {
        for (mm=nn+1; mm<N; mm++)
        {
            if (!is_rattler[nn] && !is_rattler[mm])
            {
                dy = x[2*mm+1]-x[2*nn+1];
                im = round(dy);
                dy -= im;
                dx = x[2*mm]-x[2*nn];
                dx -= round(dx);
                dx -= (round(dx-im*gam)+im*gam);
                Dnm = 0.5 * (Dn[nn]+Dn[mm]);
                if (fabs(dy)<Dnm)
                {
                    dnm = sqrt(dx*dx+dy*dy);
                    if (dnm < Dnm) contact_network[nn*N+mm] = true;
                }
            }
        }
    }
}

// STILL 2D
void unjamToJustTouch(int &N, double * Dn, double * m, double * x, double gam)
{
    // step 1, find minimum overlap
    int nn, mm;
    double dx, dy, im, Dnm, dnm;
    double maxOverlap = 0.0;
    double rsc;
    for (nn=0; nn<N; nn++)
    {
        for (mm=nn+1; mm<N; mm++)
        {
            dy = x[2 * mm + 1] - x[2 * nn + 1];
            im = round(dy);
            dy -= im;
            dx = x[2 * mm] - x[2 * nn];
            // dx -= round(dx);
            dx -= (round(dx - im * gam) + im * gam);
            Dnm = 0.5 * (Dn[nn] + Dn[mm]);
            if (fabs(dy) < Dnm)
            {
                dnm = sqrt(dx * dx + dy * dy);
                if (dnm < Dnm && Dnm-dnm > maxOverlap) 
                {
                    maxOverlap = Dnm-dnm;
                    rsc = (Dnm - maxOverlap)/Dnm;
                } 
            }
        }
    }
    // step 2, unjam by rsc
    for (int i=0; i<N; i++)
    {
        Dn[i] = Dn[i] * rsc;
        m[i] = m[i] * rsc * rsc;
    }
}

// STILL 2D
double MD_getP_DS(const int &N, double * Dn, double * x, const bool * contact_network, const double &alpha, const double &gam)
{
    // x runs x0,y0,x1,y1...xN,yN
    int nn, mm;
    double P, dx, dy, Dnm, dnm, F, im, D;
    double stress[4] = {0}; // Sxx,Sxy,Syx,Syy
    D = 1.0;
    for (nn=0; nn<N; nn++)
    {
        D = min(D, Dn[nn]);
        for (mm=nn+1; mm<N; mm++)
        {
            if (contact_network[nn*N+mm])
            {
                dy = x[2*mm+1]-x[2*nn+1];
                im = round(dy);
                dy -= im;
                dx = x[2*mm]-x[2*nn];
                // dx -= round(dx);
                dx -= (round(dx-im*gam)+im*gam);
                Dnm = 0.5 * (Dn[nn]+Dn[mm]);
                dnm = sqrt(dx*dx+dy*dy);
                double delta = 1.0-dnm/Dnm;
                if (delta < 0) // dnm > Dnm - raises ERROR in std::pow
                {
                    F = pow(-delta,alpha-1.0)/Dnm/dnm; // becomes attraction!
                }
                else F = -pow(delta,alpha-1.0)/Dnm/dnm;
                
                if (isnan(F))
                {
                    cout << "NAN F!" << endl;
                    cout << pow((1-dnm/Dnm),alpha-1.0) << endl;
                    cout << (1-dnm/Dnm) << endl;
                    cout << alpha-1.0 << endl;
                } 
                stress[0] -= F*dx*dx;
                stress[1] -= 0.5*F*(dx*dy+dy*dx);
                stress[2] -= 0.5*F*(dy*dx+dx*dy);
                stress[3] -= F*dy*dy;
            }
            
        }
    }
    for (int i=0; i<4; i++) stress[i] *= (D*D/4);
    P = (stress[0] + stress[3])/2;
    return P;
}

double MD_getFtol(const int &N, double * Dn, double * x, double * Fx, const double &alpha, const double &gam)
{
    for (int i=0; i<3*N; i++) Fx[i] = 0.0;
    int nn, mm;
    double Ftol, dx, dy, dz, Dnm, dnm, F, im, D;
    Ftol = 0;
    for (nn=0; nn<N; nn++)
    {
        D = min(D, Dn[nn]);
        for (mm=nn+1; mm<N; mm++)
        {
            dy = x[3*mm+1]-x[3*nn+1];
            im = round(dy);
            dy -= im;
            dx = x[3*mm] - x[3*nn];
            // dx -= round(dx);
            dx -= (round(dx-im*gam)+im*gam);
            dz = x[3*mm+2]-x[3*nn+2];
            dz -= round(dz);
            Dnm = 0.5 * (Dn[nn] + Dn[mm]);
            dnm = sqrt(dx * dx + dy * dy + dz * dz);
            if (dnm < Dnm)
            {
                F = -pow((1 - dnm / Dnm), alpha - 1.0) / Dnm / dnm;
                Fx[3*nn] += F*dx;
                Fx[3*mm] -= F*dx;
                Fx[3*nn+1] += F*dy;
                Fx[3*mm+1] -= F*dy;
                Fx[3*nn+2] += F*dz;
                Fx[3*mm+2] -= F*dz;
            }
        }
    }    
    for (int i=0; i<3*N; i++) Ftol += (Fx[i] * Fx[i]);
    Ftol /= (N*N);
    return Ftol;
}

bool ifUpdateNL(const int &N, double * Dn, double * x, double * Dn_skin, double * x_skin, const double &gam)
{
    double dx, dy, dz, dnm, im;
    for (int nn=0; nn<N; nn++)
    {
        dy = x[3 * nn + 1] - x_skin[3 * nn + 1];
        im = round(dy);
        dy -= im;
        dx = x[3 * nn] - x_skin[3 * nn];
        // dx -= round(dx);
        dx -= (round(dx - im * gam) + im * gam);
        dz = x[3 * nn + 2] - x_skin[3 * nn + 2];
        dnm = sqrt(dx * dx + dy * dy + dz * dz);
        if (2*dnm+Dn[nn] > Dn_skin[nn]) return true;
    }
    return false;
}

double buildNeighborList(const int &N, double * Dn, double * x, double * Dn_skin, double * x_skin, double * Fx, const double &alpha, const double &gam, const double &skin_scale, bool * NL)
{
    for (int i=0; i<3*N; i++) Fx[i] = 0.0;
    for (long i=0; i<N*N; i++) NL[i] = false;
    long nn, mm;
    double Ftol, dx, dy, dz, Dnm, dnm, F, im, D;
    Ftol = 0;
    for (nn=0; nn<N; nn++)
    {
        x_skin[nn] = x[nn];
        Dn_skin[nn] = Dn[nn] * skin_scale;
        D = min(D, Dn[nn]);
        for (mm=nn+1; mm<N; mm++)
        {
            dy = x[3*mm+1] - x[3*nn+1];
            im = round(dy);
            dy -= im;
            dx = x[3*mm] - x[3*nn];
            // dx -= round(dx);
            dx -= (round(dx - im * gam) + im * gam);
            dz = x[3*mm+2] - x[3*nn+2];
            Dnm = 0.5 * (Dn[nn] + Dn[mm]);
            if (fabs(dy) < Dnm)
            {
                dnm = sqrt(dx*dx + dy*dy + dz*dz);
                if (dnm < Dnm)
                {
                    F = -pow((1 - dnm / Dnm), alpha - 1.0) / Dnm / dnm;
                    Fx[3*nn] += F*dx;
                    Fx[3*mm] -= F*dx;
                    Fx[3*nn+1] += F*dy;
                    Fx[3*mm+1] -= F*dy;
                    Fx[3*nn+2] += F*dz;
                    Fx[3*mm+2] -= F*dz;
                }
                if (dnm < Dnm*skin_scale) NL[nn*N+mm] = true;
            }
        }
    }    
    for (int i=0; i<3*N; i++) Ftol += (Fx[i] * Fx[i]);
    Ftol /= (N*N);
    return Ftol;
}

double getFtolFromNL(const int &N, double * Dn, double * x, double * Fx, const double &alpha, const double &gam, const bool * NL)
{
    for (int i=0; i<3*N; i++) Fx[i] = 0.0;
    long nn, mm;
    double Ftol, dx, dy, dz, Dnm, dnm, F, im, D;
    Ftol = 0;
    for (nn=0; nn<N; nn++)
    {
        D = min(D, Dn[nn]);
        for (mm=nn+1; mm<N; mm++)
        {
            if (NL[nn*N+mm])
            {
                dy = x[3*mm+1] - x[3*nn+1];
                im = round(dy);
                dy -= im;
                dx = x[3*mm] - x[3*nn];
                // dx -= round(dx);
                dx -= (round(dx - im * gam) + im * gam);
                dz = x[3*mm+2] - x[3*nn+2];
                dz -= round(dz);
                Dnm = 0.5 * (Dn[nn] + Dn[mm]);
                if (fabs(dy) < Dnm)
                {
                    dnm = sqrt(dx * dx + dy * dy);
                    if (dnm < Dnm)
                    {
                        F = -pow((1 - dnm / Dnm), alpha - 1.0) / Dnm / dnm;
                        Fx[3*nn] += F*dx;
                        Fx[3*mm] -= F*dx;
                        Fx[3*nn+1] += F*dy;
                        Fx[3*mm+1] -= F*dy;
                        Fx[3*nn+2] += F*dz;
                        Fx[3*mm+2] -= F*dz;
                    }
                }
            }
        }
    }    

    for (int i=0; i<3*N; i++) Ftol += (Fx[i] * Fx[i]);
    Ftol /= (N*N);
    return Ftol;
}

double FIRE(int N, double * vx, double * Fx, double a)
{
    // norm - sqrt(sum of the squared)
    double normFx = 0;
    double normFy = 0;
    double normFz = 0;
    double normvx = 0;
    double normvy = 0;
    double normvz = 0;
    double P = 0;
    for (int i=0; i<3*N; i++) P += vx[i]*Fx[i];

    for (int i=0; i<N; i++)
    {
        normFx += (Fx[3*i]*Fx[3*i]);
        normFy += (Fx[3*i+1]*Fx[3*i+1]);
        normFz += (Fx[3*i+2]*Fx[3*i+2]);
        normvx += (vx[3*i]*vx[3*i]);
        normvy += (vx[3*i+1]*vx[3*i+1]);
        normvy += (vx[3*i+2]*vx[3*i+2]);
    }
    normFx = sqrt(normFx);
    normFy = sqrt(normFy);
    normFz = sqrt(normFz);
    normvx = sqrt(normvx);
    normvy = sqrt(normvy);
    normvz = sqrt(normvz);

    for (int i=0; i<N; i++)
    {
        vx[3*i] = (1-a) * vx[3*i] + a * (Fx[3*i]/normFx) * normvx;
        vx[3*i+1] = (1-a) * vx[3*i+1] + a * (Fx[3*i+1]/normFy) * normvy;
        vx[3*i+2] = (1-a) * vx[3*i+2] + a * (Fx[3*i+2]/normFz) * normvz;
    }

    return P;
}

double FIRE_getFv(const int &N, double * vx, double * Fx)
{
    // norm - sqrt(sum of the squared)
    double P = 0;
    for (int i=0; i<3*N; i++) P += vx[i]*Fx[i];
    return P;
}

// PENDING 3D
void FIRE_mixing(const int &N, double * vx, double * Fx, const double &a)
{
    double normFx = 0;
    double normFy = 0;
    double normvx = 0;
    double normvy = 0;
    double normFz = 0;
    double normvz = 0;
    for (int i=0; i<N; i++)
    {
        normFx += (Fx[2*i]*Fx[2*i]);
        normFy += (Fx[2*i+1]*Fx[2*i+1]);
        normvx += (vx[2*i]*vx[2*i]);
        normvy += (vx[2*i+1]*vx[2*i+1]);
    }
    normFx = sqrt(normFx);
    normFy = sqrt(normFy);
    normvx = sqrt(normvx);
    normvy = sqrt(normvy);

    for (int i=0; i<N; i++)
    {
        vx[2*i] = (1-a) * vx[2*i] + a * (Fx[2*i]/normFx) * normvx;
        vx[2*i+1] = (1-a) * vx[2*i+1] + a * (Fx[2*i+1]/normFy) * normvy;
    }
}

void MD_CMzeroing(int N, double * vx, double * m)
{
    double vx_cm, vy_cm, vz_cm, sum_m, Px, Py, Pz;
    sum_m = 0.0;
    Px = 0.0;
    Py = 0.0;
    Pz = 0.0;
    for (int i=0; i<N; i++)
    {
        Px += m[i] * vx[3*i];
        Py += m[i] * vx[3*i+1];
        Pz += m[i] * vx[3*i+2];
        sum_m += m[i];
    }
    vx_cm = Px / sum_m;
    vy_cm = Py / sum_m;
    vz_cm = Pz / sum_m;
    for (int i=0; i<N; i++)
    {
        vx[3*i] = vx[3*i] - vx_cm;
        vx[3*i+1] = vx[3*i+1] - vy_cm;
        vx[3*i+2] = vx[3*i+2] - vz_cm;
    }
}

double MD_cmpjam_minimization(int N, double * Dn, double * m, double * x, double dt, double alpha, double dphi, double gam, const double& tol)
{
    // first grow particles by dphi
    double totalarea = 0.0;
    for(int i = 0; i < N; i++) totalarea += m[i];
    if (fabs(dphi) > 1e-16)
    {
        double rsc = 1 + dphi/totalarea;
        rsc = pow(rsc, 1/3.);
        for (int i = 0; i < N; i++)
        {
            Dn[i] *= rsc;
            m[i] = M_PI * (Dn[i]*Dn[i]*Dn[i])/6.0;
        }
    }

    // x runs x0,y0,x1,y1...xN,yN, so as the other arrays with 2N
    long nt, Nt;
    double dt2 = dt * dt;
    
    nt = 0;
    Nt = 1e7;
    double * vx = new double[3*N];
    double * ax = new double[3*N];
    double * ax_old = new double[3*N];
    double * Fx = new double[3*N];

    for (int i=0; i<3*N; i++)
    {
        vx[i] = 0.0;
        ax[i] = 0.0;
        ax_old[i] = 0.0;
        Fx[i] = 0.0;
    }

    double P = 0;
    double Ftol;

    // NL playground
    bool NL_flag = false;
    double skin_scale = 1.1;
    bool * NL = new bool[N*N];
    double * x_skin = new double[3*N];
    double * Dn_skin = new double[3*N];
    if (NL_flag) Ftol = buildNeighborList(N,Dn,x,Dn_skin,x_skin,Fx,alpha,gam,skin_scale,NL);
    else Ftol = MD_getFtol(N,Dn,x,Fx,alpha,gam);
    
    // fire coefficients
    double PP = 0;
    long nfiremin = 5;
    double finc = 1.1;
    double fdec = 0.5;
    double astart = 0.1;
    double a = astart;
    double fa = 0.99;
    double dtmax = 10.0 * dt;
    long cut = nt;

    // force energy minimization
    Ftol = max(Ftol,tol*10.0);

    while (Ftol>tol && nt<Nt)
    {
        // periodic BC, warning!!! does not account for gam != 0.
        if (gam == 0)
        {
            for (int i = 0; i < 3*N; i++)
            {
                x[i] = fmod(x[i], 1);
                if (x[i] < 0)
                {
                    x[i] += 1.0;
                }
            }
        }

        // zeroing center of mass velocity
        if (nt%100 == 0) MD_CMzeroing(N,vx,m);
        
        // first step velocity verlet integration
        for (int i=0; i<3*N; i++) x[i] += vx[i]*dt + ax_old[i]*dt2/2.0;

        // get forces from position
        if (NL_flag)
        {
            if (nt > 0 && ifUpdateNL(N,Dn,x,Dn_skin,x_skin,gam)) Ftol = buildNeighborList(N,Dn,x,Dn_skin,x_skin,Fx,alpha,gam,skin_scale,NL);
            else Ftol = getFtolFromNL(N,Dn,x,Fx,alpha,gam,NL);
        }
        else
        {
            Ftol = MD_getFtol(N,Dn,x,Fx,alpha,gam);
        }
        /*
        if (false && Ftol < tol)
        {
            if (false && nt < 5000) Ftol = tol; // wait until we have a break time
            else
            {
                if (false && nt >= 5000 && nt < 0.25 * Nt) // not using
                {
                    int Nr = 0;
                    bool * is_rattler = new bool[N];
                    for (int i=0; i<N; i++) is_rattler[i] = false;
                    MD_getRattler(N, Dn, x, is_rattler, gam);
                    Ftol = 0;
                    for (int i=0; i<N; i++) 
                    {
                        if (~is_rattler[i]) 
                        {
                            Nr += 1;
                            Ftol += (Fx[2*i] * Fx[2*i] + Fx[2*i+1] * Fx[2*i+1] + Fx[2*i+2] * Fx[2*i+2]);
                        }
                    }
                    Ftol /= (Nr*Nr);
                }
            }
        }
        */

        // calculates acceleration
        for (int i=0; i<N; i++) 
        {
            ax[3*i] = Fx[3*i] / m[i];
            ax[3*i+1] = Fx[3*i+1] / m[i];
            ax[3*i+2] = Fx[3*i+2] / m[i];
        }

        // fire
        if (nt > 0)
        {
            PP = FIRE(N,vx,Fx,a);
            if (PP<0)
            {
                for (int i=0; i<3*N; i++) vx[i] = 0.0;
                cut = nt;
                dt = dt * fdec;
                a = astart;
            }
            else
            {
                if (PP>=0 && nt - cut > nfiremin)
                {
                    dt = min(dt*finc, dtmax);
                    a = a * fa;
                }
            }
            dt2 = dt * dt;
        }

        // second step velocity verlet integration
        for (int i=0; i<3*N; i++) vx[i] += (ax_old[i]+ax[i])*dt/2.0;
        for (int i=0; i<3*N; i++) ax_old[i] = ax[i];
        
        nt++;
        // if (nt > 0 && nt%10000 == 0) cout << nt << ",dFtol=" << Ftol-tol << endl;
    }
    // cout << nt << ",Ftol=" << Ftol << endl;
    if (nt == Nt) cout << "Nt too small?" << endl;
    P = MD_getP(N,Dn,x,alpha,gam);
    // cout << P << endl;
    return P;
}

double MD_cmpjam_minimization_VV_FIRE(int N, double * Dn, double * m, double * x, double dt, double alpha, double dphi, double gam, const double& tol)
{
    // first grow particles by dphi
    double totalarea = 0.0;
    for(int i = 0; i < N; i++) totalarea += m[i];
    if (fabs(dphi) > 1e-16)
    {
        double rsc = sqrt(1 + dphi/totalarea);
        for (int i = 0; i < N; i++)
        {
            Dn[i] *= rsc;
            m[i] = M_PI * (Dn[i]*Dn[i]*Dn[i])/6.0;
        }
    }

    // x runs x0,y0,x1,y1...xN,yN, so as the other arrays with 2N
    long nt, Nt;
    double dt2 = dt * dt;
    
    nt = 0;
    Nt = 1e7;
    double * vx = new double[3*N];
    double * ax = new double[3*N];
    double * ax_old = new double[3*N];
    double * Fx = new double[3*N];

    for (int i=0; i<3*N; i++)
    {
        vx[i] = 0.0;
        ax[i] = 0.0;
        ax_old[i] = 0.0;
        Fx[i] = 0.0;
    }

    double P = 0;
    double Ftol;

    // NL playground
    bool NL_flag = false;
    double skin_scale = 1.1;
    bool * NL = new bool[N*N];
    double * x_skin = new double[3*N];
    double * Dn_skin = new double[3*N];
    if (NL_flag) Ftol = buildNeighborList(N,Dn,x,Dn_skin,x_skin,Fx,alpha,gam,skin_scale,NL);
    else Ftol = MD_getFtol(N,Dn,x,Fx,alpha,gam);
    
    // fire coefficients
    double PP = 0;
    long nfiremin = 20; // delaystep
    double finc = 1.1; // dtgrow
    double fdec = 0.5; // dtshrink
    double astart = 0.25; // alpha0
    double a = astart; // alpha
    double fa = 0.99; // alphashrink
    double dtmax = 10.0 * dt; // tmax
    double dtmin = 0.02 * dt; // tmin
    long nt_PP_pos = nt;
    long nt_PP_neg = nt;
    bool initialdelay = true;
    long nt_PP_neg_max = 2000;

    int N_per_coll = 1/dt;

    while (Ftol>tol && nt<Nt)
    {
        // periodic BC, warning!!! does not account for gam != 0.
        if (gam == 0)
        {
            for (int i = 0; i < 3*N; i++)
            {
                x[i] = fmod(x[i], 1);
                if (x[i] < 0)
                {
                    x[i] += 1.0;
                }
            }
        }

        // zeroing center of mass velocity
        if (nt%N_per_coll == 0) MD_CMzeroing(N,vx,m);        

        // first step VV inegration in FIRE 2.0
        if (nt > 0)
        {
            PP = FIRE_getFv(N,vx,Fx);
            if (PP>0)
            {
                nt_PP_pos += 1;
                nt_PP_neg = 0;
                if (nt_PP_pos > nfiremin)
                {
                    dt = min(dt*finc, dtmax);
                    a = a * fa;
                }
            }
            else
            {
                nt_PP_pos = 0;
                nt_PP_neg += 1;
                
                if (nt_PP_neg > nt_PP_neg_max)
                {
                    cout << "too hard!" << endl;
                    break;
                }
                
                if (~(initialdelay && nt < nfiremin))
                {
                    if (dt*fdec >= dtmin) dt = dt * fdec;
                    a = a * fa;
                }
                for (int i=0; i<3*N; i++) x[i] -= (dt/2.0) * vx[i];
                for (int i=0; i<3*N; i++) vx[i] = 0.0;
            }
            /*
            if (PP<0)
            {
                cut = nt;
                dt = dt * fdec;
                a = astart;
                for (int i=0; i<2*N; i++) x[i] -= (dt/2.0) * vx[i];
                for (int i=0; i<2*N; i++) vx[i] = 0.0;
            }
            else
            {
                if (PP>=0 && nt - cut > nfiremin)
                {
                    dt = min(dt*finc, dtmax);
                    a = a * fa;
                }
            }
            */
            dt2 = dt * dt;
            FIRE_mixing(N,vx,Fx,a);
        }

        // second step VV integration in FIRE 2.0
        for (int i=0; i<3*N; i++) vx[i] += (dt/2.0) * ax[i];

        // third step VV inegration in FIRE 2.0
        for (int i=0; i<3*N; i++) x[i] += vx[i]*dt;
        
        // fourth step VV inegration in FIRE 2.0, get forces from position
        if (NL_flag)
        {
            if (nt > 0 && ifUpdateNL(N,Dn,x,Dn_skin,x_skin,gam)) Ftol = buildNeighborList(N,Dn,x,Dn_skin,x_skin,Fx,alpha,gam,skin_scale,NL);
            else Ftol = getFtolFromNL(N,Dn,x,Fx,alpha,gam,NL);
        }
        else
        {
            Ftol = MD_getFtol(N,Dn,x,Fx,alpha,gam);
        }

        if (Ftol < tol)
        {
            if (nt < 5000) 
            {
                Ftol = tol; // wait until we have a break time
            }
            else
            {
                if (false && nt >= 5000 && nt < 0.25 * Nt)
                {
                    int Nr = 0;
                    bool * is_rattler = new bool[N];
                    for (int i=0; i<N; i++) is_rattler[i] = false;
                    // MD_getRattler(N, Dn, x, is_rattler, gam);
                    Ftol = 0;
                    for (int i=0; i<N; i++) 
                    {
                        if (~is_rattler[i]) 
                        {
                            Nr += 1;
                            Ftol += (Fx[2*i] * Fx[2*i] + Fx[2*i+1] * Fx[2*i+1] + Fx[2*i+2] * Fx[2*i+2]);
                        }
                    }
                    Ftol /= (Nr*Nr);
                }
            }
        }

        // "fifth" step VV inegration in FIRE 2.0, get forces from position, calculates acceleration
        for (int i=0; i<N; i++) 
        {
            ax[2*i] = Fx[2*i] / m[i];
            ax[2*i+1] = Fx[2*i+1] / m[i];
            ax[2*i+2] = Fx[2*i+2] / m[i];
        }

        // second step velocity verlet integration
        for (int i=0; i<3*N; i++) vx[i] += ax[i]*dt/2.0;
        // for (int i=0; i<2*N; i++) ax_old[i] = ax[i];

        nt++;
        if (nt > 1e4 && nt%50000 == 0) cout << nt << ",Ftol=" << Ftol << endl;
    }
    if (nt == Nt) cout << "Nt too small?" << endl;
    P = MD_getP(N,Dn,x,alpha,gam);
    return P;
}

double MD_shrjam_minimization(const int &N, double * Dn, double * m, double * x, double dt, const double &alpha, const double &dG, double &gam, const double& tol)
{
    gam += dG;
    for (int i=0; i<N; i++) x[3*i] += dG * x[3*i+1];
    long nt, Nt;
    double dt2 = dt * dt;
    
    nt = 0;
    Nt = 1e7;
    // if (N > 64) Nt = Nt * N/64;

    double * vx = new double[3*N];
    double * ax = new double[3*N];
    double * ax_old = new double[3*N];
    double * Fx = new double[3*N];

    for (int i=0; i<3*N; i++)
    {
        vx[i] = 0.0;
        ax[i] = 0.0;
        ax_old[i] = 0.0;
        Fx[i] = 0.0;
    }
    double P = 0.0;
    double Ftol = MD_getFtol(N,Dn,x,Fx,alpha,gam);

    if (Ftol<tol)
    {
        cout << "bad shear!" << endl;
        Ftol = 1.0;
    } 
    // cout << Ftol << endl;

    // fire coefficients
    double PP = 0;
    long nfiremin = 5;
    double finc = 1.1;
    double fdec = 0.5;
    double astart = 0.1;
    double a = astart;
    double fa = 0.99;
    double dtmax = 10.0 * dt;
    long cut = nt;

    while (Ftol>tol && nt<Nt)
    {
        // zeroing center of mass velocity
        MD_CMzeroing(N,vx,m);

        // first step velocity verlet integration
        for (int i=0; i<3*N; i++) x[i] += vx[i]*dt + ax_old[i]*dt2/2.0;

        // get forces from position
        Ftol = MD_getFtol(N,Dn,x,Fx,alpha,gam);

        if (Ftol < tol && nt < 500)
        {
            Ftol = 1.01 * tol;
        }
        /*
        if (Ftol < tol)
        {
            if (nt < 5000) 
            {
                Ftol = tol; // wait until we have a break time
            }
            else
            {
                if (nt >= 5000 && nt < 0.5 * Nt)
                {
                    int Nr = 0;
                    bool * is_rattler = new bool[N];
                    for (int i=0; i<N; i++) is_rattler[i] = false;
                    MD_getRattler3D(N, Dn, x, is_rattler, gam);
                    Ftol = 0;
                    for (int i=0; i<N; i++) 
                    {
                        if (~is_rattler[i]) 
                        {
                            Nr += 1;
                            Ftol += (Fx[2*i] * Fx[2*i] + Fx[2*i+1] * Fx[2*i+1]);
                        }
                    }
                    Ftol /= (Nr*Nr);
                }
            }
        }
        */

        // calculates acceleration
        for (int i=0; i<N; i++) 
        {
            ax[3*i] = Fx[3*i] / m[i];
            ax[3*i+1] = Fx[3*i+1] / m[i];
            ax[3*i+2] = Fx[3*i+2] / m[i];
        }

        // fire
        if (nt > 0)
        {
            PP = FIRE(N,vx,Fx,a);
            if (PP<0)
            {
                for (int i=0; i<3*N; i++) vx[i] = 0.0;
                cut = nt;
                dt = dt * fdec;
                a = astart;
            }
            else
            {
                if (PP>=0 && nt - cut > nfiremin)
                {
                    dt = min(dt*finc, dtmax);
                    a = a * fa;
                }
            }
            dt2 = dt * dt;
        }

        // second step velocity verlet integration
        for (int i=0; i<3*N; i++) vx[i] += (ax_old[i]+ax[i])*dt/2.0;
        for (int i=0; i<3*N; i++) ax_old[i] = ax[i];
        nt++;
        // if (nt%100==0) cout << Ftol << endl;
    }
    if (nt == Nt) cout << "Nt too small?" << endl;
    P = MD_getP(N,Dn,x,alpha,gam);

    delete[] vx;
    delete[] ax;
    delete[] ax_old;
    delete[] Fx;

    return P;
}

void MD_cmpjam_main(int N, double * Dn, double * m, double * x, double dt, double alpha, double Ptol, double dphi, double gam)
{
    double tol = 1e-7;
    double ftol = 1e-29;
    const double tol0 = tol;
    const double ftol0 = ftol;
    const double dphi0 = dphi;
    const double dt0 = dt;
    
    double P = MD_getP(N,Dn,x,alpha,gam);
    double P0 = P;
    double dphi_exact = pow(Ptol,(1/alpha)) - pow(P0,(1/alpha));
    dphi_exact *= 2.0;
    // if (dphi_exact < dphi) dphi = dphi_exact;
    // cout << dphi << "," << dphi_exact << endl;
    // printf("Current P=%1.5f.\n",P);

    double * Dn_old = new double[N];
    double * m_old = new double[N];
    double * x_old = new double[3*N];
    for (int i=0; i<3*N; i++) x_old[i] = x[i];
    for (int i=0; i<N; i++)
    {
        Dn_old[i] = Dn[i];
        m_old[i] = m[i];
    }

    // save initial state
    double * Dn_ini = new double[N];
    double * m_ini = new double[N];
    double * x_ini = new double[3*N];
    for (int i=0; i<3*N; i++) x_ini[i] = x[i];
    for (int i=0; i<N; i++)
    {
        Dn_ini[i] = Dn[i];
        m_ini[i] = m[i];
    }

    long count = 0;
    long count_max = 1e7;

    bool regular_case = true;
    // bool regular_case = false;
    // if (Ptol < 1e-10) regular_case = false;

    // family scan playground
    bool family_scan = false;
    if (gam >= 1e-5) family_scan = true;
    int count_rejam = 0;

    int fire_ver = 1;

    // cout << Ptol << endl;

    time_t tstart, tend;
    tstart = time(0);

    while (P < Ptol || P > (1+tol)*Ptol)
    {
        if (P < Ptol)// || C < (N-Nr))
        {
            dphi = fabs(dphi);
            // copy current to old
            for (int i=0; i<3*N; i++) x_old[i] = x[i];
            for (int i=0; i<N; i++)
            {
                Dn_old[i] = Dn[i];
                m_old[i] = m[i];
            }
            if (fire_ver == 1) P = MD_cmpjam_minimization(N,Dn,m,x,dt,alpha,dphi,gam,ftol);
            else if (fire_ver == 2) P = MD_cmpjam_minimization_VV_FIRE(N,Dn,m,x,dt,alpha,dphi,gam,ftol);
            else
            {
                cout << "reserved!" << endl;
            }
        }
        else if (P > (1+tol)*Ptol)
        {
            if (regular_case) // regular case
            {
                // copy old to current
                for (int i=0; i<3*N; i++) x[i] = x_old[i];
                for (int i=0; i<N; i++)
                {
                    Dn[i] = Dn_old[i];
                    m[i] = m_old[i];
                }
                dphi /= 2.0;
            }
            else // difficult case
            {
                mt19937 generator(time(0));
                uniform_real_distribution<double> distribution(0.0, 1.0);
                double randscale;
                do {randscale = distribution(generator);} while (randscale < 0.1);
                dphi = -fabs(dphi) * randscale;
            }            
            if (fire_ver == 1) P = MD_cmpjam_minimization(N,Dn,m,x,dt,alpha,dphi,gam,ftol);
            else if (fire_ver == 2) P = MD_cmpjam_minimization_VV_FIRE(N,Dn,m,x,dt,alpha,dphi,gam,ftol);
            else
            {
                cout << "reserved!" << endl;
            }
        }
        // print out here
        double totalarea = 0.0;
        for(int i = 0; i < N; i++) totalarea += m[i];
        
        if (totalarea > 0.5 && count > 0 && count < count_max) printf("%ld,phi=%1.7f,dphi=%e,P=%e.\n",count,totalarea,dphi,P);
        /*
        else if (totalarea > 0.8 && count > 150 && count < 1e3 && count%100 == 0) printf("%ld,phi=%1.7f,dphi=%e,P=%e.\n",count,totalarea,dphi,P);
        else if (totalarea > 0.8 && count > 1e3 && count < 1e4 && count%1000 == 0) printf("%ld,phi=%1.7f,dphi=%e,P=%e.\n",count,totalarea,dphi,P);
        else if (totalarea > 0.8 && count > 1e4 && count%10000 == 0) printf("%ld,phi=%1.7f,dphi=%e,P=%e.\n",count,totalarea,dphi,P);
        */

        count ++;
      
        if ((P < Ptol || P > (1+tol)*Ptol) && (count > count_max || fabs(dphi) < 3.0e-16))
        {   
            count = 0;
            count_rejam += 1;
            regular_case = false;

            mt19937 generator(time(0));
            uniform_real_distribution<double> distribution(0.0, 1.0);
            double randscale;
            do {randscale = distribution(generator);} while (randscale < 0.1);
            dphi = -dphi0 * randscale;
            if (fire_ver == 1) P = MD_cmpjam_minimization(N,Dn,m,x,dt,alpha,dphi,gam,ftol);
            else if (fire_ver == 2) P = MD_cmpjam_minimization_VV_FIRE(N,Dn,m,x,dt,alpha,dphi,gam,ftol);
            else
            {
                cout << "reserved!" << endl;
            }

            double totalarea = 0.0;
            for(int i = 0; i < N; i++) totalarea += m[i];
            // printf("%ld,phi=%1.7f,dphi=%e,P=%e.\n",count,totalarea,dphi,P);
        }

        if ((P < Ptol || P > (1+tol)*Ptol) && (family_scan && count_rejam > 4))
        {
            // experimental perturbation to find new state
            mt19937 generator(time(0));
            uniform_real_distribution<double> distribution(0.0, 1.0);
            double delta = Dn[0] * 0.05; // 5% Diameter
            for (int i=0; i<N; i++)
            {
                x[3*i] += (distribution(generator)-0.5)*delta;
                x[3*i+1] += (distribution(generator)-0.5)*delta;
                x[3*i+2] += (distribution(generator)-0.5)*delta;
            }
            count_rejam = 4; // just some arbitrary number...
        }
    }
    double totalarea = 0.0;
    for(int i = 0; i < N; i++) totalarea += m[i];
    tend = time(0);
    printf("%ld,t=%1.1fs,phi=%1.7f,dphi=%e,P=%e.\n",count,difftime(tend,tstart),totalarea,dphi,P);
}

void MD_cmpjam_main_DS(int N, double * Dn, double * m, double * x, bool * contact_network, double dt, double alpha, double Ptol, double dphi, double gam)
{
    double tol = 1e-7;
    double ftol = 1e-30;
    const double tol0 = tol;
    const double ftol0 = ftol;
    const double dphi0 = dphi;
    const double dt0 = dt;

    double P = MD_getP_DS(N,Dn,x,contact_network,alpha,gam);
    double P0 = P;
    // printf("Initial P=%e.\n",P);

    double * Dn_old = new double[N];
    double * m_old = new double[N];
    double * x_old = new double[2*N];
    for (int i=0; i<2*N; i++) x_old[i] = x[i];
    for (int i=0; i<N; i++)
    {
        Dn_old[i] = Dn[i];
        m_old[i] = m[i];
    }

    // save initial state
    double * Dn_ini = new double[N];
    double * m_ini = new double[N];
    double * x_ini = new double[2*N];
    for (int i=0; i<2*N; i++) x_ini[i] = x[i];
    for (int i=0; i<N; i++)
    {
        Dn_ini[i] = Dn[i];
        m_ini[i] = m[i];
    }

    long count = 0;
    long count_max = 5e5;

    bool regular_case;

    if (P > (1+tol)*Ptol) regular_case = false; // loaded a overcompressed state
    else regular_case = true;

    if (!regular_case) cout << "Overcompressed initial state!" << endl;

    time_t tstart, tend;
    tstart = time(0);

    while (P < Ptol || P > (1+tol)*Ptol)
    {
        if (P < Ptol)
        {
            dphi = fabs(dphi);
            // copy current to old
            for (int i=0; i<2*N; i++) x_old[i] = x[i];
            for (int i=0; i<N; i++)
            {
                Dn_old[i] = Dn[i];
                m_old[i] = m[i];
            }
            P = MD_cmpjam_minimization(N,Dn,m,x,dt,alpha,dphi,gam,ftol);
        }
        else if (P > (1+tol)*Ptol)
        {
            if (regular_case) // regular case
            {
                // copy old to current
                for (int i=0; i<2*N; i++) x[i] = x_old[i];
                for (int i=0; i<N; i++)
                {
                    Dn[i] = Dn_old[i];
                    m[i] = m_old[i];
                }
                dphi /= 2.0;
            }
            else // difficult case
            {
                mt19937 generator(time(0));
                uniform_real_distribution<double> distribution(0.0, 1.0);
                double randscale;
                do {randscale = distribution(generator);} while (randscale < 0.1);
                dphi = -fabs(dphi) * randscale;
            }            
            
            P = MD_cmpjam_minimization(N,Dn,m,x,dt,alpha,dphi,gam,ftol);
        }
        // update P here based on contact network
        P = MD_getP_DS(N,Dn,x,contact_network,alpha,gam);

        // print out here
        double totalarea = 0.0;
        for(int i = 0; i < N; i++) totalarea += m[i];

        /*
        if (totalarea > 0.8 && count > 50 && count < 150) printf("%ld,phi=%1.7f,dphi=%e,P=%e.\n",count,totalarea,dphi,P);
        else if (totalarea > 0.8 && count > 150 && count < 1e3 && count%100 == 0) printf("%ld,phi=%1.7f,dphi=%e,P=%e.\n",count,totalarea,dphi,P);
        else if (totalarea > 0.8 && count > 1e3 && count < 1e4 && count%1000 == 0) printf("%ld,phi=%1.7f,dphi=%e,P=%e.\n",count,totalarea,dphi,P);
        else if (totalarea > 0.8 && count > 1e4 && count%10000 == 0) printf("%ld,phi=%1.7f,dphi=%e,P=%e.\n",count,totalarea,dphi,P);
        */

        count ++;
      
        if (count > count_max || fabs(dphi) < 3.0e-16)
        {   
            count = 0;
            regular_case = false;

            mt19937 generator(time(0));
            uniform_real_distribution<double> distribution(0.0, 1.0);
            double randscale;
            do {randscale = distribution(generator);} while (randscale < 0.1);
            dphi = -dphi0 * randscale;
            P = MD_cmpjam_minimization(N,Dn,m,x,dt,alpha,dphi,gam,ftol);

            double totalarea = 0.0;
            for(int i = 0; i < N; i++) totalarea += m[i];
            // printf("%ld,phi=%1.7f,dphi=%e,P=%e.\n",count,totalarea,dphi,P);
        }
    }
    double totalarea = 0.0;
    for(int i = 0; i < N; i++) totalarea += m[i];
    tend = time(0);
    printf("%ld,t=%1.1fs,phi=%1.7f,dphi=%e,P=%e.\n",count,difftime(tend,tstart),totalarea,dphi,P);
}

void scale(int N, double * Dn, double * m, double rescale)
{
    double D3;
    for(int i = 0; i < N; i++)
    {
        Dn[i] = Dn[i] * rescale;
        D3 = Dn[i] * Dn[i] * Dn[i];
        m[i] = M_PI * D3 / 6.0;
    }
}

//check if any particles are touching
bool anytouch(int N, double * pos, double * sc)
{
	for(int i=0; i<N-1; i++)
	{
		for(int j=i+1; j<N; j++)
		{
			//vector pointing from particle i to j
			double rij[3];
			//periodic boundary conditions
			rij[0] = pos[3*j]-pos[3*i]-round(pos[3*j]-pos[3*i]); // x
			rij[1] = pos[3*j+1]-pos[3*i+1]-round(pos[3*j+1]-pos[3*i+1]); // y
            rij[2] = pos[3*j+2]-pos[3*i+2]-round(pos[3*j+2]-pos[3*i+2]); // z
			//sum of 2 radii
			double bigR = sc[i]+sc[j];
			//they touch if the distance between them is less than bigR
			if(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2] < bigR*bigR) return true;
		}
	}
	return false;
}

// check if file exists
bool fileExist(const char * name)
{
    if (FILE *file = fopen(name, "r")) 
    {
        fclose(file);
        return true;
    } 
    else 
    {
        return false;
    }
}

// check if path exists
bool IsPathExist(const char * s)
{
  struct stat buffer;
  return (stat (s, &buffer) == 0);
}

// write to output
void writeResult(const char * path, int N, double * Dn, double * m, double * x)
{
    FILE *xyz = fopen(path, "w+");
    for (int j=0; j<N; j++) fprintf(xyz, "%1.20f ", Dn[j]);
    fprintf(xyz, "\n");
    for (int j=0; j<N; j++) fprintf(xyz, "%1.20f ", m[j]);
    fprintf(xyz, "\n");
    for (int j=0; j<N; j++) fprintf(xyz, "%1.20f ", x[3*j]);
    fprintf(xyz, "\n");
    for (int j=0; j<N; j++) fprintf(xyz, "%1.20f ", x[3*j+1]);
    fprintf(xyz, "\n");
    for (int j=0; j<N; j++) fprintf(xyz, "%1.20f ", x[3*j+2]);
    fprintf(xyz, "\n");

    fclose(xyz);
}

// write to output
void writeResultToFile(FILE * xyz, const int &N, double * Dn, double * m, double * x)
{
    for (int j=0; j<N; j++) fprintf(xyz, "%1.20f ", Dn[j]);
    fprintf(xyz, "\n");
    for (int j=0; j<N; j++) fprintf(xyz, "%1.20f ", m[j]);
    fprintf(xyz, "\n");
    for (int j=0; j<N; j++) fprintf(xyz, "%1.20f ", x[3*j]);
    fprintf(xyz, "\n");
    for (int j=0; j<N; j++) fprintf(xyz, "%1.20f ", x[3*j+1]);
    fprintf(xyz, "\n");
    for (int j=0; j<N; j++) fprintf(xyz, "%1.20f ", x[3*j+2]);
    fprintf(xyz, "\n");
}

// load result
void loadResult(const char * path, int N, double * Dn, double * m, double * x)
{
    FILE *xyz = fopen(path, "r");
    for (int j=0; j<N; j++) fscanf(xyz, "%lf", &Dn[j]);
    for (int j=0; j<N; j++) fscanf(xyz, "%lf", &m[j]);
    for (int j=0; j<N; j++) fscanf(xyz, "%lf", &x[3*j]);
    for (int j=0; j<N; j++) fscanf(xyz, "%lf", &x[3*j+1]);
    for (int j=0; j<N; j++) fscanf(xyz, "%lf", &x[3*j+2]);
    fclose(xyz);
}

void loadResultFromFile(FILE * xyz, const int &N, double * Dn, double * m, double * x)
{
    for (int j=0; j<N; j++) fscanf(xyz, "%lf", &Dn[j]);
    for (int j=0; j<N; j++) fscanf(xyz, "%lf", &m[j]);
    for (int j=0; j<N; j++) fscanf(xyz, "%lf", &x[3*j]);
    for (int j=0; j<N; j++) fscanf(xyz, "%lf", &x[3*j+1]);
    for (int j=0; j<N; j++) fscanf(xyz, "%lf", &x[3*j+2]);
}

void MD_shearModulus_main(const int &N, const double &alpha, const char * loadpath, const char * savepath, const char * saveCPpath, const double &Ptol, const bool &getCPstate, const bool &positive_shear, const bool &HQ_flag)
{
    FILE *abc = fopen(loadpath, "r");
    if (abc == NULL)
    {
        cout << loadpath << " not exists!" << endl;
        return;
    }
    // array
    double * Dn = new double[N];
    double * m = new double[N];
    double * x = new double[3*N];
    loadResultFromFile(abc, N, Dn, m, x);
    fclose(abc);

    // measurements
    double * Fx = new double[3*N];
    double stress[8] = {0};

    // simulation parameters
    double dG;
    int N_shear_steps = 20;
    int curr_steps = 0;
    double dt = 5e-4;
    double gam = 0.0;
    double ftol = 1e-29;

    if (positive_shear) dG = 1e-9;
    else dG = -1e-9;

    if (HQ_flag && Ptol < 1e-5) ftol = 1e-31;
    {
        double P = MD_cmpjam_minimization(N,Dn,m,x,dt,alpha,0.0,gam,ftol);
    }
    
    // open file to save sheared state, write gam=0 state
    FILE *xyz = fopen(savepath, "w+");
    writeResultToFile(xyz, N, Dn, m, x);
    FILE *cpstate = NULL;
    if (getCPstate)
    {
        cpstate = fopen(saveCPpath, "w+");
        writeResultToFile(cpstate, N, Dn, m, x);
    }

    while (curr_steps < N_shear_steps)
    {
        // shear at constant volume to gam+dG
        MD_shrjam_minimization(N, Dn, m, x, dt, alpha, dG, gam, ftol);
        for (int i=0; i<2*N; i++) Fx[i] = 0.0;
        double Ftol = MD_getFtol(N,Dn,x,Fx,alpha,gam);
        MD_getStressTensor(N, Dn, x, stress, alpha, gam);
        cout << "Ftol=" << Ftol << ",sigma_xy=" << -stress[1] << endl;

        // write to output file, do not close
        writeResultToFile(xyz, N, Dn, m, x);

        // find fixed pressure state if necessary
        // before compression need to save the state
        if (getCPstate)
        {
            // save initial state
            double * Dn_ini = new double[N];
            double * m_ini = new double[N];
            double * x_ini = new double[3*N];
            for (int i=0; i<3*N; i++) x_ini[i] = x[i];
            for (int i=0; i<N; i++)
            {
                Dn_ini[i] = Dn[i];
                m_ini[i] = m[i];
            }
            double dphi = 1e-5;
            MD_cmpjam_main(N, Dn, m, x, dt, alpha, Ptol, dphi, gam);

            // write fixed pressure state
            writeResultToFile(cpstate, N, Dn, m, x);
            
            // recover initial sheared state
            for (int i=0; i<3*N; i++) x[i] = x_ini[i];
            for (int i=0; i<N; i++)
            {
                Dn[i] = Dn_ini[i];
                m[i] = m_ini[i];
            }

            delete[] Dn_ini;
            delete[] m_ini;
            delete[] x_ini;
        }      
        cout << "Step=" << curr_steps << ",gamma=" << gam << ",P=" << Ptol << " completed." << endl;
        curr_steps += 1;
    }

    fclose(xyz);
    if (cpstate != NULL) fclose(cpstate);

    delete[] Dn;
    delete[] m;
    delete[] x;
    delete[] Fx;
}

// STILL PENDING
void MD_mapping_shear_func(const int &N, const double &alpha, const char * loadpath, double * Plist, const double &dG)
{
    FILE *abc = fopen(loadpath, "r");
    if (abc == NULL)
    {
        cout << loadpath << " not exists!" << endl;
        return;
    }
    // array
    double * Dn = new double[N];
    double * m = new double[N];
    double * x = new double[2*N];
    loadResultFromFile(abc, N, Dn, m, x);
    fclose(abc);

    // measurements
    double * Fx = new double[2*N];
    double stress[4] = {0};

    // step 1. shear at const. volume, to get P_list
    // always load p=0 state -> correspond to lowest pressure

    // simulation parameters
    // double dG = 1e-9;
    int N_shear_steps = 20;
    int curr_steps = 0;
    double dt = 5e-4;
    double gam = 0.0;
    double ftol = 1e-30;
    int Nstates = 21;

    // get contact network
    bool * contact_network = new bool[N*N];
    MD_getCtcNwk(N, Dn, x, contact_network, gam);

    Plist[0] = MD_getP(N, Dn, x, alpha, gam);

    while (curr_steps < N_shear_steps)
    {
        // shear at constant volume to gam+dG
        MD_shrjam_minimization(N, Dn, m, x, dt, alpha, dG, gam, ftol);
        for (int i=0; i<3*N; i++) Fx[i] = 0.0;
        double Ftol = MD_getFtol(N,Dn,x,Fx,alpha,gam);
        MD_getStressTensor(N, Dn, x, stress, alpha, gam);
        cout << "Ftol=" << Ftol << ",sigma_xy=" << -stress[1] << endl;
        
        curr_steps += 1;
        Plist[curr_steps] = MD_getP_DS(N, Dn, x, contact_network, alpha, gam);
        // cout << "Step=" << curr_steps << ",gamma=" << gam << ",Pds=" << Plist[curr_steps] << ",P=" << MD_getP(N, Dn, x, alpha, gam) << " completed." << endl;
        // cout << "Step=" << curr_steps << ",gamma=" << gam << ",P=" << Plist[curr_steps] << " completed." << endl;
    }

    delete[] Dn;
    delete[] m;
    delete[] x;
    delete[] Fx;

    cout << "Done shearing!" << endl;
}

// STILL PENDING
void MD_mapping_shearCP_func(const int &N, const double &alpha, const char * loadpath, const char * savepath, const double &Ptol, const double &dG)
{
    FILE *abc = fopen(loadpath, "r");
    if (abc == NULL)
    {
        cout << loadpath << " not exists!" << endl;
        return;
    }
    // array
    double * Dn = new double[N];
    double * m = new double[N];
    double * x = new double[2*N];
    loadResultFromFile(abc, N, Dn, m, x);
    fclose(abc);

    // Simulation parameters
    double dt = 5e-4;
    double gam = 0.0;
    double ftol = 1e-30;

    // get contact network
    bool * contact_network = new bool[N*N];
    MD_getCtcNwk(N, Dn, x, contact_network, gam);

    // step 1: compress to target pressure and write to file, at zero gamma
    double dphi = 1e-5;
    MD_cmpjam_main_DS(N, Dn, m, x, contact_network, dt, alpha, Ptol, dphi, gam);
    FILE *xyz = fopen(savepath, "w+");
    writeResultToFile(xyz, N, Dn, m, x);

    // step 2: shear and find const. pressure state
    // double dG = 1e-9;
    int N_shear_steps = 20;
    int curr_steps = 0;
    while (curr_steps < N_shear_steps)
    {
        // shear at constant volume to gam+dG
        MD_shrjam_minimization(N, Dn, m, x, dt, alpha, dG, gam, ftol);

        // get const. P state
        // MD_cmpjam_main(N, Dn, m, x, dt, alpha, Ptol, dphi, gam);
        MD_cmpjam_main_DS(N, Dn, m, x, contact_network, dt, alpha, Ptol, dphi, gam);

        // write to file
        writeResultToFile(xyz, N, Dn, m, x);

        cout << "Step=" << curr_steps << ",gamma=" << gam << ",P=" << MD_getP(N, Dn, x, alpha, gam) << " completed." << endl;
        curr_steps += 1;
    }

    fclose(xyz);

    delete[] contact_network;
}

int getNumberofLines(const char * path)
{
    int count = 0;
    string line;
 
    /* Creating input filestream */ 
    ifstream file(path);
    while (getline(file, line)) count++;
    
    cout << "Numbers of lines in the file : " << count << endl;
    file.close();
    return count;
}

void loadResultwLineNumber(const char * path, int N, double * Dn, double * m, double * x, const int &target_lines)
{
    FILE *xyz = fopen(path, "r");
    int curr_line = 0;
    int total_lines = getNumberofLines(path);
    while (curr_line < total_lines && curr_line != target_lines)
    {
        for (int j=0; j<N; j++) fscanf(xyz, "%lf", &Dn[j]);
        for (int j=0; j<N; j++) fscanf(xyz, "%lf", &m[j]);
        for (int j=0; j<N; j++) fscanf(xyz, "%lf", &x[3*j]);
        for (int j=0; j<N; j++) fscanf(xyz, "%lf", &x[3*j+1]);
        for (int j=0; j<N; j++) fscanf(xyz, "%lf", &x[3*j+2]);
        curr_line += 5;
    }
    fflush(xyz);
    fclose(xyz);
}

// Here is the const. pressure family scan
void MD_mapping_CPfamily_func(const int &N, const double &alpha, const char * loadpath, const char * savepath, const double &Ptol, const double &dG, const bool &early_termination)
{
    // step 1. load current state at a given pressure
    FILE *abc = fopen(loadpath, "r");
    if (abc == NULL)
    {
        cout << loadpath << " not exists!" << endl;
        return;
    }

    // step 2. restart if necessary
    // open the text file if exist, read number of lines to determine current gamma, then start from there
    // else start from gam = 0
    double gam = 0.0;
    bool restart_shear = false;

    // Simulation parameters
    double dt = 5e-4;
    double ftol = 1e-30;

    // array
    double * Dn = new double[N];
    double * m = new double[N];
    double * x = new double[2*N];

    FILE *xyz;

    if (!IsPathExist(savepath)) // new start
    {
        loadResultFromFile(abc, N, Dn, m, x);
        fflush(abc);
        fclose(abc);
        gam = 0.0;
        xyz = fopen(savepath, "w+");
    }
    else // restart
    {
        int countLines = getNumberofLines(savepath);
        div_t divresult = div(countLines,4);
        int completed_states;
        if (divresult.rem == 0) completed_states = divresult.quot-1;
        else completed_states = divresult.quot;
        int target_lines = completed_states * 4;
        loadResultwLineNumber(savepath,N,Dn,m,x,target_lines);
        // for (int i=0;i<N;i++) printf("%1.16f, ",x[2*i]);
        // cout << endl;
        gam = (completed_states-1) * dG;
        cout << "Restart with gamma=" << gam << endl;
        restart_shear = true;
        xyz = fopen(savepath, "r+");
        // read until reached target_lines
        int curr_line = 0;
        double * temp = new double[N];
        while (curr_line < target_lines)
        {
            for (int j=0; j<N; j++) fscanf(xyz, "%lf", &temp[j]);
            curr_line += 1;
        }
        fflush(xyz);
        fprintf(xyz, "\n");
    }
    
    // get contact network at zero strain
    bool * contact_network_0 = new bool[N*N];

    // step 2: compress to target pressure and write to file, at zero gamma
    double dphi = 1e-5;
    if (!restart_shear) // not restarting, start a new state
    {
        MD_cmpjam_main(N, Dn, m, x, dt, alpha, Ptol, dphi, gam);
        MD_getCtcNwk(N, Dn, x, contact_network_0, gam);
        writeResultToFile(xyz, N, Dn, m, x);
    }

    // early termination playground - ver. 0 -> lazy termination using max \gamma
    bool contact_network_saved = false;
    bool * contact_network = new bool[N*N];
    double gam_max = 1.0;
    if (early_termination) gam_max = 0.1;

    // step 2: shear and find const. pressure state
    // double dG = 1e-9;
    int curr_steps = 0;
    while (gam < gam_max)
    {
        // shear at constant volume to gam+dG
        double P_shear = MD_shrjam_minimization(N, Dn, m, x, dt, alpha, dG, gam, ftol);
        // cout << P_shear << endl;

        // get const. P state
        MD_cmpjam_main(N, Dn, m, x, dt, alpha, Ptol, dphi, gam);
        // MD_cmpjam_main_DS(N, Dn, m, x, contact_network, dt, alpha, Ptol, dphi, gam);

        
        if (!contact_network_saved && gam !=0 && gam < 2*dG) // sheared once
        {
            contact_network_saved = true;
            MD_getCtcNwk(N, Dn, x, contact_network_0, gam);
        }
        

        // write to file
        writeResultToFile(xyz, N, Dn, m, x);

        cout << "Step=" << curr_steps << ",gamma=" << gam << ",P=" << MD_getP(N, Dn, x, alpha, gam) << " completed." << endl;
        curr_steps += 1;
    }

    fclose(xyz);

    delete[] contact_network;
}