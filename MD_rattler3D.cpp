#include<math.h>
#include<vector>

using namespace std;

void MD_rattler3D_cleanByCn(const int N, const double * Dn, const double * x, bool * is_rattler, int * Cn, const double gam)
{
    double dx, dy, dz, im, Dnm, dnm;
    int C = 0;
    int Nr = 0;
    bool if_over = false;
    while (!if_over)
    {
        if_over = true;
        C = 0;
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
                            C += 2;
                        }
                    }
                }
                
            }
        }
        Nr = 0;
        for (int i=0; i<N; i++)
        {
            if (!is_rattler[i] && Cn[i] < 4) 
            {
                if_over = false;
                is_rattler[i] = true;
            }
            if (!is_rattler[i]) Nr += 1;
        }
    }
}

vector<double> MD_rattler3D_getNeighbors(const int N, const double * Dn, const double * x, bool * is_rattler, const double gam, const int nn)
{
    double dx, dy, dz, im, dnm, Dnm, dnm_ideal, dnm2_ideal;
    double x_temp, y_temp, z_temp;
    vector<double> neighbors;
    // neighbors are the contact points exactly at the surface, free of overlap
    // there are 27 possible positions, with exactly only one that has a seperation less than average diameter
    for (int mm=0; mm<N; mm++)
    {
        if (!is_rattler[mm] && mm != nn)
        {
            Dnm = 0.5 * (Dn[nn]+Dn[mm]);
            dy = x[3*mm+1]-x[3*nn+1];
            im = round(dy);
            dy -= im;
            dx = x[3*mm]-x[3*nn];
            dx -= (round(dx-im*gam)+im*gam);
            dz = x[3*mm+2]-x[3*nn+2];
            dz -= round(dz);
            dnm = sqrt(dx*dx+dy*dy+dz*dz);
            if (dnm < Dnm)
            {
                for (int i=-1; i<2; i++)
                {
                    for (int j=-1; j<2; j++)
                    {
                        for (int k=-1; k<2; k++)
                        {
                            z_temp = x[3*mm+2]+k;
                            y_temp = x[3*mm+1]+j;
                            x_temp = x[3*mm]+i+gam*j;
                            dy = y_temp-x[3*nn+1];
                            dx = x_temp-x[3*nn];
                            dz = z_temp-x[3*nn+2];
                            dnm2_ideal = dx*dx+dy*dy+dz*dz;
                            if (dnm2_ideal < Dnm*Dnm)
                            {
                                dnm_ideal = sqrt(dnm2_ideal);
                                neighbors.push_back(dx/dnm_ideal);
                                neighbors.push_back(dy/dnm_ideal);
                                neighbors.push_back(dz/dnm_ideal);
                            }
                        }
                    }
                }
            }
        }
    }
    return neighbors;
}

double euclid3_norm(double * n)
{
    double n_norm = 0;
    for (int dim=0; dim<3; dim++) n_norm += (n[dim] * n[dim]);
    n_norm = sqrt(n_norm);
    return n_norm;
}

double get_acos3(double * v1, double * v2)
{
    double v1_norm = euclid3_norm(v1);
    double v2_norm = euclid3_norm(v2);
    double inner = 0;
    for (int dim=0; dim<3; dim++) inner += (v1[dim] * v2[dim]);
    inner /= (v1_norm * v2_norm);
    return acos(inner);
}

bool MD_rattler3D_check(vector<double> neighbors)
{
    // ***************************************** //
    // **** this only works for 4 contacts! **** //
    // ***************************************** //

    int Np = neighbors.size() / 3; // number of neighboring particles
    // preparation: enumerates indices
    // vector<int> all_idx;
    // for (int i=0; i<Np; i++) all_idx.push_back(i);
    // preparation: use the combination table for C(4,3)
    int pairs[12] = {0,1,2,0,1,3,0,2,3,1,2,3};
    int candidates[4] = {3,2,1,0};

    for (int pp=0; pp<Np; pp++)
    {
        int idx[3];
        double v0[3], v1[3], v2[3], n[3], vp[3];
        for (int nei=0; nei<3; nei++) idx[nei] = pairs[pp*3+nei]; // 0,1,2
        for (int dim=0; dim<3; dim++)
        {
            v0[dim] = -neighbors[idx[0]*3+dim];
            v1[dim] = neighbors[idx[1]*3+dim] - neighbors[idx[0]*3+dim];
            v2[dim] = neighbors[idx[2]*3+dim] - neighbors[idx[0]*3+dim];
            vp[dim] = -neighbors[candidates[pp]*3+dim];
        }
        // cross product to get n
        n[0] = v1[1] * v2[2] - v1[2] * v2[1];
        n[1] = -(v1[0] * v2[2] - v1[2] * v2[0]);
        n[2] = v1[0] * v2[1] - v1[1] * v2[0];
        double n_norm = euclid3_norm(n);
        for (int dim=0; dim<3; dim++) n[dim] /= n_norm;
        double angle;
        angle = get_acos3(n,v0);
        if (angle > M_PI/2.0)
        {
            for (int dim=0; dim<3; dim++) n[dim] = -n[dim];
        }
        angle = get_acos3(n,vp);
        if (angle < M_PI/2.0) return true;
    }

    return false;
}

void MD_getRattler3D(const int N, const double * Dn, const double * x, bool * is_rattler, const double gam)
{
    // step 1. reset is_rattlers
    for (int i=0; i<N; i++) is_rattler[i] = false;

    // prepare Cn
    int * Cn = new int[N];
    for (int i=0; i<N; i++) Cn[i] = 0;

    bool done = false;
    while (!done)
    {
        MD_rattler3D_cleanByCn(N,Dn,x,is_rattler,Cn,gam);
        done = true;
        for (int nn=0; nn<N; nn++)
        {
            if (!is_rattler[nn] && Cn[nn] == 4)
            {
                // step 2. prepare neighbors
                vector<double> neighbors = MD_rattler3D_getNeighbors(N,Dn,x,is_rattler,gam,nn);

                // step 3. determine if this is a rattler
                if (MD_rattler3D_check(neighbors)) // this is rattler
                {
                    is_rattler[nn] = true;
                    done = false;
                    break;
                }
            }
        }
    }
    delete[] Cn;
}