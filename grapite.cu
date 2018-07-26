#include <cstdio>
#include <algorithm>
#include <thrust/pair.h>
#include <thrust/host_vector.h>
#include <thrust/partition.h>
#include "../etics/src/common.hpp"
#include "../etics/src/scf.hpp"

#define SCFSTEP 0.125
#define FORCEHISTORY 4

extern "C" {
    void __grape_g6_open(int clusterid);
    void __grape_g6_set_ti(int clusterid, double ti);
    int  __grape_g6_set_j_particle(int clusterid, int address, int index, double tj, double dtj, double mass,
                                   double a2by18[3], double a1by6[3], double aby2[3], double v[3], double x[3]);
    void __grape_g6calc_firsthalf(int clusterid, int nj, int ni, int index[], double xi[][3], double vi[][3],
                                  double fold[][3], double jold[][3], double phiold[], double eps2, double h2[]);
    int  __grape_g6calc_lasthalf(int clusterid, int nj, int ni, int index[], double xi[][3], double vi[][3],
                                 double eps2, double h2[], double acc[][3], double jerk[][3], double pot[]);
}

namespace grapite {
    int t_scf = 0;
}

struct PartitoningPredicate {
    __host__ __device__ bool operator()(const thrust::pair<int,int> &Element) {return (Element.second == 0);}
};

class TranslatorClass {
  public:
    TranslatorClass() {N=0;}
    TranslatorClass(int _N) {
        N = _N;
        Dictionary.resize(N);
        ReverseDictionary.resize(N);
        StatusByOriginalIndex.resize(N);
    }
    ~TranslatorClass() {} // You should free the dictionaries.
    void AddParticle(int Index, int Status) {
        Dictionary[Index] = thrust::pair<int,int>(Index, Status);
        StatusByOriginalIndex[Index] = Status;
    }
    int Partition(int Size) {
        int NewSize =  thrust::partition(Dictionary.begin(), Dictionary.begin() + Size, PartitoningPredicate()) - Dictionary.begin();
        for (int Index = 0; Index < Size; Index++) {
            ReverseDictionary[Dictionary[Index].first] = Index;
        }
        return NewSize;
        // After using this, please don't add new particles.
    }
    int PartedToOrig(int Index) {
        return Dictionary[Index].first;
    }
    int OrigToParted(int Index) {
        return ReverseDictionary[Index];
    }
  private:
    int N; // Needed?
    thrust::host_vector< thrust::pair<int,int> > Dictionary;
    thrust::host_vector<int> ReverseDictionary;
    thrust::host_vector<int> StatusByOriginalIndex;
};

template <int width> class MemoryBufferClass {
  public:
    MemoryBufferClass() {}
    MemoryBufferClass(int _max_size, int _N) {
        assert(_max_size >= 2*_N); // Actual minimum is N+1, but it's only effective when the max_size is larger
        max_size = _max_size;
        N = _N;
        size = 0;
        prop = new double[width*max_size];
        address = new int[max_size];
        // GPU arrays
        cudaMalloc((void**)&buffer_address_gpu, N*sizeof(int)); // This one though will never be larger than N
        cudaMalloc((void**)&buffer_prop_gpu, width*N*sizeof(double));
        // The following three arrays are needed for the shrink function
        prop_new = new double[width*max_size];
        address_new = new int[max_size];
        already_updated = new bool[N];
    }
    ~MemoryBufferClass() {}
    void add_particle(int _address, double *_prop) {
        if (size == max_size) this->shrink();
        address[size] = _address;
        std::copy(_prop, _prop+width, prop+size*width);
        address[size] = _address;
        size++;
    }
    void reset() {
        size = 0;
    }
    void shrink() {
        memset(already_updated, 0, N*sizeof(bool));
        int index_new = 0;
        for (int i = size-1; i >= 0; i--) {
            int current_address = address[i];
            if (already_updated[current_address]) continue;
            address_new[index_new] = current_address;
            std::copy(prop+i*width, prop+(i+1)*width, prop_new+index_new*width);
            index_new++;
            already_updated[current_address] = true;
        }
        size = index_new;
        // swap the new arrays with the old ones, nothing to free.
        int *tmp_int_pointer = address;
        address = address_new;
        address_new = tmp_int_pointer;
        double *tmp_double_pointer = prop;
        prop = prop_new;
        prop_new = tmp_double_pointer;
    }
    void h2d() {
        this->shrink();
        cudaMemcpy(buffer_address_gpu, address, sizeof(int)*size,    cudaMemcpyHostToDevice);
        cudaMemcpy(buffer_prop_gpu,    prop,    sizeof(double)*width*size, cudaMemcpyHostToDevice);
    }
    int max_size, N;
    int size;
    int *address;
    double *prop;
    bool *already_updated;
    int *address_new;
    double *prop_new;
    // The following are on the GPU
    int *buffer_address_gpu;
    double *buffer_prop_gpu;
};

#define INTERPOLYEPOCHS FORCEHISTORY
#define INTERPOLYFUNCNUM 3
// it would be nicer if this class uses templates instead of preprocessor macros.
class InterPoly {
  public:
    InterPoly() {
        memset(data, 0, INTERPOLYFUNCNUM*INTERPOLYEPOCHS*sizeof(double));
        memset(x, 0, INTERPOLYEPOCHS*sizeof(double));
        pos = 0;
    }
    ~InterPoly() {}
    void push(double x_, double y[INTERPOLYEPOCHS]) {
        pos = (pos-1+INTERPOLYEPOCHS) % INTERPOLYEPOCHS; // we push from top down
        x[pos] = x_;
        for (int i=0; i<INTERPOLYFUNCNUM; i++) data[i][pos] = y[i];
    }
    void add_to_current(double y[INTERPOLYEPOCHS]) {
        for (int i=0; i<INTERPOLYFUNCNUM; i++) data[i][pos] += y[i];
    }
    double xnow() {
        return x[pos];
    }
    void derivnow(double *result) {
        double y[INTERPOLYEPOCHS];
        for (int func = 0; func < INTERPOLYFUNCNUM; func++) {
            memcpy(y, data[func], INTERPOLYEPOCHS*sizeof(double));
            result[func] = 0;
            for (int j=1; j<INTERPOLYEPOCHS; j++) {
                for (int i=0; i<INTERPOLYEPOCHS-j; i++) y[(pos+i)%INTERPOLYEPOCHS] = (y[(pos+i+1)%INTERPOLYEPOCHS]-y[(pos+i)%INTERPOLYEPOCHS])/(x[(pos+i+j)%INTERPOLYEPOCHS]-x[(pos+i)%INTERPOLYEPOCHS]);
                double term = 1.0;
                for (int i=1; i<j; i++) term *= (x[pos] - x[(pos+i)%INTERPOLYEPOCHS]);
                result[func] += (y[pos]*term);
            }
        }
    }
#if 0
// This calculates the actual value of the polynomial at a point; was only needed for testing the algorithm.
    void value(double arg, double *result) {
        double y[INTERPOLYEPOCHS];
        for (int func = 0; func < INTERPOLYFUNCNUM; func++) {
            memcpy(y, data[func], INTERPOLYEPOCHS*sizeof(double));
            result[func] = y[pos];
            for (int j=1; j<INTERPOLYEPOCHS; j++) {
                for (int i=0; i<INTERPOLYEPOCHS-j; i++) y[(pos+i)%INTERPOLYEPOCHS] = (y[(pos+i+1)%INTERPOLYEPOCHS]-y[(pos+i)%INTERPOLYEPOCHS])/(x[(pos+i+j)%INTERPOLYEPOCHS]-x[(pos+i)%INTERPOLYEPOCHS]);
                double term = 1.0;
                for (int i=0; i<j; i++) term *= (arg - x[(pos+i)%INTERPOLYEPOCHS]);
                result[func] += (y[pos]*term);
            }
        }
    }
#endif
    double data[INTERPOLYFUNCNUM][INTERPOLYEPOCHS];
    double x[INTERPOLYEPOCHS];
    int pos;
};
#undef INTERPOLYEPOCHS
#undef INTERPOLYFUNCNUM

class JerkCalculatorClass {
  public:
    JerkCalculatorClass() {}
    JerkCalculatorClass(int N) {
        ForcePoly = new InterPoly[N];
        PastSteps = new int[N](); // initialized to zero
    }
    ~JerkCalculatorClass() {}
    void load(int j, double t, vec3 F) {
        if (t <= ForcePoly[j].xnow()) return; // maybe best to raise an error. We only go forward; it's possible to make a special provision to rewrite the last force if t == ForcePoly[j].xnow(), but we'll ignore this case.
        double Array[3];
        Array[0] = F.x;
        Array[1] = F.y;
        Array[2] = F.z;
        ForcePoly[j].push(t, Array);
        PastSteps[j]++;
    }
    void update_add(int j, vec3 F) {
        double Array[3];
        Array[0] = F.x;
        Array[1] = F.y;
        Array[2] = F.z;
        ForcePoly[j].add_to_current(Array);
    }
    void add_jerk(int j, double jerk[3]) {
        if (PastSteps[j] < (FORCEHISTORY+1)) return;
        double Array[3];
        ForcePoly[j].derivnow(Array);
        jerk[0] += Array[0];
        jerk[1] += Array[1];
        jerk[2] += Array[2];
    }

    InterPoly *ForcePoly;
    int *PastSteps;
};

namespace grapite {
    // Global variables in this namespace
    int rank; // only needed for debugging purposes
    double ti, t_exp;
    TranslatorClass Translator_j;
    TranslatorClass Translator_i;
    int *status;
    MemoryBufferClass<3> PositionBuffer;
    MemoryBufferClass<1> MassBuffer;
    Particle *Particles_j;
    Particle *Particles_i;
    Particle *Particles_ih;
    etics::scf::scfclass SCFObject_cor;
    etics::scf::scfclass SCFObject_hal;
    int nj_total, nj_core, ni_total, ni_core;
    double (*xi_parted)[3], (*vi_parted)[3], (*acc_parted)[3], (*jerk_parted)[3], *pot_parted;
    int *index_parted;
    double *Potential, *Potential_h;
    vec3 *F, *F_h;
    JerkCalculatorClass JerkCalculator;

    void initialize_globals(int N, int n_loc);
    void generate_angular_momentum_mask(int N, double x[][3], double v[][3], int skip, double fraction, int output[]);
    void set_j_particle_pos(int j, double *x);
    void set_j_particle_mass(int j, double m);
    int partition_arrays_i(int n_act, int ind_act[], double x_act_array[][3], double v_act_array[][3], int ind_act_parted[], double x_act_parted[][3], double v_act_parted[][3]);
    void revert_arrays_i(double acc_parted[][3], double jerk_parted[][3], double pot_parted[], double acc[][3], double jerk[][3], double pot[]);
    __global__ void ScatterPositionsInMemory(int *address, double *pos, int size, Particle *ParitcleList);
    __global__ void ScatterMassesInMemory(int *address, double *mass, int size, Particle *ParitcleList);
    void flush_memory_pos();
    void flush_memory_mass();
    void calculate_expansion();
    void gravity(int ni_tot, int ni_cor, int index[], double xi[][3], double acc[][3], double jerk[][3], double pot[]);
}
using namespace grapite;

void grapite::initialize_globals(int N, int n_loc) {
    t_exp = 0;
    Translator_j = TranslatorClass(n_loc);
    Translator_i = TranslatorClass(N);

    status = new int[N];

    PositionBuffer = MemoryBufferClass<3>(2*N, N); // Minimum max size is N+1 but please use much larger to prevent frequent reorganization!
    MassBuffer = MemoryBufferClass<1>(2*N, N);

    cudaMalloc((void**)&Particles_j, n_loc*sizeof(Particle));
    cudaMalloc((void**)&Particles_i, N*sizeof(Particle));
    Particles_ih = new Particle[N];

    SCFObject_cor.Init(N, 0, 0, 0, 0);
    SCFObject_hal.Init(N, 0, 0, 0, 0);

    cudaMalloc((void**)&Potential, N*sizeof(double));
    cudaMalloc((void**)&F, N*sizeof(vec3));
    Potential_h = new double[N];
    F_h = new vec3[N];
    JerkCalculator = JerkCalculatorClass(N);

    index_parted = new int[N];
    xi_parted    = (double (*)[3])(new double[N*3]);
    vi_parted    = (double (*)[3])(new double[N*3]);
    acc_parted   = (double (*)[3])(new double[N*3]);
    jerk_parted  = (double (*)[3])(new double[N*3]);
    pot_parted   = new double[N];
}

void grapite::generate_angular_momentum_mask(int N, double x[][3], double v[][3], int skip, double fraction, int output[]) {
    assert(N >= skip);
    double *L        = new double[N-skip];
    double *L_sorted = new double[N-skip];
    int N2 = (int)round((N-skip)*fraction);
    int N1 = N - N2;

    for (int j = skip; j < N; j++) {
        double x_ = x[j][0];
        double y_ = x[j][1];
        double z_ = x[j][2];
        double vx_ = v[j][0];
        double vy_ = v[j][1];
        double vz_ = v[j][2];
        L[j-skip] = sqrt(pow(y_*vz_-z_*vy_, 2) + pow(z_*vx_-x_*vz_, 2) + pow(x_*vy_-y_*vx_, 2));
    }
    std::copy(L, L+N-skip, L_sorted);
    thrust::sort(L_sorted, L_sorted + N-skip);
    double L_crit;
    if (fraction!=0)
        L_crit = 0.5*(L_sorted[N1-skip-1] + L_sorted[N1-skip]);
    else {
        N1 = N;
        L_crit = 2*L_sorted[N1-skip-1];
    }
    for (int j = 0; j < skip; j++) {
        output[j] = 0;
    }
    for (int j = skip; j < N; j++) {
        output[j] = L[j-skip] > L_crit;
    }
}

void grapite::set_j_particle_pos(int j, double *x) {
    PositionBuffer.add_particle(j, x);
}

void grapite::set_j_particle_mass(int j, double m) {
    MassBuffer.add_particle(j, &m);
}

int grapite::partition_arrays_i(int n_act, int ind_act[], double x_act_array[][3], double v_act_array[][3], int ind_act_parted[], double x_act_parted[][3], double v_act_parted[][3]) {
    for (int i = 0; i < n_act; i++) {
        Translator_i.AddParticle(i, status[ind_act[i]]);
    }
    int N1 = Translator_i.Partition(n_act);
    for (int i = 0; i < n_act; i++) {
        int parted_index = Translator_i.PartedToOrig(i);
        #pragma unroll
        for (int k = 0; k < 3; k++) {
            x_act_parted[i][k] = x_act_array[parted_index][k];
            v_act_parted[i][k] = v_act_array[parted_index][k];
        }
        ind_act_parted[i] = ind_act[parted_index];
    }

    return N1;
}

void grapite::revert_arrays_i(double acc_parted[][3], double jerk_parted[][3], double pot_parted[], double acc[][3], double jerk[][3], double pot[]) {
    for (int i = 0; i < ni_total; i++) {
        #pragma unroll
        for (int k = 0; k < 3; k++) {
            acc[Translator_i.PartedToOrig(i)][k] = acc_parted[i][k];
            jerk[Translator_i.PartedToOrig(i)][k] = jerk_parted[i][k];
        }
        pot[Translator_i.PartedToOrig(i)] = pot_parted[i];
    }
}

__global__ void grapite::ScatterPositionsInMemory(int *address, double *pos, int size, Particle *ParitcleList) {
    int i = threadIdx.x + blockIdx.x *  blockDim.x;
    while (i < size) {
        ParitcleList[address[i]].pos.x = pos[i*3];
        ParitcleList[address[i]].pos.y = pos[i*3+1];
        ParitcleList[address[i]].pos.z = pos[i*3+2];

        i += blockDim.x * gridDim.x;
    }
}

__global__ void grapite::ScatterMassesInMemory(int *address, double *mass, int size, Particle *ParitcleList) {
    int i = threadIdx.x + blockIdx.x *  blockDim.x;
    while (i < size) {
        ParitcleList[address[i]].m = mass[i];

        i += blockDim.x * gridDim.x;
    }
}

void grapite::flush_memory_pos() {
    PositionBuffer.h2d();
    ScatterPositionsInMemory<<<128,128>>>(PositionBuffer.buffer_address_gpu, PositionBuffer.buffer_prop_gpu, PositionBuffer.size, Particles_j);
    PositionBuffer.reset();
}

void grapite::flush_memory_mass() {
    MassBuffer.h2d();
    ScatterMassesInMemory<<<128,128>>>(MassBuffer.buffer_address_gpu, MassBuffer.buffer_prop_gpu, MassBuffer.size, Particles_j);
    MassBuffer.reset();
}

void grapite::calculate_expansion() {
    flush_memory_pos();
    flush_memory_mass();

    SCFObject_cor.GetGpuLock();
    SCFObject_cor.SendCachePointersToGPU();
    SCFObject_cor.LoadParticlesToCache(Particles_j, nj_core);
    SCFObject_cor.CalculateCoefficients();
    SCFObject_cor.ReleaseGpuLock();

    SCFObject_hal.GetGpuLock();
    SCFObject_hal.SendCachePointersToGPU();
    SCFObject_hal.LoadParticlesToCache(Particles_j+nj_core, nj_total-nj_core);
    SCFObject_hal.CalculateCoefficients();
    SCFObject_hal.ReleaseGpuLock();
}

void grapite::gravity(int ni_tot, int ni_cor, int index[], double xi[][3], double acc[][3], double jerk[][3], double pot[]) {
    // Put i-particles on the GPU in the ETICS data format
    // maybe should be a new subroutine?
    for (int i = 0; i < ni_tot; i++) {
        Particle p;
        p.m = 1;
        p.pos = vec3(xi[i][0], xi[i][1], xi[i][2]);
        Particles_ih[i] = p;
    }
    cudaMemcpy(Particles_i, Particles_ih, ni_tot*sizeof(Particle), cudaMemcpyHostToDevice);
    int ni_hal = ni_tot - ni_cor;

    // Calculate force and potential on i particles from the halo distribution
    SCFObject_hal.GetGpuLock();
    SCFObject_hal.SendCoeffsToGPU();
    SCFObject_hal.SendCachePointersToGPU();
    SCFObject_hal.LoadParticlesToCache(Particles_i, ni_tot);
    SCFObject_hal.CalculateGravityFromCoefficients(Potential, F);
    SCFObject_hal.ReleaseGpuLock();

    // Move potential and force to CPU
    cudaMemcpy(Potential_h, Potential, ni_tot*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(F_h, F, ni_tot*sizeof(vec3), cudaMemcpyDeviceToHost);
    // Add potential and force to the relevant arrays
    for (int i = 0; i < ni_cor; i++) {
        // halo --> core
        pot[i]  += Potential_h[i];
        acc[i][0] += F_h[i].x;
        acc[i][1] += F_h[i].y;
        acc[i][2] += F_h[i].z;
    }
    for (int i = ni_cor; i < ni_tot; i++) { // Halo particles were not previously initialized, so += changes to =
        // halo --> halo
        pot[i]  = Potential_h[i];
        acc[i][0] = F_h[i].x;
        acc[i][1] = F_h[i].y;
        acc[i][2] = F_h[i].z;
        jerk[i][0] = 0;
        jerk[i][1] = 0;
        jerk[i][2] = 0;
    }

    for (int i = 0; i < ni_tot; i++) JerkCalculator.load(index[i], ti, F_h[i]);

    if (ni_hal > 0) {
        // core --> halo (obviously, only needed if there are halo particles in the batch)
        SCFObject_cor.GetGpuLock();
        SCFObject_cor.SendCoeffsToGPU();
        SCFObject_cor.SendCachePointersToGPU();
        SCFObject_cor.LoadParticlesToCache(Particles_i+ni_cor, ni_hal);
        SCFObject_cor.CalculateGravityFromCoefficients(Potential, F);
        SCFObject_cor.ReleaseGpuLock();

        cudaMemcpy(Potential_h, Potential, ni_hal*sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(F_h, F, ni_hal*sizeof(vec3), cudaMemcpyDeviceToHost);
        for (int i = ni_cor; i < ni_tot; i++) {
            // core --> halo
            pot[i] += Potential_h[i-ni_cor];
            acc[i][0] += F_h[i-ni_cor].x;
            acc[i][1] += F_h[i-ni_cor].y;
            acc[i][2] += F_h[i-ni_cor].z;
        }

        for (int i = ni_cor; i < ni_tot; i++) JerkCalculator.update_add(index[i], F_h[i-ni_cor]);

    }
    for (int i = 0; i < ni_tot; i++) JerkCalculator.add_jerk(index[i], jerk[i]);
}

extern "C"
void g6_open(int clusterid) {
    printf("=== This is Grapite (intermediate layer) ===\n");
    __grape_g6_open(clusterid);
    // Only the GRAPE is opened, the ETICS part will be initialized when we know how many particles we have.
}

extern "C"
void g6_set_ti(int clusterid, double ti) {
    __grape_g6_set_ti(clusterid, ti);
    grapite::ti = ti;
}

extern "C"
int grapite_tag_particles(int N, double x[][3], double v[][3], int skip, double fraction, int rank, int n_loc) {
    grapite::rank = rank;
    initialize_globals(N, n_loc);
    nj_total = n_loc;
    generate_angular_momentum_mask(N, x, v, skip, fraction, status);
    int j_start = n_loc*rank;
    int j_end   = n_loc*(rank+1);
    for (int j = j_start; j < j_end; j++) Translator_j.AddParticle(j-j_start, status[j]);
    nj_core = Translator_j.Partition(n_loc);
    return nj_core;
}

#warning we do not differentiate between address and index here
extern "C"
int g6_set_j_particle(int clusterid, int address, int index, double tj, double dtj, double mass,
                      double a2by18[3], double a1by6[3], double aby2[3], double v[3], double x[3]) {
    int new_address = Translator_j.OrigToParted(address);
    if (new_address < nj_core) __grape_g6_set_j_particle(clusterid, new_address, index, tj, dtj, mass, a2by18, a1by6, aby2, v, x);



    set_j_particle_mass(new_address, mass);
    set_j_particle_pos(new_address, x);
    return 0; // we don't know if there is any meaning to the returned value.
}

extern "C"
void g6calc_firsthalf(int clusterid, int nj, int ni, int index[], double xi[][3], double vi[][3],
                      double fold[][3], double jold[][3], double phiold[], double eps2, double h2[]) {
    ni_total = ni;
    ni_core = partition_arrays_i(ni, index, xi, vi, index_parted, xi_parted, vi_parted);
    __grape_g6calc_firsthalf(clusterid, nj_core, ni_core, index_parted, xi_parted, vi_parted, fold, jold, phiold, eps2, h2);
    if (ti >= t_exp) {
        calculate_expansion();
        t_exp += SCFSTEP;
    }
    /// NOTICE! we don't partition fold, jold, phiold -- is it OK?? Those values may be used for predict of i particles... but we actually predict them outside; xi and vi are already the predicted values, so what gives?
    // if fold, jold, phiold are indeed ignored by the grape, or at least yebisu and sapporo, we should pass the grapes zero filled arrays, not the un-partitioned arrays we get from phiGRAPE
}

extern "C"
int g6calc_lasthalf(int clusterid, int nj, int ni, int index[], double xi[][3], double vi[][3],
                    double eps2, double h2[], double acc[][3], double jerk[][3], double pot[]) {
    __grape_g6calc_lasthalf(clusterid, nj_core, ni_core, index_parted, xi_parted, vi_parted, eps2, h2, acc_parted, jerk_parted, pot_parted);
    gravity(ni_total, ni_core, index_parted, xi_parted, acc_parted, jerk_parted, pot_parted);
    revert_arrays_i(acc_parted, jerk_parted, pot_parted, acc, jerk, pot);
    return 0;
}
