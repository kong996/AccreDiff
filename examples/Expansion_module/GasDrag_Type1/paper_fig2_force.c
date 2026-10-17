/* Compiled REBOUND additional force for the paper-specific Fig. 2 model.
 *
 * Units: AU, Msun, G=1. Time is the REBOUND code time (year/(2*pi)).
 * This mirrors paper_fig2_drag.py paper_2016 and the precomputed radial
 * force in paper_fig2_gravity.py. The radial gas profile supports the paper
 * power law and Gas_DiskModel's logarithmically curved surface density.
 * The flow/drag convention remains paper_2016.
 *
 * The parameters hold borrowed pointers to NumPy arrays. The Python bridge
 * keeps those arrays and this shared library alive through integration.
 * Do not run concurrent simulations with different configurations in one
 * process: these parameters are process-global.
 */

#include <math.h>
#include "rebound.h"

struct fig2_params {
    int n_table;
    const double *table_r;
    const double *table_fr;
    const double *table_dfr;
    double tau_code;
    double sigma0_cgs;
    double alpha1;
    double alpha2;
    double h1_au;
    double height_exponent;
    double mu_code;
    double drag_coefficient;
    double solid_density_cgs;
    double radius_cm;
    double au_cm;
    int first_test;
    int use_disk_gravity;
    int use_drag;
    int use_boundaries;
    double inner_boundary;
    double outer_boundary;
    int require_positive_energy;
    int *active;
    int *event_code;
    double *event_time;
    double *event_radius;
    double *event_position;
    double *event_velocity;
    double *event_radial_velocity;
    double *event_kepler_energy;
};

static struct fig2_params params;
static int domain_violations = 0;
static int first_violation_particle = -1;
static int first_violation_kind = 0;
static double first_violation_time = NAN;
static double first_violation_radius = NAN;
static double first_violation_kepler_energy = NAN;
static double minimum_violation_radius = INFINITY;
static double maximum_violation_radius = -INFINITY;

static void record_violation(int particle, int kind, double time, double radius,
                             double kepler_energy) {
    if (domain_violations == 0) {
        first_violation_particle = particle;
        first_violation_kind = kind;
        first_violation_time = time;
        first_violation_radius = radius;
        first_violation_kepler_energy = kepler_energy;
    }
    domain_violations++;
    if (isfinite(radius)) {
        if (radius < minimum_violation_radius) minimum_violation_radius = radius;
        if (radius > maximum_violation_radius) maximum_violation_radius = radius;
    }
}

void fig2_set_params(const struct fig2_params *incoming) {
    params = *incoming;
    domain_violations = 0;
    first_violation_particle = -1;
    first_violation_kind = 0;
    first_violation_time = NAN;
    first_violation_radius = NAN;
    first_violation_kepler_energy = NAN;
    minimum_violation_radius = INFINITY;
    maximum_violation_radius = -INFINITY;
}

int fig2_get_domain_violations(void) {
    return domain_violations;
}

int fig2_get_first_violation_particle(void) { return first_violation_particle; }
int fig2_get_first_violation_kind(void) { return first_violation_kind; }
double fig2_get_first_violation_time(void) { return first_violation_time; }
double fig2_get_first_violation_radius(void) { return first_violation_radius; }
double fig2_get_first_violation_kepler_energy(void) { return first_violation_kepler_energy; }
double fig2_get_minimum_violation_radius(void) { return minimum_violation_radius; }
double fig2_get_maximum_violation_radius(void) { return maximum_violation_radius; }

static double radial_force(double r, int particle, double time, double kepler_energy) {
    if (r < params.table_r[0] || r > params.table_r[params.n_table-1]) {
        record_violation(particle, 1, time, r, kepler_energy);
        return 0.0;
    }
    int lo = 0;
    int hi = params.n_table-1;
    while (hi-lo > 1) {
        int mid = (lo+hi)/2;
        if (params.table_r[mid] <= r) lo = mid;
        else hi = mid;
    }
    const double h = params.table_r[lo+1]-params.table_r[lo];
    const double s = (r-params.table_r[lo])/h;
    const double s2 = s*s;
    const double s3 = s2*s;
    return (2*s3-3*s2+1)*params.table_fr[lo]
         + (s3-2*s2+s)*h*params.table_dfr[lo]
         + (-2*s3+3*s2)*params.table_fr[lo+1]
         + (s3-s2)*h*params.table_dfr[lo+1];
}

static void record_particle_event(int j, int code, double time, double radius,
                                  double x, double y, double z,
                                  double vx, double vy, double vz,
                                  double radial_velocity, double energy) {
    if (params.event_code[j] != 0) return;
    params.event_code[j] = code;
    params.event_time[j] = time;
    params.event_radius[j] = radius;
    params.event_position[3*j] = x;
    params.event_position[3*j+1] = y;
    params.event_position[3*j+2] = z;
    params.event_velocity[3*j] = vx;
    params.event_velocity[3*j+1] = vy;
    params.event_velocity[3*j+2] = vz;
    params.event_radial_velocity[j] = radial_velocity;
    params.event_kepler_energy[j] = energy;
    params.active[j] = 0;
}

void fig2_additional_forces(struct reb_simulation *sim) {
    struct reb_particle *sun = &sim->particles[0];
    const double decay = exp(-sim->t/params.tau_code);
    for (int j=1; j<sim->N; ++j) {
        struct reb_particle *p = &sim->particles[j];
        const double x = p->x-sun->x;
        const double y = p->y-sun->y;
        const double z = p->z-sun->z;
        const double r = hypot(x,y);
        const double vx = p->vx-sun->vx;
        const double vy = p->vy-sun->vy;
        const double vz = p->vz-sun->vz;
        if (params.use_boundaries && j >= params.first_test) {
            if (!params.active[j]) continue;
            if (!isfinite(r) || r <= 0.0 || !isfinite(vx) ||
                !isfinite(vy) || !isfinite(vz)) {
                record_particle_event(j, 4, sim->t, r, x, y, z, vx, vy, vz,
                                      NAN, NAN);
                continue;
            }
            const double energy = 0.5*(vx*vx+vy*vy+vz*vz)-sun->m/r;
            const double vr = (x*vx+y*vy+z*vz)/r;
            if (r <= params.inner_boundary) {
                record_particle_event(j, 1, sim->t, r, x, y, z, vx, vy, vz,
                                      vr, energy);
                continue;
            }
            if (r >= params.outer_boundary) {
                const int ejected = vr > 0.0 &&
                    (energy >= 0.0 || !params.require_positive_energy);
                record_particle_event(j, ejected ? 3 : 2, sim->t, r,
                                      x, y, z, vx, vy, vz, vr, energy);
                continue;
            }
        }
        if (!isfinite(r) || r <= 0.0 || !isfinite(vx) ||
            !isfinite(vy) || !isfinite(vz)) {
            record_violation(j, 2, sim->t, r, NAN);
            continue;
        }
        const double kepler_energy = 0.5*(vx*vx+vy*vy+vz*vz)-sun->m/r;
        if (params.use_disk_gravity) {
            const double fr = radial_force(r, j, sim->t, kepler_energy)*decay;
            p->ax += fr*x/r;
            p->ay += fr*y/r;
        }
        if (params.use_drag && j >= params.first_test) {
            const double h_au = params.h1_au*pow(r,params.height_exponent);
            const double logr = log(r);
            const double sigma = params.sigma0_cgs
                               *exp(-(params.alpha1*logr+params.alpha2)*logr)*decay;
            const double rho = sigma/(h_au*params.au_cm);
            const double eta = pow(h_au/(4.0*r),2);
            const double vg = sqrt(params.mu_code/r)*(1.0-eta);
            const double ux = (p->vx-sun->vx) + vg*y/r;
            const double uy = (p->vy-sun->vy) - vg*x/r;
            const double uabs = hypot(ux,uy);
            const double c = params.drag_coefficient*rho*params.au_cm
                           /(params.solid_density_cgs*params.radius_cm);
            p->ax -= c*uabs*ux;
            p->ay -= c*uabs*uy;
        }
    }
}
