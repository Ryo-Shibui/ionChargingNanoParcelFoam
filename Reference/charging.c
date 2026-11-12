#include "charging.h"

#include "constants.h"
#include "util.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

// -----------------------------------------------------------------------------
// Collision coefficient beta with ion velocity distribution.
//
// Model selection:
//   PsiE > 0  : Coulomb attraction. Use corrected Eq. (4-37)
//               H = exp(mu) * H_HS.
//   PsiE <= 0 : neutral or Coulomb-repulsive interaction. Retain the legacy
//               piecewise H model (HS / Eq.18 / free-molecular).
//
// Reported applicability of Eq. (4-37):
//   0 < PsiE < 60 and KnD < 2000.
// Outside this range, Eq. (4-37) is still used when it is mathematically real,
// but the result is marked EXTRAPOLATION in eta_table.csv.
// A legacy fallback is used only if Eq. (4-37) cannot be evaluated as a finite
// real number, so the simulation does not silently receive beta=0.
//
// The Eq. (4-37) implementation follows the executable cells and numerical
// output of the supplied notebook:
//   y  = base^(-1/k)
//   mu = (C/A) * base^(-1/k - 1) * exp(-y)
// and B uses ln(1 + 1.537*PsiE).
// -----------------------------------------------------------------------------

enum {
    ETA_C_NGAUSS = 800,
    ETA_FM_NLAG  = 80,
    BCRIT_NR     = 50000
};

static const double B1_EQ18 = 1.598;
static const double B2_EQ18 = 1.1709;

// -----------------------------
// 1) Gauss-Legendre nodes/weights on [-1,1]
//    Standard Newton iteration on Legendre polynomial roots.
// -----------------------------
static void gauss_legendre(int n, double *x, double *w)
{
    // Based on classic algorithm (e.g., Numerical Recipes).
    const double eps = 1e-15;
    const int m = (n + 1) / 2;

    for (int i = 0; i < m; i++) {
        // Initial guess
        double z = cos(M_PI * ((double)i + 0.75) / ((double)n + 0.5));
        double z1;

        // Newton iteration
        for (int it = 0; it < 100; it++) {
            double p1 = 1.0;
            double p2 = 0.0;
            for (int j = 1; j <= n; j++) {
                double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * (double)j - 1.0) * z * p2 - ((double)j - 1.0) * p3) / (double)j;
            }
            // p1 = P_n(z), p2 = P_{n-1}(z)
            double pp = ((double)n) * (z * p1 - p2) / (z * z - 1.0); // P'_n(z)
            z1 = z;
            z = z1 - p1 / pp;
            if (fabs(z - z1) < eps) break;
        }

        // Recompute derivative for weight
        double p1 = 1.0;
        double p2 = 0.0;
        for (int j = 1; j <= n; j++) {
            double p3 = p2;
            p2 = p1;
            p1 = ((2.0 * (double)j - 1.0) * z * p2 - ((double)j - 1.0) * p3) / (double)j;
        }
        double pp = ((double)n) * (z * p1 - p2) / (z * z - 1.0);
        double wi = 2.0 / ((1.0 - z * z) * pp * pp);

        x[i] = -z;
        x[n - 1 - i] = z;
        w[i] = wi;
        w[n - 1 - i] = wi;
    }
}

// -----------------------------
// 2) Gauss-Laguerre nodes/weights for alpha=0 on [0,inf)
//    We compute roots of L_n(x) via Newton with initial guesses (NR-style)
//    and weights using: w_i = x_i / ((n+1)^2 * (L_{n+1}(x_i))^2)
// -----------------------------
static void laguerre_eval_Ln_Lnm1(int n, double x, double *Ln, double *Lnm1)
{
    // Recurrence for ordinary Laguerre L_k(x):
    // L_0=1, L_1=1-x
    // (k+1)L_{k+1} = (2k+1-x)L_k - k L_{k-1}
    double L0 = 1.0;
    if (n == 0) {
        *Ln = L0;
        *Lnm1 = 0.0;
        return;
    }
    double L1 = 1.0 - x;
    if (n == 1) {
        *Ln = L1;
        *Lnm1 = L0;
        return;
    }
    for (int k = 1; k < n; k++) {
        double L2 = ((2.0 * (double)k + 1.0 - x) * L1 - (double)k * L0) / ((double)k + 1.0);
        L0 = L1;
        L1 = L2;
    }
    *Ln = L1;
    *Lnm1 = L0;
}

static double laguerre_Ln(int n, double x)
{
    double Ln, Lnm1;
    laguerre_eval_Ln_Lnm1(n, x, &Ln, &Lnm1);
    return Ln;
}

static void gauss_laguerre_alpha0(int n, double *x, double *w)
{
    const double eps = 1e-14;

    for (int i = 1; i <= n; i++) {
        // Initial guesses (Numerical Recipes style)
        // NOTE: Keep this explicit to avoid "may be used uninitialized" warnings.
        // The update rule is written in terms of the previously converged roots.
        double z = 0.0;
        if (i == 1) {
            z = 3.0 / (1.0 + 2.4 * (double)n);
        } else if (i == 2) {
            z = x[0] + 15.0 / (1.0 + 2.5 * (double)n);
        } else {
            double ai = (double)i - 2.0;
            z = x[i - 2] + ((1.0 + 2.55 * ai) / (1.9 * ai)) * (x[i - 2] - x[i - 3]);
        }

        // Newton iterations
        for (int it = 0; it < 50; it++) {
            double Ln, Lnm1;
            laguerre_eval_Ln_Lnm1(n, z, &Ln, &Lnm1);

            // derivative: x L_n'(x) = n(L_n - L_{n-1})
            // => L_n'(x) = n(L_n - L_{n-1})/x
            double dLn;
            if (z == 0.0) {
                dLn = - (double)n; // limiting slope near 0 (rare for roots)
            } else {
                dLn = ((double)n) * (Ln - Lnm1) / z;
            }

            double dz = Ln / dLn;
            z -= dz;
            if (fabs(dz) < eps) break;
        }

        x[i - 1] = z;

        // Weight using L_{n+1}(x_i)
        double Lnp1 = laguerre_Ln(n + 1, z);
        double denom = ((double)(n + 1) * (double)(n + 1)) * (Lnp1 * Lnp1);
        w[i - 1] = (denom > 0.0) ? (z / denom) : 0.0;
    }
}

// -----------------------------
// 3) Paper definitions: PsiE, PsiI, phi(rhat)
// -----------------------------
static double PsiE_paper(int zp, int zi, double T, double ap, double eps_r_g)
{
    const double eps_g = EPS0 * eps_r_g;
    return -((double)zp * (double)zi) * (E_CH * E_CH) / (4.0 * M_PI * eps_g * KB * T * ap);
}

static double PsiI_paper(int zi, double T, double ap, double eps_r_particle, double eps_r_g)
{
    double pref;
    if (!util_is_finite(eps_r_particle)) {
        pref = 1.0;
    } else {
        // Treat very large as conductor
        if (eps_r_particle > 1e20) pref = 1.0;
        else pref = (eps_r_particle - 1.0) / (eps_r_particle + 1.0);
    }
    const double eps_g = EPS0 * eps_r_g;
    return pref * ((double)zi * (double)zi) * (E_CH * E_CH) / (4.0 * M_PI * eps_g * KB * T * ap);
}

static double phi_paper(double rhat, double PsiE, double PsiI)
{
    // rhat = r/ap, must be > 1
    if (!(rhat > 1.0)) rhat = 1.0 + 1e-15;

    const double r2 = rhat * rhat;
    const double denom = 2.0 * r2 * (r2 - 1.0);
    if (!(denom > 0.0)) {
        return -1e300; // exp -> 0
    }
    return -PsiE / rhat - PsiI / denom;
}

// -----------------------------
// 4) eta_c from Eq.(3b): eta_c = 1 / ∫_0^1 exp(phi(1/u)) du
// -----------------------------
static double eta_c_from_phi(double PsiE, double PsiI)
{
    const int n = ETA_C_NGAUSS;
    double *x = (double*)malloc((size_t)n * sizeof(double));
    double *w = (double*)malloc((size_t)n * sizeof(double));
    if (!x || !w) {
        free(x);
        free(w);
        return 0.0;
    }

    gauss_legendre(n, x, w);

    double integral = 0.0;
    for (int i = 0; i < n; i++) {
        // Map [-1,1] -> [0,1]
        const double u = 0.5 * (x[i] + 1.0);
        const double wu = 0.5 * w[i];
        if (u <= 0.0) continue;
        const double rhat = 1.0 / u;
        const double ph = phi_paper(rhat, PsiE, PsiI);
        integral += wu * util_exp_clip(ph);
    }

    free(x);
    free(w);

    if (!(integral > 0.0) || !util_is_finite(integral)) return 0.0;
    const double eta_c = 1.0 / integral;
    if (!util_is_finite(eta_c) || eta_c < 0.0) return 0.0;
    return eta_c;
}

// -----------------------------
// 5) eta_fm (Ouyang 2012) with bcrit cache
// -----------------------------
typedef struct {
    double PsiE;
    double PsiI;
    int n_r;
    double *r_grid;   // rhat grid
    double *phi_grid; // phi(rhat)
} BcritCache;

static void bcrit_cache_free(BcritCache *c)
{
    if (!c) return;
    free(c->r_grid);
    free(c->phi_grid);
    c->r_grid = NULL;
    c->phi_grid = NULL;
    c->n_r = 0;
}

static int bcrit_cache_init(BcritCache *c, double PsiE, double PsiI, double rmax, int n_r)
{
    if (!c) return -1;
    c->PsiE = PsiE;
    c->PsiI = PsiI;
    c->n_r = n_r;
    c->r_grid = (double*)malloc((size_t)n_r * sizeof(double));
    c->phi_grid = (double*)malloc((size_t)n_r * sizeof(double));
    if (!c->r_grid || !c->phi_grid) {
        bcrit_cache_free(c);
        return -1;
    }

    const double r1 = 1.0 + 1e-10;
    const double log_r1 = log(r1);
    const double log_rmax = log(rmax);
    const double denom = (double)(n_r - 1);

    for (int i = 0; i < n_r; i++) {
        double t = (denom > 0.0) ? ((double)i / denom) : 0.0;
        double lr = log_r1 + t * (log_rmax - log_r1);
        double rhat = exp(lr);
        c->r_grid[i] = rhat;
        c->phi_grid[i] = phi_paper(rhat, PsiE, PsiI);
    }
    return 0;
}

static double b_of_r(double rhat, double nu, double PsiE, double PsiI)
{
    const double ph = phi_paper(rhat, PsiE, PsiI);
    const double inside = 1.0 + fabs(ph) / (nu * nu);
    if (!(inside > 0.0) || !util_is_finite(inside)) return NAN;
    const double val = rhat * sqrt(inside);
    if (!util_is_finite(val)) return NAN;
    return val;
}

static double golden_section_minimize_b(double a, double b, double nu, double PsiE, double PsiI)
{
    // returns f(xopt)
    const double gr = (sqrt(5.0) - 1.0) / 2.0;
    double c = b - gr * (b - a);
    double d = a + gr * (b - a);
    double fc = b_of_r(c, nu, PsiE, PsiI);
    double fd = b_of_r(d, nu, PsiE, PsiI);

    for (int it = 0; it < 80; it++) {
        if (!util_is_finite(fc) && !util_is_finite(fd)) break;
        if (!util_is_finite(fd) || (util_is_finite(fc) && fc < fd)) {
            b = d;
            d = c;
            fd = fc;
            c = b - gr * (b - a);
            fc = b_of_r(c, nu, PsiE, PsiI);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + gr * (b - a);
            fd = b_of_r(d, nu, PsiE, PsiI);
        }
    }

    double xopt = 0.5 * (a + b);
    double fopt = b_of_r(xopt, nu, PsiE, PsiI);
    return util_is_finite(fopt) ? fopt : NAN;
}

static double bcrit_nu(const BcritCache *c, double nu)
{
    if (!c || nu <= 0.0) return 0.0;
    double best = INFINITY;
    int best_i = -1;

    const double inv_nu2 = 1.0 / (nu * nu);
    for (int i = 0; i < c->n_r; i++) {
        double inside = 1.0 + fabs(c->phi_grid[i]) * inv_nu2;
        if (!(inside > 0.0) || !util_is_finite(inside)) continue;
        double bval = c->r_grid[i] * sqrt(inside);
        if (!util_is_finite(bval)) continue;
        if (bval < best) {
            best = bval;
            best_i = i;
        }
    }

    if (best_i < 0 || !util_is_finite(best)) return 0.0;

    // Refine around best index using golden section on [j-2, j+2]
    int j = best_i;
    if (j < 2) j = 2;
    if (j > c->n_r - 3) j = c->n_r - 3;
    double a = c->r_grid[j - 2];
    double b = c->r_grid[j + 2];
    if (!(b > a)) return best;

    double refined = golden_section_minimize_b(a, b, nu, c->PsiE, c->PsiI);
    if (!util_is_finite(refined)) refined = best;
    return refined;
}

static double eta_fm_ouyang(double PsiE, double PsiI)
{
    // Build cache
    BcritCache cache = {0};
    if (bcrit_cache_init(&cache, PsiE, PsiI, 5e3, BCRIT_NR) != 0) return 0.0;

    // Gauss-Laguerre nodes/weights
    const int n = ETA_FM_NLAG;
    double x[n];
    double w[n];
    gauss_laguerre_alpha0(n, x, w);

    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        double xi = x[i];
        if (!(xi > 0.0)) continue;
        double nu = sqrt(xi);
        double bc = bcrit_nu(&cache, nu);
        sum += w[i] * (xi * bc * bc);
    }

    bcrit_cache_free(&cache);

    if (!util_is_finite(sum) || sum < 0.0) return 0.0;
    return sum;
}

// -----------------------------
// 6) H models
// -----------------------------
static double H_HS(double KnD)
{
    const double num = 4.0 * M_PI * KnD * KnD
        + 25.836 * KnD * KnD * KnD
        + sqrt(8.0 * M_PI) * 11.211 * KnD * KnD * KnD * KnD;
    const double den = 1.0 + 3.502 * KnD + 7.211 * KnD * KnD + 11.211 * KnD * KnD * KnD;
    if (!(den > 0.0)) return 0.0;
    return num / den;
}

static double H_continuum(double KnD)
{
    return 4.0 * M_PI * KnD * KnD;
}

static double H_free_molecular(double KnD)
{
    return sqrt(8.0 * M_PI) * KnD;
}

static double KnPsi_from_KnD(double KnD, double PsiE)
{
    const double psi = fabs(PsiE);
    if (psi == 0.0) return INFINITY;
    return 3.0 * KnD / (2.0 * psi);
}

static double H_eq18(double KnD, double KnPsi)
{
    double m = (KnD < KnPsi) ? KnD : KnPsi;
    return H_continuum(KnD) / (1.0 + B1_EQ18 * pow(m, B2_EQ18));
}

static double H_piecewise_legacy(double KnD, double PsiE)
{
    const double psi = fabs(PsiE);
    const double KnPsi = KnPsi_from_KnD(KnD, PsiE);

    if (KnD <= KnPsi / 3.0) {
        return H_HS(KnD);
    }

    if (psi >= 0.5 && KnPsi < 2.5 && KnD < 2.5 && KnD > KnPsi / 3.0) {
        return H_eq18(KnD, KnPsi);
    }

    if (KnPsi >= 2.5 && KnD >= 2.5 && KnD > KnPsi / 3.0) {
        return H_free_molecular(KnD);
    }

    return H_HS(KnD);
}

// -----------------------------
// 7) Eq. (4-37) correction for PsiE > 0
// -----------------------------
static int correction_parameters_437(double PsiE,
                                     double *A, double *B,
                                     double *C, double *k)
{
    if (A) *A = 2.5;
    if (B) *B = NAN;
    if (C) *C = NAN;
    if (k) *k = NAN;

    if (!(PsiE > 0.0) || !util_is_finite(PsiE)) return -1;

    const double a = 2.5;
    const double b = 4.528 * exp(-1.088 * PsiE)
                   + 0.7091 * log(1.0 + 1.537 * PsiE);
    const double c = 11.36 * pow(PsiE, 0.272) - 10.33;
    const double kk = -0.003533 * PsiE + 0.05971;

    if (!util_is_finite(b) || !util_is_finite(c) || !util_is_finite(kk)) return -1;

    if (A) *A = a;
    if (B) *B = b;
    if (C) *C = c;
    if (k) *k = kk;
    return 0;
}

static int mu_correction_437(double KnD, double PsiE,
                             double *A_out, double *B_out,
                             double *C_out, double *k_out,
                             double *mu_out)
{
    double A = NAN, B = NAN, C = NAN, k = NAN;
    if (A_out) *A_out = NAN;
    if (B_out) *B_out = NAN;
    if (C_out) *C_out = NAN;
    if (k_out) *k_out = NAN;
    if (mu_out) *mu_out = NAN;

    if (!(KnD > 0.0) || !util_is_finite(KnD)) return -1;
    if (correction_parameters_437(PsiE, &A, &B, &C, &k) != 0) return -1;

    if (A_out) *A_out = A;
    if (B_out) *B_out = B;
    if (C_out) *C_out = C;
    if (k_out) *k_out = k;

    const double x = (log(KnD) - B) / A;
    double mu;

    if (fabs(k) < 1e-12) {
        // lim_{k->0} (1+k*x)^(-1/k) = exp(-x)
        const double y = exp(-x);
        mu = (C / A) * y * exp(-y);
    } else {
        const double base = 1.0 + k * x;
        if (!(base > 0.0) || !util_is_finite(base)) return -1;

        const double y = pow(base, -1.0 / k);
        const double factor = pow(base, -1.0 / k - 1.0);
        if (!util_is_finite(y) || !util_is_finite(factor)) return -1;

        mu = (C / A) * factor * exp(-y);
    }

    if (!util_is_finite(mu)) return -1;
    if (mu_out) *mu_out = mu;
    return 0;
}

static int H_corrected_437(double KnD, double PsiE,
                           double *A, double *B, double *C, double *k,
                           double *mu, double *H_out)
{
    if (H_out) *H_out = NAN;

    double mu_local = NAN;
    if (mu_correction_437(KnD, PsiE, A, B, C, k, &mu_local) != 0) {
        if (mu) *mu = NAN;
        return -1;
    }
    if (mu) *mu = mu_local;

    const double Hhs = H_HS(KnD);
    if (!(Hhs > 0.0) || !util_is_finite(Hhs)) return -1;

    // Evaluate through log(H) to catch overflow explicitly.
    const double logH = log(Hhs) + mu_local;
    if (!util_is_finite(logH) || logH > log(DBL_MAX)) return -1;

    const double H = (logH < log(DBL_MIN)) ? 0.0 : exp(logH);
    if (!util_is_finite(H) || H < 0.0) return -1;

    if (H_out) *H_out = H;
    return 0;
}

// -----------------------------
// 8) KnD and beta (eta)
// -----------------------------
static double friction_from_ion_mobility(int zi, double Zi)
{
    return (fabs((double)zi) * E_CH) / Zi;
}

static double KnD_particle_ion(double T, double mi, double fi, double ap, double eta_c, double eta_fm)
{
    if (!(fi > 0.0) || !(ap > 0.0) || !(eta_fm > 0.0)) return 0.0;
    return sqrt(mi * KB * T) * eta_c / (fi * ap * eta_fm);
}

static double beta_from_H(double H, double mi, double fi, double ap, double eta_c, double eta_fm)
{
    if (!(mi > 0.0) || !(eta_c > 0.0)) return 0.0;
    return H * fi * ap * ap * ap * (eta_fm * eta_fm) / (mi * eta_c);
}

static double compute_beta_one(int zp, int zi, double T, double ap,
                               double eps_r_particle, double Zi, double mi,
                               CollisionDiagnostic *diag)
{
    CollisionDiagnostic d;
    memset(&d, 0, sizeof(d));
    d.A = d.B = d.C = d.k = d.mu = NAN;
    d.H = d.beta = 0.0;
    d.model_used = COLLISION_MODEL_LEGACY_PIECEWISE;

    const double eps_r_g = 1.0;
    d.PsiE = PsiE_paper(zp, zi, T, ap, eps_r_g);
    d.PsiI = PsiI_paper(zi, T, ap, eps_r_particle, eps_r_g);

    d.eta_c = eta_c_from_phi(d.PsiE, d.PsiI);
    if (!(d.eta_c > 0.0)) {
        if (diag) *diag = d;
        return 0.0;
    }

    d.eta_fm = eta_fm_ouyang(d.PsiE, d.PsiI);
    if (!(d.eta_fm > 0.0)) {
        if (diag) *diag = d;
        return 0.0;
    }

    const double fi = friction_from_ion_mobility(zi, Zi);
    d.KnD = KnD_particle_ion(T, mi, fi, ap, d.eta_c, d.eta_fm);
    if (!(d.KnD > 0.0)) {
        if (diag) *diag = d;
        return 0.0;
    }

    d.H_HS = H_HS(d.KnD);
    d.eq437_requested = (d.PsiE > 0.0) ? 1 : 0;
    d.psiE_in_range = (d.PsiE > 0.0 && d.PsiE < 60.0) ? 1 : 0;
    d.KnD_in_range = (d.KnD < 2000.0) ? 1 : 0;
    d.eq437_valid_range = (d.psiE_in_range && d.KnD_in_range) ? 1 : 0;

    if (d.eq437_requested) {
        double H437 = NAN;
        if (H_corrected_437(d.KnD, d.PsiE,
                            &d.A, &d.B, &d.C, &d.k, &d.mu, &H437) == 0) {
            d.H = H437;
            d.eq437_formula_real = 1;
            d.model_used = COLLISION_MODEL_EQ437;
        } else {
            // Numerical safety only. This path is explicitly exposed in CSV.
            d.H = H_piecewise_legacy(d.KnD, d.PsiE);
            d.eq437_formula_real = 0;
            d.model_used = COLLISION_MODEL_EQ437_FALLBACK_LEGACY;
        }
    } else {
        d.H = H_piecewise_legacy(d.KnD, d.PsiE);
        d.eq437_formula_real = 0;
        d.model_used = COLLISION_MODEL_LEGACY_PIECEWISE;
    }

    if (!(d.H >= 0.0) || !util_is_finite(d.H)) {
        if (diag) *diag = d;
        return 0.0;
    }

    d.beta = beta_from_H(d.H, mi, fi, ap, d.eta_c, d.eta_fm);
    if (!util_is_finite(d.beta) || d.beta < 0.0) d.beta = 0.0;

    if (diag) *diag = d;
    return d.beta;
}

// -----------------------------
// Public API
// -----------------------------
const char *collision_model_name(CollisionModel model)
{
    switch (model) {
        case COLLISION_MODEL_EQ437:
            return "EQ437_ATTRACTION";
        case COLLISION_MODEL_EQ437_FALLBACK_LEGACY:
            return "EQ437_FAILED_LEGACY_FALLBACK";
        case COLLISION_MODEL_LEGACY_PIECEWISE:
        default:
            return "LEGACY_PIECEWISE";
    }
}

const char *collision_validity_name(const CollisionDiagnostic *d)
{
    if (!d) return "UNKNOWN";
    if (!(d->PsiE > 0.0)) return "NOT_APPLICABLE_PsiE_LE_0";
    if (d->model_used == COLLISION_MODEL_EQ437_FALLBACK_LEGACY || !d->eq437_formula_real) {
        return "EQ437_NOT_REAL_FALLBACK";
    }
    if (d->eq437_valid_range) return "VALID";
    if (!d->psiE_in_range && !d->KnD_in_range) return "EXTRAPOLATION_PsiE_AND_KnD";
    if (!d->psiE_in_range) return "EXTRAPOLATION_PsiE";
    if (!d->KnD_in_range) return "EXTRAPOLATION_KnD";
    return "EXTRAPOLATION";
}

void eta_table_free(EtaTable *t)
{
    if (!t) return;
    free(t->eta_minus);
    free(t->eta_plus);
    free(t->diag_minus);
    free(t->diag_plus);
    t->eta_minus = NULL;
    t->eta_plus = NULL;
    t->diag_minus = NULL;
    t->diag_plus = NULL;
    t->m = 0;
}

int eta_table_build(EtaTable *t,
                    double dp,
                    double T,
                    int qmin, int qmax,
                    double eps_r_particle,
                    double ion_mobility_minus,
                    double ion_mobility_plus)
{
    if (!t) return -1;
    if (!(dp > 0.0) || !(T > 0.0)) return -1;
    if (qmax < qmin) return -1;
    if (!(ion_mobility_minus > 0.0) || !(ion_mobility_plus > 0.0)) return -1;

    const int m = qmax - qmin + 1;
    t->qmin = qmin;
    t->qmax = qmax;
    t->m = m;

    t->eta_minus = (double*)calloc((size_t)m, sizeof(double));
    t->eta_plus  = (double*)calloc((size_t)m, sizeof(double));
    t->diag_minus = (CollisionDiagnostic*)calloc((size_t)m, sizeof(CollisionDiagnostic));
    t->diag_plus  = (CollisionDiagnostic*)calloc((size_t)m, sizeof(CollisionDiagnostic));
    if (!t->eta_minus || !t->eta_plus || !t->diag_minus || !t->diag_plus) {
        eta_table_free(t);
        return -1;
    }

    const double ap = 0.5 * dp;
    const double mi = M_ION / NA_CONST; // ion mass [kg] (same convention as Python)

    for (int i = 0; i < m; i++) {
        const int zp = qmin + i;

        // Stored as eta_*(q), but the quantity is beta [m^3/s].
        t->eta_minus[i] = compute_beta_one(
            zp, -1, T, ap, eps_r_particle, ion_mobility_minus, mi,
            &t->diag_minus[i]);
        t->eta_plus[i] = compute_beta_one(
            zp, +1, T, ap, eps_r_particle, ion_mobility_plus, mi,
            &t->diag_plus[i]);
    }

    return 0;
}

double eta_lam_max(const EtaTable *t, double n_minus, double n_plus)
{
    if (!t || !t->eta_minus || !t->eta_plus) return 0.0;
    double best = 0.0;
    for (int i = 0; i < t->m; i++) {
        const double lam = t->eta_minus[i] * n_minus + t->eta_plus[i] * n_plus;
        if (lam > best) best = lam;
    }
    return best;
}
