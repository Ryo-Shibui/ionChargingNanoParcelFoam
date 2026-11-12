#ifndef CHARGING_H
#define CHARGING_H

typedef enum {
    COLLISION_MODEL_LEGACY_PIECEWISE = 0,
    COLLISION_MODEL_EQ437 = 1,
    COLLISION_MODEL_EQ437_FALLBACK_LEGACY = 2
} CollisionModel;

typedef struct {
    double PsiE;
    double PsiI;
    double eta_c;
    double eta_fm;
    double KnD;

    /* Eq. (4-37) parameters. NaN when PsiE <= 0. */
    double A;
    double B;
    double C;
    double k;
    double mu;

    double H_HS;
    double H;
    double beta;

    int eq437_requested;       /* PsiE > 0 */
    int eq437_formula_real;    /* Eq. (4-37) produced a finite real H */
    int psiE_in_range;         /* 0 < PsiE < 60 */
    int KnD_in_range;          /* KnD < 2000 */
    int eq437_valid_range;     /* both reported validity conditions */

    CollisionModel model_used;
} CollisionDiagnostic;

typedef struct {
    int qmin;
    int qmax;
    int m;
    double *eta_minus;
    double *eta_plus;
    CollisionDiagnostic *diag_minus;
    CollisionDiagnostic *diag_plus;
} EtaTable;

void eta_table_free(EtaTable *t);

int eta_table_build(EtaTable *t,
                    double dp,
                    double T,
                    int qmin, int qmax,
                    double eps_r_particle,
                    double ion_mobility_minus,
                    double ion_mobility_plus);

double eta_lam_max(const EtaTable *t, double n_minus, double n_plus);

const char *collision_model_name(CollisionModel model);
const char *collision_validity_name(const CollisionDiagnostic *d);

#endif
