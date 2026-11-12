# Collision-coefficient model update

The collision coefficient calculation is selected from the sign of `PsiE`:

- `PsiE > 0`: the particle and ion have opposite charge signs and Coulomb attraction acts. The code uses the corrected Eq. (4-37):

  `H = exp(mu) * H_HS`

- `PsiE <= 0`: the pair is neutral with respect to direct Coulomb interaction or Coulomb-repulsive. The previous piecewise model is retained.

The reported range of Eq. (4-37) is:

- `0 < PsiE < 60`
- `KnD < 2000`

When `PsiE > 0` but the point is outside that range, the code still uses Eq. (4-37) as an extrapolation and marks it in `eta_table.csv`. If Eq. (4-37) cannot be evaluated as a finite real number, the code uses the old piecewise model as an explicit safety fallback and reports this in the same CSV.

## Output

The first six columns of `eta_table.csv` are unchanged. Additional columns report, separately for negative and positive ions:

- model used
- validity status
- `PsiE`, `PsiI`, `eta_c`, `eta_fm`, `KnD`
- Eq. (4-37) parameters `A`, `B`, `C`, `k`, `mu`
- `H_HS`, selected `H`
- mathematical-real and applicability flags

Validity status values include:

- `VALID`
- `EXTRAPOLATION_PsiE`
- `EXTRAPOLATION_KnD`
- `EXTRAPOLATION_PsiE_AND_KnD`
- `NOT_APPLICABLE_PsiE_LE_0`
- `EQ437_NOT_REAL_FALLBACK`

The implementation follows the executable calculation and numerical output of the supplied notebook, including the negative exponent form and the coefficient `1.537` in `B`.
