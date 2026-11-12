#ifndef UTIL_H
#define UTIL_H

#include <stddef.h>

int util_is_finite(double x);
double util_clamp(double x, double lo, double hi);
double util_exp_clip(double x);   // exp(x) with overflow/underflow clipping

int util_ensure_dir(const char *path); // mkdir -p (0755)
int util_snprintf(char *dst, size_t n, const char *fmt, ...);

#endif
