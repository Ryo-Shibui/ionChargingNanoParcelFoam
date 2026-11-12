#include "util.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

int util_is_finite(double x) {
    return isfinite(x) ? 1 : 0;
}

double util_clamp(double x, double lo, double hi) {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

double util_exp_clip(double x) {
    if (x > 700.0) x = 700.0;
    if (x < -745.0) x = -745.0;
    return exp(x);
}

int util_ensure_dir(const char *path)
{
    if (!path || path[0] == '\0') return -1;

    struct stat st;
    if (stat(path, &st) == 0) {
        if (S_ISDIR(st.st_mode)) return 0;
        return -1;
    }

    if (mkdir(path, 0755) == 0) return 0;

    // simple recursive mkdir -p
    char tmp[512];
    size_t len = strlen(path);
    if (len >= sizeof(tmp)) return -1;
    strcpy(tmp, path);

    for (size_t i = 1; i < len; i++) {
        if (tmp[i] == '/') {
            tmp[i] = '\0';
            if (tmp[0] != '\0') (void)mkdir(tmp, 0755);
            tmp[i] = '/';
        }
    }

    if (mkdir(tmp, 0755) == 0) return 0;
    if (stat(path, &st) == 0 && S_ISDIR(st.st_mode)) return 0;
    return -1;
}

int util_snprintf(char *dst, size_t n, const char *fmt, ...)
{
    if (!dst || n == 0 || !fmt) return -1;
    va_list ap;
    va_start(ap, fmt);
    int ret = vsnprintf(dst, n, fmt, ap);
    va_end(ap);
    if (ret < 0) return -1;
    if ((size_t)ret >= n) return -1;
    return 0;
}
