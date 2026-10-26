/* ipoly.c -- Poly.eval of certificate_check.py in C, operation for operation.
 *
 * Every product and sum is taken in IEEE double (round to nearest) and then moved one
 * step outward (the next double towards -inf for a lower end, towards +inf for an upper
 * end), exactly as down/up/imul/iadd/ipow_list/Poly.eval do with numpy; the monomials
 * are visited in the same order for every box.  The result is therefore the same
 * enclosure, bit for bit, computed without numpy's temporary arrays.  The boxes are
 * taken in blocks so that the inner loops run over boxes and can use vector
 * instructions.  Compile without -ffast-math and with -ffp-contract=off, so that no
 * product and sum is fused.
 */
#include <stdint.h>
#include <string.h>

#define MAXDEG 64
#define BLK 512

typedef union { double d; int64_t i; } dbits;

/* np.nextafter(x, +inf) and np.nextafter(x, -inf) for finite x, without branches on
 * the sign: one step of the bit pattern away from or towards zero; +-0 go to the
 * smallest subnormal of the right sign. */
static inline double up1(double x) {
    dbits u; u.d = x;
    int64_t s = (u.i >> 63) | 1;              /* +1 for x >= +0, -1 for x <= -0 */
    dbits r; r.i = u.i + s;
    dbits z; z.i = 1;                          /* +4.9e-324 */
    return x == 0.0 ? z.d : r.d;
}
static inline double dn1(double x) {
    dbits u; u.d = x;
    int64_t s = (u.i >> 63) | 1;
    dbits r; r.i = u.i - s;
    dbits z; z.i = (int64_t)0x8000000000000001ULL;   /* -4.9e-324 */
    return x == 0.0 ? z.d : r.d;
}

/* np.minimum / np.maximum for non-NaN operands, written so that the compiler may use
 * minpd / maxpd: (b < a) ? b : a equals (a <= b) ? a : b, also for equal zeros. */
static inline double mn2(double a, double b) { return b < a ? b : a; }
static inline double mx2(double a, double b) { return b > a ? b : a; }

#define IMUL(a0, a1, b0, b1, r0, r1) do { \
    double p0_ = (a0) * (b0), p1_ = (a0) * (b1), p2_ = (a1) * (b0), p3_ = (a1) * (b1); \
    r0 = dn1(mn2(mn2(mn2(p0_, p1_), p2_), p3_)); \
    r1 = up1(mx2(mx2(mx2(p0_, p1_), p2_), p3_)); } while (0)

/* L, H: (3, n) row-major (box lower and upper ends); exps: (nm, 3); clo, chi: (nm,);
 * deg: (3,) largest exponent per variable; out lo, hi: (n,). Returns 0, or -1 if a
 * degree exceeds MAXDEG. */
int ipoly_eval(int64_t n, const double *L, const double *H, int64_t nm, const int64_t *exps,
               const double *clo, const double *chi, const int64_t *deg, double *olo, double *ohi) {
    for (int v = 0; v < 3; v++)
        if (deg[v] >= MAXDEG) return -1;
    static double pwl[3][MAXDEG][BLK], pwh[3][MAXDEG][BLK];
    double lo[BLK], hi[BLK];
    for (int64_t s = 0; s < n; s += BLK) {
        int64_t m = n - s < BLK ? n - s : BLK;
        for (int v = 0; v < 3; v++) {
            const double *a0 = L + v * n + s, *a1 = H + v * n + s;
            for (int64_t i = 0; i < m; i++) { pwl[v][0][i] = 1.0; pwh[v][0][i] = 1.0; }
            for (int k = 1; k <= deg[v]; k++) {
                double *pl = pwl[v][k - 1], *ph = pwh[v][k - 1], *ql = pwl[v][k], *qh = pwh[v][k];
                for (int64_t i = 0; i < m; i++) {
                    double r0, r1;
                    IMUL(pl[i], ph[i], a0[i], a1[i], r0, r1);
                    ql[i] = r0; qh[i] = r1;
                }
            }
        }
        for (int64_t i = 0; i < m; i++) { lo[i] = 0.0; hi[i] = 0.0; }
        for (int64_t j = 0; j < nm; j++) {
            const double *xl = pwl[0][exps[3 * j]], *xh = pwh[0][exps[3 * j]];
            const double *yl = pwl[1][exps[3 * j + 1]], *yh = pwh[1][exps[3 * j + 1]];
            const double *zl = pwl[2][exps[3 * j + 2]], *zh = pwh[2][exps[3 * j + 2]];
            const double cl = clo[j], ch = chi[j];
            for (int64_t i = 0; i < m; i++) {
                double m0, m1, t0, t1;
                IMUL(xl[i], xh[i], yl[i], yh[i], m0, m1);
                IMUL(m0, m1, zl[i], zh[i], t0, t1);
                IMUL(t0, t1, cl, ch, m0, m1);
                lo[i] = dn1(lo[i] + m0); hi[i] = up1(hi[i] + m1);
            }
        }
        for (int64_t i = 0; i < m; i++) { olo[s + i] = lo[i]; ohi[s + i] = hi[i]; }
    }
    return 0;
}
