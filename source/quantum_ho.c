/* qte.quantumho.c – Quantum Harmonic Oscillator external for Max/MSP
 *
 * This external computes the Hamiltonian matrix H = 0.5*(P^2 + Q^2)
 * using:
 *    - PImpulse = diag(0, 1, ..., n-1)
 *    - F is the Fourier matrix: F[k][l] = (1/sqrt(n)) * exp(2πi*k*l/n)
 *    - Finv is the conjugate transpose of F
 *    - P = F * diag(PImpulse) * Finv  (note the new order)
 *    - Q = diag( a * ( -((n - 1)/2) + i ) )  (multiply entire expression by a)
 *    - Finally, H = 0.5 * (P^2 + Q^2)
 *
 * The result is output as a flat list of 2*n*n numbers:
 * (Re[M[0][0]], Im[M[0][0]], Re[M[0][1]], Im[M[0][1]], ...)
 */

#include "ext.h"
#include "ext_obex.h"
#include <math.h>
#include <stdlib.h>
#include <complex.h>

/* ------------------------------------------------------------
   Our object structure
   ------------------------------------------------------------ */
typedef struct _qte_quantumho {
    t_object ob;
    long n;         // Matrix dimension
    double a;       // Potential parameter
    void *out;      // Outlet pointer
} t_qte_quantumho;

/* Global class pointer */
static t_class *qte_quantumho_class = NULL;

/* ------------------------------------------------------------
   Function prototypes
   ------------------------------------------------------------ */
void *qte_quantumho_new(t_symbol *s, long argc, t_atom *argv);
void qte_quantumho_free(t_qte_quantumho *x);
void qte_quantumho_assist(t_qte_quantumho *x, void *b, long m, long a, char *s);
void qte_quantumho_bang(t_qte_quantumho *x);

/* Helper: Allocate an n x n matrix of double complex numbers. */
static double complex **alloc_complex_matrix(long n) {
    double complex **matrix = (double complex **)malloc(n * sizeof(double complex *));
    if (!matrix)
        return NULL;
    for (long i = 0; i < n; i++) {
        matrix[i] = (double complex *)calloc(n, sizeof(double complex));
        if (!matrix[i]) {
            for (long j = 0; j < i; j++)
                free(matrix[j]);
            free(matrix);
            return NULL;
        }
    }
    return matrix;
}

/* Helper: Free an n x n complex matrix. */
static void free_complex_matrix(double complex **matrix, long n) {
    if (matrix) {
        for (long i = 0; i < n; i++)
            free(matrix[i]);
        free(matrix);
    }
}

/* Helper: Compute the Fourier matrix F (n x n) */
static void compute_fourier_matrix(double complex **F, long n) {
    double norm = 1.0 / sqrt((double)n);
    for (long k = 0; k < n; k++) {
        for (long l = 0; l < n; l++) {
            double angle = 2.0 * M_PI * k * l / n;
            F[k][l] = norm * (cos(angle) + I * sin(angle));
        }
    }
}

/* Helper: Conjugate transpose Finv[i][j] = conj(F[j][i]). */
static void compute_conjugate_transpose(double complex **F, double complex **Finv, long n) {
    for (long i = 0; i < n; i++) {
        for (long j = 0; j < n; j++) {
            Finv[i][j] = conj(F[j][i]);
        }
    }
}

/* Helper: Multiply diag(D) with a matrix M. diag(D) * M means
   result[i][j] = D[i] * M[i][j]. */
static void multiply_diag_matrix(double *D, double complex **M, double complex **result, long n) {
    for (long i = 0; i < n; i++) {
        for (long j = 0; j < n; j++) {
            result[i][j] = D[i] * M[i][j];
        }
    }
}

/* Helper: Multiply two n x n complex matrices A and B, store in result. */
static void multiply_complex_matrices(double complex **A, double complex **B, double complex **result, long n) {
    for (long i = 0; i < n; i++) {
        for (long j = 0; j < n; j++) {
            double complex sum = 0.0 + 0.0 * I;
            for (long k = 0; k < n; k++) {
                sum += A[i][k] * B[k][j];
            }
            result[i][j] = sum;
        }
    }
}

/* Helper: Round a double to 5 decimal places. */
static double round5(double x) {
    return round(x * 100000.0) / 100000.0;
}

/* Compute the Hamiltonian matrix H = 0.5 * (P^2 + Q^2),
 * with new definitions:
 *   P = F * diag(PImpulse) * Finv,
 *   Q[i] = a * ( -((n - 1)/2) + i )
 */
static double complex **compute_hamiltonian(t_qte_quantumho *x) {
    long n = x->n;
    double a = x->a;
    
    // Allocate the final matrix H and temp matrices
    double complex **H = alloc_complex_matrix(n);
    double complex **F = alloc_complex_matrix(n);
    double complex **Finv = alloc_complex_matrix(n);
    double complex **P_temp = alloc_complex_matrix(n);
    double complex **P = alloc_complex_matrix(n);
    double complex **P2 = alloc_complex_matrix(n);
    if (!H || !F || !Finv || !P_temp || !P || !P2) {
        // Clean up if any allocation fails
        if (H)    free_complex_matrix(H, n);
        if (F)    free_complex_matrix(F, n);
        if (Finv) free_complex_matrix(Finv, n);
        if (P_temp) free_complex_matrix(P_temp, n);
        if (P)    free_complex_matrix(P, n);
        if (P2)   free_complex_matrix(P2, n);
        return NULL;
    }
    
    // Build the Fourier transform and its conjugate transpose
    compute_fourier_matrix(F, n);
    compute_conjugate_transpose(F, Finv, n);
    
    // Create the diagonal vector for the momentum impulse: [0, 1, ..., n-1]
    double *PImpulse = (double *)malloc(n * sizeof(double));
    if (!PImpulse) {
        free_complex_matrix(H, n);
        free_complex_matrix(F, n);
        free_complex_matrix(Finv, n);
        free_complex_matrix(P_temp, n);
        free_complex_matrix(P, n);
        free_complex_matrix(P2, n);
        return NULL;
    }
    for (long i = 0; i < n; i++) {
        PImpulse[i] = (double)i;
    }
    
    // P_temp = diag(PImpulse) * Finv
    // This yields diag(...) * Finv in P_temp
    multiply_diag_matrix(PImpulse, Finv, P_temp, n);
    
    // P = F * P_temp = F * (diag(...) * Finv)
    // => F diag(...) Finv
    multiply_complex_matrices(F, P_temp, P, n);
    
    // P^2 = P * P
    multiply_complex_matrices(P, P, P2, n);
    
    // Q diagonal: Q[i] = a * ( -((n - 1)/2) + i )
    double *Q_diag = (double *)malloc(n * sizeof(double));
    for (long i = 0; i < n; i++) {
        Q_diag[i] = a * ( -((n - 1) / 2.0) + i );
    }
    
    // Finally, build H = 0.5*( P^2 + Q^2 ) on the diagonal.
    // Off-diagonal = 0.5 * P^2[i][j].
    // Diagonal gets an added Q^2.
    for (long i = 0; i < n; i++) {
        for (long j = 0; j < n; j++) {
            double complex qterm = 0.0;
            if (i == j) {
                double q2 = Q_diag[i] * Q_diag[i];
                qterm = q2 + 0.0 * I;
            }
            double complex val = 0.5 * (P2[i][j] + qterm);
            
            // Round real & imaginary parts
            double rpart = round5(creal(val));
            double ipart = round5(cimag(val));
            H[i][j] = rpart + I * ipart;
        }
    }
    
    // Cleanup
    free(PImpulse);
    free(Q_diag);
    free_complex_matrix(F, n);
    free_complex_matrix(Finv, n);
    free_complex_matrix(P_temp, n);
    free_complex_matrix(P, n);
    free_complex_matrix(P2, n);
    
    return H;
}

/* -------------------------------------------------------------------
   Max External Methods
   ------------------------------------------------------------------- */

/* ext_main: Called by Max at load time */
void ext_main(void *r) {
    t_class *c;
    // Name the Max object "qte.quantumho"
    c = class_new("qte.quantumho",
                  (method)qte_quantumho_new,
                  (method)qte_quantumho_free,
                  sizeof(t_qte_quantumho),
                  0L,
                  A_GIMME,
                  0);
    
    class_addmethod(c, (method)qte_quantumho_bang, "bang", 0);
    class_addmethod(c, (method)qte_quantumho_assist, "assist", A_CANT, 0);
    class_register(CLASS_BOX, c);
    qte_quantumho_class = c;
}

/* Constructor */
void *qte_quantumho_new(t_symbol *s, long argc, t_atom *argv) {
    t_qte_quantumho *x = (t_qte_quantumho *)object_alloc(qte_quantumho_class);
    if (x) {
        x->n = 8;   // default dimension
        x->a = 1.0; // default potential parameter
        if (argc >= 1) {
            if (atom_gettype(argv) == A_LONG) {
                x->n = atom_getlong(argv);
            } else if (atom_gettype(argv) == A_FLOAT) {
                x->n = (long)atom_getfloat(argv);
            }
        }
        if (argc >= 2) {
            x->a = atom_getfloat(argv + 1);
        }
        x->out = outlet_new(x, NULL);
    }
    return x;
}

/* Destructor */
void qte_quantumho_free(t_qte_quantumho *x) {
    // Nothing special to free
}

/* Assist: Provide inlet/outlet assistance */
void qte_quantumho_assist(t_qte_quantumho *x, void *b, long m, long a, char *s) {
    if (m == 1) // inlet
        sprintf(s, "Bang to compute Hamiltonian");
    else        // outlet
        sprintf(s, "Outputs real,imag pairs of H as a list");
}

/* Bang method: compute & output Hamiltonian as real/imag pairs */
void qte_quantumho_bang(t_qte_quantumho *x) {
    long n = x->n;
    double complex **H = compute_hamiltonian(x);
    if (!H) {
        object_error((t_object *)x, "Failed to compute Hamiltonian (out of memory?)");
        return;
    }
    // 2 floats (real, imag) per matrix entry
    long list_size = 2 * n * n;
    t_atom *out_list = (t_atom *)sysmem_newptr(list_size * sizeof(t_atom));
    if (!out_list) {
        object_error((t_object *)x, "Failed to allocate memory for output list");
        free_complex_matrix(H, n);
        return;
    }
    
    // Flatten real & imaginary parts
    for (long i = 0; i < n; i++) {
        for (long j = 0; j < n; j++) {
            long idx = 2 * (i * n + j);
            atom_setfloat(out_list + idx,      creal(H[i][j]));
            atom_setfloat(out_list + idx + 1,  cimag(H[i][j]));
        }
    }
    outlet_list(x->out, gensym("list"), list_size, out_list);
    
    sysmem_freeptr(out_list);
    free_complex_matrix(H, n);
}

