#include "ext.h"
#include "ext_obex.h"
#include <math.h>
#include <stdlib.h>
#include <complex.h>  // C99 complex support

// Our object structure.
typedef struct _simplemax {
    t_object ob;
    long n;      // Matrix dimension
    double a;    // Potential parameter
    void *out;   // Outlet pointer
} t_simplemax;

// Global class pointer.
static t_class *simplemax_class = NULL;

/* Function prototypes */
void *simplemax_new(t_symbol *s, long argc, t_atom *argv);
void simplemax_free(t_simplemax *x);
void simplemax_assist(t_simplemax *x, void *b, long m, long a, char *s);
void simplemax_bang(t_simplemax *x);

// Helper: Compute Fourier matrix F (n x n) with entries F[k][l] = (1/sqrt(n)) * exp(2πi*k*l/n)
void compute_fourier_matrix(double complex **F, long n) {
    double norm = 1.0 / sqrt((double)n);
    for (long k = 0; k < n; k++) {
        for (long l = 0; l < n; l++) {
            double angle = 2.0 * M_PI * k * l / n;
            F[k][l] = norm * (cos(angle) + I * sin(angle));
        }
    }
}

// Helper: Compute the Hamiltonian matrix H = 0.5*(P^2 + Q^2)
void compute_harmonic_oscillator(t_simplemax *x, double complex **H) {
    long n = x->n;
    double a = x->a;
    
    // Allocate Fourier matrices F and Finv.
    double complex **F = (double complex **)malloc(n * sizeof(double complex *));
    double complex **Finv = (double complex **)malloc(n * sizeof(double complex *));
    for (long i = 0; i < n; i++) {
        F[i] = (double complex *)malloc(n * sizeof(double complex));
        Finv[i] = (double complex *)malloc(n * sizeof(double complex));
    }
    compute_fourier_matrix(F, n);
    
    // Finv is the conjugate transpose of F.
    for (long i = 0; i < n; i++) {
        for (long j = 0; j < n; j++) {
            Finv[i][j] = conj(F[j][i]);
        }
    }
    
    // Create diagonal PImpulse = diag(0, 1, 2, ... n-1)
    double *PImpulse = (double *)malloc(n * sizeof(double));
    for (long i = 0; i < n; i++)
        PImpulse[i] = i;
    
    // Multiply: compute P_temp = (PImpulse * F) (diagonal multiplication)
    double complex **P_temp = (double complex **)malloc(n * sizeof(double complex *));
    for (long i = 0; i < n; i++) {
        P_temp[i] = (double complex *)malloc(n * sizeof(double complex));
        for (long j = 0; j < n; j++) {
            P_temp[i][j] = PImpulse[i] * F[i][j];
        }
    }
    
    // Compute P = Finv * P_temp.
    double complex **P = (double complex **)malloc(n * sizeof(double complex *));
    for (long i = 0; i < n; i++) {
        P[i] = (double complex *)malloc(n * sizeof(double complex));
        for (long j = 0; j < n; j++) {
            P[i][j] = 0;
            for (long k = 0; k < n; k++)
                P[i][j] += Finv[i][k] * P_temp[k][j];
        }
    }
    
    // Create Q diagonal: Q[i] = -((n-1)*a/2) + i
    double *Q_diag = (double *)malloc(n * sizeof(double));
    for (long i = 0; i < n; i++) {
        Q_diag[i] = -((n - 1) * a / 2.0) + i;
    }
    
    // Compute P^2 = P * P
    double complex **P2 = (double complex **)malloc(n * sizeof(double complex *));
    for (long i = 0; i < n; i++) {
        P2[i] = (double complex *)malloc(n * sizeof(double complex));
        for (long j = 0; j < n; j++) {
            P2[i][j] = 0;
            for (long k = 0; k < n; k++)
                P2[i][j] += P[i][k] * P[k][j];
        }
    }
    
    // Q^2 is diagonal: (Q_diag[i])^2
    double *Q2_diag = (double *)malloc(n * sizeof(double));
    for (long i = 0; i < n; i++)
        Q2_diag[i] = Q_diag[i] * Q_diag[i];
    
    // H = 0.5*(P^2 + Q^2)
    for (long i = 0; i < n; i++) {
        for (long j = 0; j < n; j++) {
            if (i == j) {
                H[i][j] = 0.5 * (P2[i][j] + Q2_diag[i]);
            } else {
                H[i][j] = 0.5 * P2[i][j];
            }
        }
    }
    
    // Free temporary allocations
    for (long i = 0; i < n; i++) {
        free(F[i]);
        free(Finv[i]);
        free(P_temp[i]);
        free(P[i]);
        free(P2[i]);
    }
    free(F);
    free(Finv);
    free(P_temp);
    free(P);
    free(PImpulse);
    free(Q_diag);
    free(Q2_diag);
}

// ext_main: Called when the external is loaded.
void ext_main(void *r) {
    t_class *c;
    // Change "quantum_ho" to "simplemax" so it matches the old Info.plist if desired
    c = class_new("simplemax",               // If you want an object named "simplemax"
                  (method)simplemax_new,
                  (method)simplemax_free,
                  sizeof(t_simplemax),
                  0L,
                  A_GIMME,
                  0);

    class_addmethod(c, (method)simplemax_bang, "bang", 0);
    class_addmethod(c, (method)simplemax_assist, "assist", A_CANT, 0);
    class_register(CLASS_BOX, c);
    
    post("simplemax external registered successfully");
    
    simplemax_class = c; // Save the class pointer
}

// Constructor: Create a new instance.
void *simplemax_new(t_symbol *s, long argc, t_atom *argv) {
    t_simplemax *x = (t_simplemax *)object_alloc(simplemax_class);
    if (x) {
        x->n = 8;
        x->a = 1.0;
        if (argc >= 1)
            x->n = atom_getlong(argv);
        if (argc >= 2)
            x->a = atom_getfloat(argv + 1);
        
        x->out = outlet_new(x, NULL);
    }
    return x;
}

// Destructor: Free the object.
void simplemax_free(t_simplemax *x) {
    // No special resources to free
}

// Assist: Provide messages for inlets/outlets.
void simplemax_assist(t_simplemax *x, void *b, long m, long a, char *s) {
    if (m == 1) {
        sprintf(s, "Bang to compute Hamiltonian matrix");
    } else {
        sprintf(s, "Output: Flattened Hamiltonian matrix as a list");
    }
}

// Method: Compute and output the Hamiltonian matrix when receiving a bang.
void simplemax_bang(t_simplemax *x) {
    long n = x->n;
    double complex **H = (double complex **)malloc(n * sizeof(double complex *));
    for (long i = 0; i < n; i++) {
        H[i] = (double complex *)malloc(n * sizeof(double complex));
    }
    
    compute_harmonic_oscillator(x, H);
    
    // Output real parts as a flat list
    t_atom *out_list = (t_atom *)sysmem_newptr(n * n * sizeof(t_atom));
    for (long i = 0; i < n; i++) {
        for (long j = 0; j < n; j++) {
            atom_setfloat(out_list + (i * n + j), creal(H[i][j]));
        }
    }
    outlet_list(x->out, gensym("list"), n * n, out_list);
    sysmem_freeptr(out_list);
    
    for (long i = 0; i < n; i++) {
        free(H[i]);
    }
    free(H);
}

