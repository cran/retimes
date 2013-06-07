void skewcheck(double *x, int *n, int *check);
void histestim(double *x, int *n, double *smoothing, double *value);
void distestim(double *x, int *n, double *smoothing, double *value);
void kernestim(double *x, int *n, double *smoothing, double *value);
double smoothing(double *x, int n);
double kernel(double *x, int n, int j);
double quantile(double *x, int n, double p);
