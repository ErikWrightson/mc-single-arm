#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "F1F209Wrapper.hh"

static const double kDEG = 3.14159265358979323846 / 180.0;

int main()
{
    const double M = 0.9383;

    double F1n, F2n, F1p, F2p, F1d, F2d, r;
    double F1be, F2be, rbe;
    double F1dqe, F2dqe;
    double Q2, W2, nu, x;

    F1F209Wrapper *model = new F1F209Wrapper();

    for (int i = 1; i <= 2; i++) {
        Q2 = 1.0 * i;
        for (int j = 1; j <= 5; j++) {
            W2 = 1. + 0.5 * j;
            nu = (W2 - M * M + Q2) / 2. / M;
            x = Q2 / 2. / M / nu;
            model->GetF1F2IN09(0., 1., Q2, W2, F1n, F2n, r);
            model->GetF1F2IN09(1., 1., Q2, W2, F1p, F2p, r);
            model->GetF1F2IN09(1., 2., Q2, W2, F1d, F2d, r);
            model->GetF1F2IN09(4., 9., Q2, W2, F1be, F2be, rbe);
            printf("%5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf\n", Q2, W2, x, F2n, F2p, F2d, F1p, F1be, r, rbe);
        }
    }

    for (int i = 1; i <= 3; i++) {
        Q2 = 1.0 * i;
        W2 = M*M;
        model->GetF1F2QE09(1., 2., Q2, W2, F1dqe, F2dqe);
        printf("%5.3lf %5.3lf %5.3lf\n", Q2, F1dqe, F2dqe);
    }
}
