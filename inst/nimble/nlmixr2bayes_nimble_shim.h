#ifndef NLMIXR2BAYES_NIMBLE_SHIM_H
#define NLMIXR2BAYES_NIMBLE_SHIM_H

#ifdef __cplusplus
extern "C" {
#endif

int nlmixr2bayes_nimble_cond_batch_theta(const double *theta, int ntheta,
                                          const double *eta, int nid,
                                          int neta, double *value,
                                          double *gradEta,
                                          double *gradTheta);

#ifdef __cplusplus
}
#endif

#endif
