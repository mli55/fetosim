/* Hill EP, Power GG, Longo LD. A mathematical model of placental O2 transfer with consideration of hemoglobin
   reaction rates. Am J Physiol 222(3): 721-729, 1972. doi:10.1152/ajplegacy.1972.222.3.721
   Reimplemented from the paper. Tags: H(n) = equation n, H-T1 = Table 1, H-AIV = Appendix IV.

   Concurrent exchange between a maternal and a fetal placental capillary: O2 leaves maternal blood and enters fetal
   blood at the rate Dp (Pm - Pf), with 1/Dp = 1/(theta_m Vm) + 1/DM + 1/(theta_f Vf), H(5). Hemoglobin-bound O2
   only (the paper adds dissolved O2 only above 100 mmHg, Appendix V).

   The paper integrates H(1), H(2) over the fetal transit time with RK4. Here the same equations are integrated
   along the capillary (fraction x of its length, which is the paper's time divided by the transit time; H-AIV
   keeps both bloods at the same x) with an exponential step: over one step the pressure difference decays at the
   local rate lambda = Dp (dPm/dCm / Qm + dPf/dCf / Qf). This stays stable when a flow becomes very small, which
   happens during contractions; at the paper's standard values it agrees with RK4 (tests/test_hill1972.py).
   Parameters HP[] use the order of model/parameters/placental_o2_hill1972.json. */
#include <math.h>

enum {
  H_MODE,                                   /* used by core.c: 0 = Jongen A(16) diffusion, 1 = this exchanger */
  H_K1M, H_K2M, H_K3M, H_K1F, H_K2F, H_K3F, /* H(3) dissociation constants, maternal and fetal */
  H_HBM, H_HBF, H_PHM, H_PHF,               /* hemoglobin g/100 ml, pH */
  H_VM, H_VF, H_DM,                         /* capillary blood volumes ml, membrane diffusing capacity ml/(min mmHg) */
  H_KCM, H_KCF, H_KCREF,                    /* O2 - red cell reaction rate constants, 1/(mM s), H(6) */
  H_C1, H_C2, H_C3, H_C4,                   /* H(6) */
  H_NSTEP, H_QMIN,                          /* integration steps; flows below H_QMIN ml/s give no exchange */
  H_NPAR
};

static double h_sat(double P, double k1, double k2, double k3, double pH) {                /* H(3) solved for S */
  double x = pow(10.0, (log10(P) - k1 + k2*(pH - 7.4))/k3); return 100*x/(1 + x);
}
static double h_p(double S, double k1, double k2, double k3, double pH) {                  /* H(3) */
  return pow(10.0, k1 - k2*(pH - 7.4) + k3*log10(S/(100 - S)));
}
static double h_theta(const double *H, double S, double hb, double kc) {                   /* H(6) */
  return (H[H_C1]*(1 - exp(-H[H_C2]*(100 - S))) + H[H_C3]*S - H[H_C4])*(1.34*hb/100/0.2)*(kc/H[H_KCREF]);
}

/* O2 transfer rate (ml O2/s) for inflowing maternal and fetal PO2 (mmHg) and capillary flows qm, qf (ml/s).
   end[0..3] (may be NULL) receive the end-capillary maternal PO2, fetal PO2 and the mean Dp, ml/(min mmHg). */
static double hill_exchange(const double *H, double pm_in, double pf_in, double qm, double qf, double *end) {
  if (end) { end[0] = pm_in; end[1] = pf_in; end[2] = 0; }
  if (qm < H[H_QMIN] || qf < H[H_QMIN] || pm_in <= pf_in) return 0.0;
  const double capm = 1.34*H[H_HBM]/100, capf = 1.34*H[H_HBF]/100;                        /* ml O2/ml, H(4) */
  double sm = h_sat(pm_in, H[H_K1M], H[H_K2M], H[H_K3M], H[H_PHM]);
  double sf = h_sat(pf_in, H[H_K1F], H[H_K2F], H[H_K3F], H[H_PHF]), sf0 = sf, sdp = 0;
  const int n = (int)lround(H[H_NSTEP]); const double dx = 1.0/n;
  for (int i = 0; i < n; i++) {
    double pm = h_p(sm, H[H_K1M], H[H_K2M], H[H_K3M], H[H_PHM]);
    double pf = h_p(sf, H[H_K1F], H[H_K2F], H[H_K3F], H[H_PHF]);
    double dp = 1/(1/(h_theta(H, sm, H[H_HBM], H[H_KCM])*H[H_VM]) + 1/H[H_DM]
                   + 1/(h_theta(H, sf, H[H_HBF], H[H_KCF])*H[H_VF]))/60;                   /* H(5), ml/(s mmHg) */
    sdp += dp;
    double dpdcm = pm*H[H_K3M]*100/(sm*(100 - sm))*100/capm;                               /* dP/dC from H(3), H(4) */
    double dpdcf = pf*H[H_K3F]*100/(sf*(100 - sf))*100/capf;
    double lam = dp*(dpdcm/qm + dpdcf/qf)*dx;
    double j = dp*(pm - pf)*dx*(lam > 1e-9 ? (1 - exp(-lam))/lam : 1.0);                   /* ml O2/s over this step */
    sm -= j/qm*100/capm; sf += j/qf*100/capf;                                              /* H(1), H(2), H-AIV */
    if (sm < 1e-6) sm = 1e-6;
    if (sf > 100 - 1e-6) sf = 100 - 1e-6;
  }
  if (end) { end[0] = h_p(sm, H[H_K1M], H[H_K2M], H[H_K3M], H[H_PHM]);
             end[1] = h_p(sf, H[H_K1F], H[H_K2F], H[H_K3F], H[H_PHF]); end[2] = sdp/n*60; }
  return qf*(sf - sf0)/100*capf;
}
