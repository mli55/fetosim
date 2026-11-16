/* Ursino M. Interaction between carotid baroregulation and the pulsating heart: a mathematical model.
   Am J Physiol 275 (Heart Circ Physiol 44): H1733-H1747, 1998. doi:10.1152/ajpheart.1998.275.5.H1733
   Reimplemented from the paper. Tags: U(n) = appendix equation n, U-T1/T2/T3 = Tables 1-3.

   Closed-loop pulsatile circulation (left and right heart with time-varying elastance, systemic arteries, splanchnic
   and extrasplanchnic peripheral and venous beds, pulmonary arteries, peripheral bed and veins) with the carotid
   baroreflex acting on peripheral resistances, venous unstressed volumes, end-systolic elastances and heart period.

   Integration: fixed-step RK4. Pure delays use the step-start history of f_es and f_ev (held over one step).
   Parameters UP[] and state US[] use the order of model/parameters/maternal_circulation_ursino1998.json and of
   USTATE in model/ursino1998.py. */
#include <math.h>
#include <stdlib.h>
#include <string.h>

enum {  /* parameters */
  /* U-T1 vascular system */
  U_CSA, U_CSP, U_CEP, U_CSV, U_CEV, U_CPA, U_CPP, U_CPV,
  U_VUSA, U_VUSP, U_VUEP, U_VUPA, U_VUPP, U_VUPV,
  U_RSA, U_RSV, U_REV, U_RPA, U_RPP, U_RPV, U_LSA, U_LPA, U_VT,
  /* U-T2 heart */
  U_CLA, U_VULA, U_RLA, U_P0LV, U_KELV, U_VULV, U_KRLV,
  U_CRA, U_VURA, U_RRA, U_P0RV, U_KERV, U_VURV, U_KRRV, U_TSYS0, U_KSYS,
  /* U-T3 afferent and efferent pathways */
  U_PN, U_FMIN, U_FMAX, U_KA, U_TAUZ, U_TAUP,
  U_FESINF, U_FES0, U_KES, U_FESMIN, U_FEV0, U_FEVINF, U_KEV, U_FCS0,
  /* U-T3 effectors: gain, time constant, delay, value without sympathetic drive */
  U_GELV, U_TAUELV, U_DELV, U_ELV0,
  U_GERV, U_TAUERV, U_DERV, U_ERV0,
  U_GRSP, U_TAURSP, U_DRSP, U_RSP0,
  U_GREP, U_TAUREP, U_DREP, U_REP0,
  U_GVSV, U_TAUVSV, U_DVSV, U_VUSV0,
  U_GVEV, U_TAUVEV, U_DVEV, U_VUEV0,
  U_GTS, U_TAUTS, U_DTS, U_GTV, U_TAUTV, U_DTV, U_T0,
  /* experiment switches (not in the paper's tables; all 0 = the published closed-loop model) */
  U_PCS_MODE,   /* 0: P_cs = P_sa (closed loop); 1: P_cs from the input pcs (isolated carotid sinus) */
  U_PACE_T,     /* > 0: heart period fixed at this value, s (artificial pacing, Fig. 6) */
  U_VAGOTOMY,   /* 1: f_ev = 0 (Fig. 7 dashed lines) */
  U_HOLD_E, U_HOLD_R, U_HOLD_V, U_HOLD_T,   /* 1: effector group frozen at its current value (Fig. 9 bottom) */
  U_NPAR
};

enum {  /* state */
  X_PPA, X_FPA, X_PPP, X_PPV, X_PSA, X_FSA, X_PSP, X_PEV, X_PLA, X_VLV, X_PRA, X_VRV, X_XI, X_PT,
  X_DELV, X_DERV, X_DRSP, X_DREP, X_DVSV, X_DVEV, X_DTS, X_DTV, X_VT,
  U_NSTATE
};

#define U_NOUT 24

typedef struct {       /* algebraic quantities of one evaluation */
  double psv, plv, prv, pmaxlv, pmaxrv, fol, for_, fil, fir, T, phi, fcs, fes, fev;
  double rsp, rep, elv, erv, vusv, vuev, fsp, fep, fsv, fev_flow, dpsa;
} ualg;

static double u_fcs(const double *P, double pt) {                                        /* U(31) */
  double e = exp((pt - P[U_PN])/P[U_KA]); return (P[U_FMIN] + P[U_FMAX]*e)/(1 + e);
}
static double u_fes(const double *P, double fcs) {                                       /* U(33) */
  return P[U_FESINF] + (P[U_FES0] - P[U_FESINF])*exp(-P[U_KES]*fcs);
}
static double u_fev(const double *P, double fcs) {                                       /* U(34) */
  if (P[U_VAGOTOMY] > 0.5) return 0.0;
  double e = exp((fcs - P[U_FCS0])/P[U_KEV]); return (P[U_FEV0] + P[U_FEVINF]*e)/(1 + e);
}
static double u_sig(const double *P, double g, double fes_d) {                           /* U(35), U(38) */
  return fes_d >= P[U_FESMIN] ? g*log(fes_d - P[U_FESMIN] + 1) : 0.0;
}

/* dS/dt and algebraic outputs. fd[] = delayed firing rates: f_es(t - D) for Elv, Erv, Rsp, Rep, Vusv, Vuev, Ts and
   f_ev(t - D_Tv). pcs, dpcs: external carotid sinus pressure and its rate (used when U_PCS_MODE = 1). inflow = I(t),
   ml/s into the total blood volume, U(11). q_xa = extra flow leaving the systemic arterial node P_sa and q_xv = extra
   flow entering the extrasplanchnic veins P_ev (both 0 in the paper; the uterine branch when coupled to the fetal
   model, where inflow = q_xv - q_xa keeps U(9) consistent). */
static void u_deriv(const double *P, const double *S, const double *fd, double pcs, double dpcs, double inflow,
                    double q_xa, double q_xv, double *D, ualg *a) {
  a->elv = P[U_ELV0] + S[X_DELV]; a->erv = P[U_ERV0] + S[X_DERV];                        /* U(37) */
  a->rsp = P[U_RSP0] + S[X_DRSP]; a->rep = P[U_REP0] + S[X_DREP];
  a->vusv = P[U_VUSV0] + S[X_DVSV]; a->vuev = P[U_VUEV0] + S[X_DVEV];
  a->T = P[U_PACE_T] > 0 ? P[U_PACE_T] : S[X_DTS] + S[X_DTV] + P[U_T0];                  /* U(42) */
  double tsys = P[U_TSYS0] - P[U_KSYS]/a->T;                                             /* U(22) */
  double u = S[X_XI] - floor(S[X_XI]);                                                   /* U(21) */
  double s = sin(M_PI*a->T/tsys*u);
  a->phi = u <= tsys/a->T ? s*s : 0.0;                                                   /* U(19) */

  double vu = P[U_VUSA] + P[U_VUSP] + P[U_VUEP] + a->vusv + a->vuev + P[U_VURA]
            + P[U_VUPA] + P[U_VUPP] + P[U_VUPV] + P[U_VULA];                             /* U(10) */
  a->psv = (S[X_VT] - P[U_CSA]*S[X_PSA] - (P[U_CSP] + P[U_CEP])*S[X_PSP] - P[U_CEV]*S[X_PEV] - P[U_CRA]*S[X_PRA]
            - S[X_VRV] - P[U_CPA]*S[X_PPA] - P[U_CPP]*S[X_PPP] - P[U_CPV]*S[X_PPV] - P[U_CLA]*S[X_PLA] - S[X_VLV]
            - vu)/P[U_CSV];                                                              /* U(9) */

  /* left heart, U(12)-U(18) */
  a->pmaxlv = a->phi*a->elv*(S[X_VLV] - P[U_VULV]) + (1 - a->phi)*P[U_P0LV]*(exp(P[U_KELV]*S[X_VLV]) - 1);
  double rlv = P[U_KRLV]*a->pmaxlv;
  a->fol = a->pmaxlv > S[X_PSA] ? (a->pmaxlv - S[X_PSA])/rlv : 0.0;
  a->plv = a->pmaxlv - rlv*a->fol;
  a->fil = S[X_PLA] > a->plv ? (S[X_PLA] - a->plv)/P[U_RLA] : 0.0;
  /* right heart, U(23)-U(29) */
  a->pmaxrv = a->phi*a->erv*(S[X_VRV] - P[U_VURV]) + (1 - a->phi)*P[U_P0RV]*(exp(P[U_KERV]*S[X_VRV]) - 1);
  double rrv = P[U_KRRV]*a->pmaxrv;
  a->for_ = a->pmaxrv > S[X_PPA] ? (a->pmaxrv - S[X_PPA])/rrv : 0.0;
  a->prv = a->pmaxrv - rrv*a->for_;
  a->fir = S[X_PRA] > a->prv ? (S[X_PRA] - a->prv)/P[U_RRA] : 0.0;

  a->fsp = (S[X_PSP] - a->psv)/a->rsp; a->fep = (S[X_PSP] - S[X_PEV])/a->rep;
  a->fsv = (a->psv - S[X_PRA])/P[U_RSV]; a->fev_flow = (S[X_PEV] - S[X_PRA])/P[U_REV];
  double fpp = (S[X_PPP] - S[X_PPV])/P[U_RPP], fpv = (S[X_PPV] - S[X_PLA])/P[U_RPV];

  /* effectors, U(35)-U(41) */
  double hE = P[U_HOLD_E] > 0.5, hR = P[U_HOLD_R] > 0.5, hV = P[U_HOLD_V] > 0.5, hT = P[U_HOLD_T] > 0.5;
  D[X_DELV] = hE ? 0 : (u_sig(P, P[U_GELV], fd[0]) - S[X_DELV])/P[U_TAUELV];
  D[X_DERV] = hE ? 0 : (u_sig(P, P[U_GERV], fd[1]) - S[X_DERV])/P[U_TAUERV];
  D[X_DRSP] = hR ? 0 : (u_sig(P, P[U_GRSP], fd[2]) - S[X_DRSP])/P[U_TAURSP];
  D[X_DREP] = hR ? 0 : (u_sig(P, P[U_GREP], fd[3]) - S[X_DREP])/P[U_TAUREP];
  D[X_DVSV] = hV ? 0 : (u_sig(P, P[U_GVSV], fd[4]) - S[X_DVSV])/P[U_TAUVSV];
  D[X_DVEV] = hV ? 0 : (u_sig(P, P[U_GVEV], fd[5]) - S[X_DVEV])/P[U_TAUVEV];
  D[X_DTS] = hT ? 0 : (u_sig(P, P[U_GTS], fd[6]) - S[X_DTS])/P[U_TAUTS];
  D[X_DTV] = hT ? 0 : (P[U_GTV]*fd[7] - S[X_DTV])/P[U_TAUTV];

  /* circulation, U(1)-U(8), U(11), U(12), U(14), U(23), U(25) */
  D[X_PPA] = (a->for_ - S[X_FPA])/P[U_CPA];
  D[X_FPA] = (S[X_PPA] - S[X_PPP] - P[U_RPA]*S[X_FPA])/P[U_LPA];
  D[X_PPP] = (S[X_FPA] - fpp)/P[U_CPP];
  D[X_PPV] = (fpp - fpv)/P[U_CPV];
  a->dpsa = (a->fol - S[X_FSA] - q_xa)/P[U_CSA];
  D[X_PSA] = a->dpsa;
  D[X_FSA] = (S[X_PSA] - S[X_PSP] - P[U_RSA]*S[X_FSA])/P[U_LSA];
  D[X_PSP] = (S[X_FSA] - a->fsp - a->fep)/(P[U_CSP] + P[U_CEP]);
  D[X_PEV] = (a->fep - a->fev_flow + q_xv - D[X_DVEV])/P[U_CEV];
  D[X_PLA] = (fpv - a->fil)/P[U_CLA];
  D[X_VLV] = a->fil - a->fol;
  D[X_PRA] = (a->fsv + a->fev_flow - a->fir)/P[U_CRA];
  D[X_VRV] = a->fir - a->for_;
  D[X_XI] = 1/a->T;
  D[X_VT] = inflow;

  /* afferent pathway, U(30) */
  double pc = P[U_PCS_MODE] > 0.5 ? pcs : S[X_PSA], dpc = P[U_PCS_MODE] > 0.5 ? dpcs : a->dpsa;
  D[X_PT] = (pc + P[U_TAUZ]*dpc - S[X_PT])/P[U_TAUP];
  a->fcs = u_fcs(P, S[X_PT]); a->fes = u_fes(P, a->fcs); a->fev = u_fev(P, a->fcs);
}

typedef struct { int n, head; double *fes, *fev; } uhist;   /* ring buffers of step-start f_es, f_ev */

static void u_hist_init(uhist *h, const double *P, const double *S, double dt) {
  double dmax = P[U_DVSV] > P[U_DVEV] ? P[U_DVSV] : P[U_DVEV];
  const int ds[] = {U_DELV, U_DERV, U_DRSP, U_DREP, U_DTS, U_DTV};
  for (int i = 0; i < 6; i++) if (P[ds[i]] > dmax) dmax = P[ds[i]];
  h->n = (int)lround(dmax/dt) + 2; h->head = 0;
  h->fes = malloc(h->n*sizeof(double)); h->fev = malloc(h->n*sizeof(double));
  double fcs = u_fcs(P, S[X_PT]), fes = u_fes(P, fcs), fev = u_fev(P, fcs);
  for (int i = 0; i < h->n; i++) { h->fes[i] = fes; h->fev[i] = fev; }      /* history = start value */
}
static void u_hist_free(uhist *h) { free(h->fes); free(h->fev); }
static double u_lag(const double *buf, const uhist *h, double D, double dt) {
  int k = (int)lround(D/dt); return buf[((h->head - k) % h->n + h->n) % h->n];
}

/* One RK4 step of size dt. Returns the algebraic quantities at the step start in *a0. */
static void u_step(const double *P, double *S, uhist *h, double dt, double pcs0, double pcs1, double inflow,
                   double q_xa, double q_xv, ualg *a0) {
  double k1[U_NSTATE], k2[U_NSTATE], k3[U_NSTATE], k4[U_NSTATE], Y[U_NSTATE], fd[8]; ualg a;
  double fcs = u_fcs(P, S[X_PT]);
  h->fes[h->head] = u_fes(P, fcs); h->fev[h->head] = u_fev(P, fcs);
  const int ds[] = {U_DELV, U_DERV, U_DRSP, U_DREP, U_DVSV, U_DVEV, U_DTS};
  for (int i = 0; i < 7; i++) fd[i] = u_lag(h->fes, h, P[ds[i]], dt);
  fd[7] = u_lag(h->fev, h, P[U_DTV], dt);
  h->head = (h->head + 1) % h->n;
  double dpcs = (pcs1 - pcs0)/dt, pm = 0.5*(pcs0 + pcs1);
  u_deriv(P, S, fd, pcs0, dpcs, inflow, q_xa, q_xv, k1, a0);
  for (int i = 0; i < U_NSTATE; i++) Y[i] = S[i] + 0.5*dt*k1[i];
  u_deriv(P, Y, fd, pm, dpcs, inflow, q_xa, q_xv, k2, &a);
  for (int i = 0; i < U_NSTATE; i++) Y[i] = S[i] + 0.5*dt*k2[i];
  u_deriv(P, Y, fd, pm, dpcs, inflow, q_xa, q_xv, k3, &a);
  for (int i = 0; i < U_NSTATE; i++) Y[i] = S[i] + dt*k3[i];
  u_deriv(P, Y, fd, pcs1, dpcs, inflow, q_xa, q_xv, k4, &a);
  for (int i = 0; i < U_NSTATE; i++) S[i] += dt*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6;
}

#ifndef URSINO_NO_MAIN
static double u_interp(const double *a, int n, double t, double du) {
  if (!a) return 0.0;
  double x = t/du; int i = (int)floor(x);
  if (i < 0) return a[0];
  if (i >= n - 1) return a[n - 1];
  double w = x - i; return a[i]*(1 - w) + a[i + 1]*w;
}

/* Standalone run: n_steps of size dt. pcs (mmHg) and inflow (ml/s) are sampled every du s (n_in samples, linear
   interpolation; NULL = 0). Records U_NOUT outputs every out_every steps. Returns 0, or 1 if non-finite. */
int ursino_simulate(const double *P, double *S, int n_steps, double dt, const double *pcs, const double *inflow,
                    int n_in, double du, int out_every, double *out) {
  uhist h; u_hist_init(&h, P, S, dt); ualg a; int rec = 0, bad = 0;
  for (int k = 0; k < n_steps; k++) {
    double t = k*dt, S0[U_NSTATE]; memcpy(S0, S, sizeof(S0));
    u_step(P, S, &h, dt, u_interp(pcs, n_in, t, du), u_interp(pcs, n_in, t + dt, du), u_interp(inflow, n_in, t, du),
           0.0, 0.0, &a);
    if (out_every > 0 && k % out_every == 0) {
      double *o = out + (size_t)rec*U_NOUT;
      o[0] = S0[X_PSA]; o[1] = a.plv; o[2] = S0[X_VLV]; o[3] = a.fol; o[4] = S0[X_PRA]; o[5] = a.T; o[6] = a.rsp;
      o[7] = a.rep; o[8] = a.elv; o[9] = a.erv; o[10] = a.vusv; o[11] = a.vuev; o[12] = a.fcs; o[13] = a.fes;
      o[14] = a.fev; o[15] = a.psv; o[16] = S0[X_PEV]; o[17] = S0[X_PSP]; o[18] = a.fsp; o[19] = a.fep;
      o[20] = S0[X_PPA]; o[21] = S0[X_PLA]; o[22] = S0[X_VRV]; o[23] = S0[X_XI];
      rec++;
    }
    if (!isfinite(S[X_PSA]) || !isfinite(S[X_VLV])) { bad = 1; break; }
  }
  u_hist_free(&h);
  return bad;
}

int ursino_n_par(void) { return U_NPAR; }
int ursino_n_state(void) { return U_NSTATE; }
int ursino_n_out(void) { return U_NOUT; }
#endif
