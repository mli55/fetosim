/* Jongen, van der Hout-van der Jagt et al. (J Biomech 2016, Part A/B; thesis Ch. 4 noise) CTG model.
   Reimplemented from the papers; equation tags A(n)/B(n) = Part A / Part B equation n.
   Forward Euler. Gap-filling choices are listed in the parameter notes (model/parameters/jongen2016.json).
   Parameters P[] and state S[] use the index order of model/parameters/jongen2016.json and STATE in model/simulate.py. */
#include <math.h>
#include <string.h>
#include <stdlib.h>

enum {  /* parameters, same order as model/parameters/jongen2016.json */
  /* fetal hemodynamics */
  F_TREF, F_VART0, F_VVEN0, F_VTIS0, F_VCER0, F_VVIL0, F_CART, F_CVEN, F_CTIS, F_CCER, F_CVIL,
  F_RTISA, F_RTISV, F_RCERA, F_RCERV, F_RUMA, F_RUMV, F_EMAX, F_EMIN, F_V0MAX, F_V0MIN,
  /* maternal hemodynamics */
  M_T, M_VART0, M_VVEN0, M_VTIS0, M_CART, M_CVEN, M_CTIS, M_RTISA, M_RTISV, M_RUTA, M_RUTV, M_R0UTV, M_P0UTV, M_P1UTV,
  M_EMAX, M_EMIN, M_V0MAX, M_V0MIN, IVS_VMAX, IVS_P0, IVS_P1,
  /* umbilical (B-T1) */
  U_R0UMV, U_P0UMV, U_P1UMV, U_R0UMA, U_P0UMA, U_P1UMA,
  /* oxygen (A-T2) */
  O_DPLAC, O_DLUNG, O_QTISF, O_QCERF, O_QTISM, O_THCER, O_THTIS, O_THTISM, O_PAIR, O_ALPHA, O_BETA, O_HBF, O_HBM,
  O_C1F, O_C1M, O_C2F, O_C2M,
  /* regulation (A-T3) */
  G_CREF, G_GAMMA, G_TAURC, G_ETA1B, G_ETA2B, G_RBMIN, G_RBMAX, G_TAURB, G_YREFB,
  G_ETA1C, G_ETA2C, G_RCMIN, G_RCMAX, G_TAURCH, G_YREFC, G_WBS, G_WCS, G_WBV, G_WCV,
  G_DTV, G_TAUTV, G_KTV, G_DTS, G_TAUTS, G_KTS, G_DE, G_TAUE, G_KE, G_DR, G_TAUR, G_KR, G_DV, G_TAUV, G_KV,
  /* scenario */
  S_PREST, S_WUM, S_GAIN, S_MRT, S_MRB, S_MEV,
  /* maternal heart-rate reflex (Randall et al. 2019; model/parameters/maternal_hr_randall2019.json, appended) */
  MR_MODE, MR_SBP_OFF, MR_PTH_REST,
  MR_A, MR_B, MR_KB, MR_KPB, MR_KPR, MR_KS, MR_TAUB, MR_TAUPB, MR_TAUPR, MR_TAUS, MR_TAUH,
  MR_QW, MR_QPB, MR_QPR, MR_QS, MR_SW, MR_SPB, MR_SPR, MR_SS, MR_HI, MR_HPB, MR_HPR, MR_HS, MR_DS,
  NPAR
};

enum {  /* state */
  VM_ART, VM_VEN, VM_TIS, VM_IVS, MM_ART, MM_VEN, MM_TIS, MM_IVS, MM_H,
  VF_ART, VF_VEN, VF_TIS, VF_CER, VF_VIL, MF_ART, MF_VEN, MF_TIS, MF_CER, MF_VIL, MF_H,
  HM_PREV, HF_PREV, RCERA, RB, RC, XTV, XTS, XE, XR, XV,
  PG0, /* 11 pO2 Newton warm starts: m art ven tis ivs h, f art ven tis cer vil h */
  R_EBC = PG0 + 11, R_EBA, R_TPB, R_TPR, R_TS, R_H,   /* maternal HR reflex states (Randall 2019) */
  NSTATE
};

#define NOUT 26

static double pO2_of(double c, double ahb, double c1, double c2, double beta, double guess) {
  double p = guess > 0.5 ? guess : 20.0;
  for (int it = 0; it < 8; it++) {
    double g = p*p*p + c2*p, gp = 3*p*p + c2;
    double f = ahb*g/(g + c1) + beta*p - c;
    double d = ahb*c1*gp/((g + c1)*(g + c1)) + beta;
    double np = p - f/d;
    if (np < 0.01) np = 0.5*p;
    if (fabs(np - p) < 1e-10) { p = np; break; }
    p = np;
  }
  return p;
}

static double cO2_of(double p, double ahb, double c1, double c2, double beta) {
  double g = p*p*p + c2*p; return ahb*g/(g + c1) + beta*p;
}

static double fatan(double ptm, double p0, double p1) { return 0.5 + atan((ptm - p0)/p1)/M_PI; }

static double rvar(double ptm, double r0, double p0, double p1, double rref) {  /* A(6) */
  double f = fatan(ptm, p0, p1); double r = r0/(f*f); return r > rref ? r : rref;
}

static double gompertz(double y, double eta1, double eta2, double rmin, double rmax, double yref) {  /* A(21) */
  return rmin + (rmax - rmin)*exp(eta1*exp(eta2*(y - yref)));
}

static double interp(const double *a, int n, double t, double du) {
  if (!a) return 0.0;
  double x = t/du; int i = (int)floor(x);
  if (i < 0) return a[0];
  if (i >= n - 1) return a[n - 1];
  double w = x - i; return a[i]*(1 - w) + a[i + 1]*w;
}

static double wall_strain(double p, double A, double qw, double sw) {  /* Randall 2019 arterial wall strain */
  double e = exp(-qw*(p - sw)); return 1 - sqrt((1 + e)/(A + e));
}

/* Simulate n_steps of size dt. Inputs sampled every du seconds (n_in samples), linear interpolation:
   put = uterine pressure [mmHg]; pum = extra umbilical pressure [mmHg] added to w_um*(p_ut-p_rest) (may be NULL);
   ruta = multiplier of R_uta (NULL = 1); noise = unit-SD noise N(t) (NULL = none);
   pth = maternal thoracic pressure [mmHg] for the HR reflex (NULL = MR_PTH_REST).
   Maternal HR reflex (Randall 2019), MR_MODE: 0 = off (Jongen 2016 exactly, maternal period M_T fixed);
   1 = open loop (MHR computed from the maternal arterial pressure, M_T still fixed);
   2 = closed loop (maternal heart period = 60 / MHR).
   Records NOUT outputs every out_every steps into out (row-major). Returns 0, or 1 if the state became non-finite. */
int simulate(const double *P, double *S, int n_steps, double dt,
             const double *put, const double *pum, const double *ruta, const double *noise, const double *pth,
             int n_in, double du, int out_every, double *out) {
  const int mr = (int)lround(P[MR_MODE]);
  int nDs = mr ? (int)lround(P[MR_DS]/dt) : 0, nbT = nDs + 1, headT = 0;
  double *bts = malloc(nbT*sizeof(double));   /* Ts history; before the run it is held at the start value */
  for (int i = 0; i < nbT; i++) bts[i] = S[R_TS];
  const double ahbf = P[O_ALPHA]*P[O_HBF], ahbm = P[O_ALPHA]*P[O_HBM];
  int nTv = (int)lround(P[G_DTV]/dt), nTs = (int)lround(P[G_DTS]/dt), nE = (int)lround(P[G_DE]/dt);
  int nR = (int)lround(P[G_DR]/dt), nV = (int)lround(P[G_DV]/dt);
  int nb = 1 + (nV > nTs ? nV : nTs); if (nb < 1 + nTv) nb = 1 + nTv; if (nb < 1 + nE) nb = 1 + nE; if (nb < 1 + nR) nb = 1 + nR;
  double *bv = calloc(nb, sizeof(double)), *bs = calloc(nb, sizeof(double));  /* delay buffers, steady-state history e=0 */
  int head = 0, rec = 0, bad = 0;
  double *pg = S + PG0;
  for (int k = 0; k < n_steps; k++) {
    double t = k*dt;
    double p_ut = interp(put, n_in, t, du);
    double p_um = P[S_WUM]*(p_ut - P[S_PREST]) + (pum ? interp(pum, n_in, t, du) : 0.0);    /* B(2) */
    double N = noise ? interp(noise, n_in, t, du) : 0.0;
    double r_uta = P[M_RUTA]*(ruta ? interp(ruta, n_in, t, du) : 1.0);

    /* ---------- mother ---------- */
    double pa_m = (S[VM_ART] - P[M_VART0])/P[M_CART];
    double pv_m = (S[VM_VEN] - P[M_VVEN0])/P[M_CVEN];
    double pt_m = (S[VM_TIS] - P[M_VTIS0])/P[M_CTIS];
    double fr = S[VM_IVS]/P[IVS_VMAX]; if (fr < 1e-6) fr = 1e-6; if (fr > 1 - 1e-6) fr = 1 - 1e-6;
    double pivs = P[IVS_P0] + P[IVS_P1]*tan(M_PI*(fr - 0.5)) + p_ut;                         /* A(3) inverse, absolute */
    double ptm_utv = 0.5*(pivs + pv_m) - p_ut;                                                 /* A(7) */
    double rutv = rvar(ptm_utv, P[M_R0UTV], P[M_P0UTV], P[M_P1UTV], P[M_RUTV]);
    double ved_m = P[M_V0MIN] + pv_m/P[M_EMIN], vee_m = P[M_V0MAX] + pa_m/P[M_EMAX];
    double t_m = (mr == 2) ? 60.0/S[R_H] : P[M_T];                                             /* maternal heart period */
    double hm = 0.5*(ved_m + vee_m), qco_m = (ved_m - vee_m)/t_m;
    double qv_m = (hm - S[HM_PREV])/dt, qin_m = qco_m + qv_m;
    double qtisa_m = (pa_m - pt_m)/P[M_RTISA], qtisv_m = (pt_m - pv_m)/P[M_RTISV];
    double quta = (pa_m - pivs)/r_uta, qutv = (pivs - pv_m)/rutv;

    /* ---------- fetus (transmural; p_ext = p_ut cancels in flows) ---------- */
    double xTv = S[XTV], xTs = S[XTS], xE = S[XE], xR = S[XR], xV = S[XV], g = P[S_GAIN];
    double zTv = g*P[G_KTV]*xTv, zTs = g*P[G_KTS]*xTs, zE = g*P[G_KE]*xE, zR = g*P[G_KR]*xR, zV = g*P[G_KV]*xV;
    double T = P[F_TREF]*(1 + zTv - zTs);                                                       /* A(25) */
    double emax = P[F_EMAX]*(1 + zE), vven0 = P[F_VVEN0]*(1 + zV);
    double rtisa = P[F_RTISA]*(1 + zR) + P[S_MRT]*N*P[F_RTISA];                                 /* A(24), Th 4.2 */
    double pa = (S[VF_ART] - P[F_VART0])/P[F_CART];
    double pv = (S[VF_VEN] - vven0)/P[F_CVEN];
    double pt = (S[VF_TIS] - P[F_VTIS0])/P[F_CTIS];
    double pc = (S[VF_CER] - P[F_VCER0])/P[F_CCER];
    double pvl = (S[VF_VIL] - P[F_VVIL0])/P[F_CVIL];
    double ruma = rvar(0.5*(pvl + pa) - p_um, P[U_R0UMA], P[U_P0UMA], P[U_P1UMA], P[F_RUMA]);   /* B(1) */
    double rumv = rvar(0.5*(pvl + pv) - p_um, P[U_R0UMV], P[U_P0UMV], P[U_P1UMV], P[F_RUMV]);
    double ved = P[F_V0MIN] + pv/P[F_EMIN], vee = P[F_V0MAX] + pa/emax;                         /* A(9) */
    double hf = 0.5*(ved + vee), qco = (ved - vee)/T;                                            /* A(10) */
    double qv = (hf - S[HF_PREV])/dt, qin = qco + qv;                                            /* A(11) */
    double qtisa = (pa - pt)/rtisa, qtisv = (pt - pv)/P[F_RTISV];
    double qcera = (pa - pc)/S[RCERA], qcerv = (pc - pv)/P[F_RCERV];
    double quma = (pa - pvl)/ruma, qumv = (pvl - pv)/rumv;

    /* ---------- oxygen: concentrations and pO2 ---------- */
    double c_ma = S[MM_ART]/S[VM_ART], c_mv = S[MM_VEN]/S[VM_VEN], c_mt = S[MM_TIS]/S[VM_TIS];
    double c_mi = S[MM_IVS]/S[VM_IVS], c_mh = S[MM_H]/S[HM_PREV];
    double c_fa = S[MF_ART]/S[VF_ART], c_fv = S[MF_VEN]/S[VF_VEN], c_ft = S[MF_TIS]/S[VF_TIS];
    double c_fc = S[MF_CER]/S[VF_CER], c_fl = S[MF_VIL]/S[VF_VIL], c_fh = S[MF_H]/S[HF_PREV];
    pg[0] = pO2_of(c_ma, ahbm, P[O_C1M], P[O_C2M], P[O_BETA], pg[0]);
    pg[1] = pO2_of(c_mv, ahbm, P[O_C1M], P[O_C2M], P[O_BETA], pg[1]);
    pg[2] = pO2_of(c_mt, ahbm, P[O_C1M], P[O_C2M], P[O_BETA], pg[2]);
    pg[3] = pO2_of(c_mi, ahbm, P[O_C1M], P[O_C2M], P[O_BETA], pg[3]);
    pg[4] = pO2_of(c_mh, ahbm, P[O_C1M], P[O_C2M], P[O_BETA], pg[4]);
    pg[5] = pO2_of(c_fa, ahbf, P[O_C1F], P[O_C2F], P[O_BETA], pg[5]);
    pg[6] = pO2_of(c_fv, ahbf, P[O_C1F], P[O_C2F], P[O_BETA], pg[6]);
    pg[7] = pO2_of(c_ft, ahbf, P[O_C1F], P[O_C2F], P[O_BETA], pg[7]);
    pg[8] = pO2_of(c_fc, ahbf, P[O_C1F], P[O_C2F], P[O_BETA], pg[8]);
    pg[9] = pO2_of(c_fl, ahbf, P[O_C1F], P[O_C2F], P[O_BETA], pg[9]);
    pg[10] = pO2_of(c_fh, ahbf, P[O_C1F], P[O_C2F], P[O_BETA], pg[10]);
    /* upwind convection: flow q from X to Y carries c_X if q>0, else c_Y */
#define UP(q, cx, cy) ((q) > 0 ? (q)*(cx) : (q)*(cy))
    double qd_pl = P[O_DPLAC]*(pg[3] - pg[9]), qd_lu = P[O_DLUNG]*(P[O_PAIR] - pg[1]);           /* A(16) */
    double qm_tm = pg[2] >= P[O_THTISM] ? P[O_QTISM] : P[O_QTISM]*pg[2]/P[O_THTISM];               /* A(17) */
    double qm_ft = pg[7] >= P[O_THTIS] ? P[O_QTISF] : P[O_QTISF]*pg[7]/P[O_THTIS];
    double qm_fc = pg[8] >= P[O_THCER] ? P[O_QCERF] : P[O_QCERF]*pg[8]/P[O_THCER];
    double fm_h_a = UP(qco_m, c_mh, c_ma), fm_v_h = UP(qin_m, c_mv, c_mh);
    double fm_a_t = UP(qtisa_m, c_ma, c_mt), fm_t_v = UP(qtisv_m, c_mt, c_mv);
    double fm_a_i = UP(quta, c_ma, c_mi), fm_i_v = UP(qutv, c_mi, c_mv);
    double ff_h_a = UP(qco, c_fh, c_fa), ff_v_h = UP(qin, c_fv, c_fh);
    double ff_a_t = UP(qtisa, c_fa, c_ft), ff_t_v = UP(qtisv, c_ft, c_fv);
    double ff_a_c = UP(qcera, c_fa, c_fc), ff_c_v = UP(qcerv, c_fc, c_fv);
    double ff_a_l = UP(quma, c_fa, c_fl), ff_l_v = UP(qumv, c_fl, c_fv);

    /* ---------- regulation ---------- */
    double ca = cO2_of(pg[5], ahbf, P[O_C1F], P[O_C2F], P[O_BETA]);
    double rstar = P[F_RCERA]/((1 - P[G_GAMMA]) + P[G_GAMMA]*P[G_CREF]/ca);                       /* A(19) */
    if (rstar < 0.5*P[F_RCERA]) rstar = 0.5*P[F_RCERA];
    double rbs = gompertz(pa, P[G_ETA1B], P[G_ETA2B], P[G_RBMIN], P[G_RBMAX], P[G_YREFB]);
    double rcs = gompertz(pg[5], P[G_ETA1C], P[G_ETA2C], P[G_RCMIN], P[G_RCMAX], P[G_YREFC]);
    double rb = P[G_TAURB] > 0 ? S[RB] : rbs;                                                      /* Th Ch.4: no baro LPF */
    double rbN = rb + P[S_MRB]*N;
    double ev = P[G_WBV]*rbN + P[G_WCV]*S[RC] + P[S_MEV]*N;                                       /* A(22) */
    double es = -P[G_WBS]*rbN + P[G_WCS]*S[RC];
    bv[head] = ev; bs[head] = es;
#define DEL(buf, n) (buf[(head - (n) + nb) % nb])
    double eTv = DEL(bv, nTv), eTs = DEL(bs, nTs), eE = DEL(bs, nE), eR = DEL(bs, nR), eV = DEL(bs, nV);
    head = (head + 1) % nb;

    /* ---------- maternal HR reflex: Randall 2019 (driver_baroreflex.f FCN) ---------- */
    double sbp = pa_m + P[MR_SBP_OFF], pthv = P[MR_PTH_REST], d_r[6] = {0, 0, 0, 0, 0, 0};
    if (mr) {
      if (pth) pthv = interp(pth, n_in, t, du);
      double ewc = wall_strain(sbp, P[MR_A], P[MR_QW], P[MR_SW]);              /* carotid: SBP */
      double ewa = wall_strain(sbp - pthv, P[MR_A], P[MR_QW], P[MR_SW]);       /* aortic: SBP - Pth */
      double nb_r = P[MR_B]*(ewc - S[R_EBC]) + (1 - P[MR_B])*(ewa - S[R_EBA]);
      double gpb = 1/(1 + exp(-P[MR_QPB]*(nb_r - P[MR_SPB]))), gs = 1/(1 + exp(P[MR_QS]*(nb_r - P[MR_SS])));
      double gpr = 1/(1 + exp(P[MR_QPR]*(pthv - P[MR_SPR])));
      double htil = P[MR_HI]*(1 - P[MR_HPB]*S[R_TPB] + P[MR_HPR]*S[R_TPR] + P[MR_HS]*S[R_TS]);
      bts[headT] = S[R_TS];
      double ts_lag = bts[(headT - nDs + nbT) % nbT];                             /* Ts(t - Ds) */
      headT = (headT + 1) % nbT;
      d_r[0] = (-S[R_EBC] + P[MR_KB]*ewc)/P[MR_TAUB];
      d_r[1] = (-S[R_EBA] + P[MR_KB]*ewa)/P[MR_TAUB];
      d_r[2] = (-S[R_TPB] + P[MR_KPB]*gpb)/P[MR_TAUPB];
      d_r[3] = (-S[R_TPR] + P[MR_KPR]*gpr)/P[MR_TAUPR];
      d_r[4] = (-ts_lag + P[MR_KS]*gs)/P[MR_TAUS];
      d_r[5] = (-S[R_H] + htil)/P[MR_TAUH];
    }

    if (out_every > 0 && k % out_every == 0) {
      double *o = out + (size_t)rec*NOUT;
      o[0] = 60.0/T; o[1] = pa; o[2] = pg[5]; o[3] = quta; o[4] = qutv; o[5] = S[VM_IVS]; o[6] = pg[7]; o[7] = pg[8];
      o[8] = qtisa; o[9] = qcera; o[10] = quma; o[11] = qumv; o[12] = S[VF_VIL]; o[13] = rb; o[14] = S[RC];
      o[15] = p_ut; o[16] = es; o[17] = ev; o[18] = zTv; o[19] = zTs; o[20] = qco;
      o[21] = mr ? S[R_H] : 60.0/P[M_T]; o[22] = pa_m; o[23] = qco_m; o[24] = sbp; o[25] = pthv;
      rec++;
    }

    /* ---------- Euler update ---------- */
    S[VM_ART] += dt*(qco_m - qtisa_m - quta);
    S[VM_TIS] += dt*(qtisa_m - qtisv_m);
    S[VM_IVS] += dt*(quta - qutv);
    S[VM_VEN] += dt*(qtisv_m + qutv - qin_m);
    S[MM_ART] += dt*(fm_h_a - fm_a_t - fm_a_i);
    S[MM_TIS] += dt*(fm_a_t - fm_t_v - qm_tm);
    S[MM_IVS] += dt*(fm_a_i - fm_i_v - qd_pl);
    S[MM_VEN] += dt*(fm_t_v + fm_i_v - fm_v_h + qd_lu);
    S[MM_H] += dt*(fm_v_h - fm_h_a);
    S[VF_ART] += dt*(qco - qtisa - qcera - quma);
    S[VF_TIS] += dt*(qtisa - qtisv);
    S[VF_CER] += dt*(qcera - qcerv);
    S[VF_VIL] += dt*(quma - qumv);
    S[VF_VEN] += dt*(qtisv + qcerv + qumv - qin);
    S[MF_ART] += dt*(ff_h_a - ff_a_t - ff_a_c - ff_a_l);
    S[MF_TIS] += dt*(ff_a_t - ff_t_v - qm_ft);
    S[MF_CER] += dt*(ff_a_c - ff_c_v - qm_fc);
    S[MF_VIL] += dt*(ff_a_l - ff_l_v + qd_pl);
    S[MF_VEN] += dt*(ff_t_v + ff_c_v + ff_l_v - ff_v_h);
    S[MF_H] += dt*(ff_v_h - ff_h_a);
    S[HM_PREV] = hm; S[HF_PREV] = hf;
    S[RCERA] += dt*(rstar - S[RCERA])/P[G_TAURC];                                                   /* A(18) */
    if (P[G_TAURB] > 0) S[RB] += dt*(rbs - S[RB])/P[G_TAURB];                                       /* A(20) */
    S[RC] += dt*(rcs - S[RC])/P[G_TAURCH];
    S[XTV] += dt*(eTv - xTv)/P[G_TAUTV];                                                            /* A(23) */
    S[XTS] += dt*(eTs - xTs)/P[G_TAUTS];
    S[XE] += dt*(eE - xE)/P[G_TAUE];
    S[XR] += dt*(eR - xR)/P[G_TAUR];
    S[XV] += dt*(eV - xV)/P[G_TAUV];
    if (mr) for (int i = 0; i < 6; i++) S[R_EBC + i] += dt*d_r[i];
    if (!isfinite(S[VF_ART]) || !isfinite(T) || T <= 0 || S[VF_VIL] <= 0 || S[VM_IVS] <= 0) { bad = 1; break; }
    if (mr && (!isfinite(S[R_H]) || S[R_H] <= 0)) { bad = 1; break; }
  }
  free(bv); free(bs); free(bts);
  return bad;
}

int n_par(void) { return NPAR; }
int n_state(void) { return NSTATE; }
int n_out(void) { return NOUT; }
