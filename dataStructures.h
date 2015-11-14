#ifndef __DATASTRUCTURES_H__
#define __DATASTRUCTURES_H__

typedef struct TrackInfo{
  double px, py, pz;   // Momentum
  double r, z;         // Impact parameters
  int rz_svd, rphi_svd;// SVD hits
  double atckpi;       // Identification
} TrkInfo;

typedef struct BEvent{
  int exp,run,evtn;
  int charge;
  double de,mbc; //
  double costhBcms;

  double md0,md0_raw;
  double mds,mds_raw;
  double mphi,mphi_raw;

  double pcm_d0, pcm_ds;

  double cos_hel_d0, cos_hel_ds, cos_hel_phi;

// Final state particles info
  TrkInfo k_d0;
  TrkInfo pi_d0;
  TrkInfo pi_ds;
  TrkInfo kp_phi;
  TrkInfo km_phi;

  double chi2_fit_d0,chi2_fit_ds,chi2_fit_phi;

  double ipx,ipy,ipz; //
  double boostx,boosty,boostz; //

/////////////////////////////
//  Continuum Suppression  //
/////////////////////////////
  double cos_thr;
  double thr_sig;
  double thr_oth;

  double k0mm2;
  double k0et;
  double k0hso00;
  double k0hso01;
  double k0hso02;
  double k0hso03;
  double k0hso04;
  double k0hso10;
  double k0hso12;
  double k0hso14;
  double k0hso20;
  double k0hso22;
  double k0hso24;
  double k0hoo0;
  double k0hoo1;
  double k0hoo2;
  double k0hoo3;
  double k0hoo4;

  double k1mm2;
  double k1et;
  double k1hso00;
  double k1hso01;
  double k1hso02;
  double k1hso03;
  double k1hso04;
  double k1hso10;
  double k1hso12;
  double k1hso14;
  double k1hso20;
  double k1hso22;
  double k1hso24;
  double k1hoo0;
  double k1hoo1;
  double k1hoo2;
  double k1hoo3;
  double k1hoo4;
}bEvnt;

typedef struct MCGenEvent{
  int d0_chain[9];
  int ds_chain[9];
  int phi_chain[9];
  int pi_d0_chain[9];
  int pi_ds_chain[9];
  int k_d0_chain[9];
  int kp_phi_chain[9];
  int km_phi_chain[9];
  int bpid,bpf;
  int d0id,d0f;
  int dsid,dsf;
  int phiid,phif;
  int ngen;
  int idhep[1000];
  int daF[1000];
  int daL[1000];
  int mo[1000];
} MCGenEvnt;

#endif

