// Used in conjunction with the particle class to store various useful numbers
// as part of particle.

#if !defined(USERINFO_H_INCLUDED)
#define USERINFO_H_INCLUDED

#include "CLHEP/Vector/LorentzVector.h"

#include "belle.h"
#include "particle/ParticleUserInfo.h"

//#include "helix/Helix.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif 

class TrkUserInfo : public ParticleUserInfo
{
public:
  TrkUserInfo();
  TrkUserInfo(const TrkUserInfo &);
  virtual ~TrkUserInfo();
  TrkUserInfo * clone(void) const;
  TrkUserInfo & operator = (const TrkUserInfo &);
public:
// charged pi
  void rz_svd_hits(const int v) {m_rz_svd_hits = v;}
  int rz_svd_hits(void) const {return m_rz_svd_hits;}

  void rphi_svd_hits(const int v) {m_rphi_svd_hits = v;}
  int rphi_svd_hits(void) const {return m_rphi_svd_hits;}

  void atckpi(const double& v) {m_atckpi = v;}
  double atckpi(void) const {return m_atckpi;}

  void r(const double v) {m_r = v;}
  double r(void) const {return m_r;}

  void z(const double v) {m_z = v;}
  double z(void) const {return m_z;}
private:
// charged pions
  int m_rz_svd_hits;
  int m_rphi_svd_hits;
  double m_atckpi;
  double m_r,m_z;
};

class KsUserInfo : public ParticleUserInfo
{
public:
  KsUserInfo();
  KsUserInfo(const KsUserInfo &);
  virtual ~KsUserInfo();
  KsUserInfo * clone(void) const;
  KsUserInfo & operator = (const KsUserInfo &);
public:
  void dr(const double& v) {m_dr = v;}
  double dr(void) const {return m_dr;}

  void dz(const double& v) {m_dz = v;}
  double dz(void) const {return m_dz;}

  void dphi(const double& v) {m_dphi = v;}
  double dphi(void) const {return m_dphi;}

  void fl(const double& v) {m_fl = v;}
  double fl(void) const {return m_fl;}

private:
  double m_dr,m_dz,m_dphi,m_fl;
};

class PhiUserInfo : public ParticleUserInfo
{
public:
  PhiUserInfo();
  PhiUserInfo(const PhiUserInfo &);
  virtual ~PhiUserInfo();
  PhiUserInfo * clone(void) const;
  PhiUserInfo & operator = (const PhiUserInfo &);
public:
  void Mass(const double& v) {m_mass = v;}
  double Mass(void) const {return m_mass;}

  void Chi2(const double& v) {m_chi2 = v;}
  double Chi2(void) const {return m_chi2;}

  void mcflag(const int v) {m_mcflag = v;}
  int mcflag(void) const {return m_mcflag;}

private:
  double m_mass;
  double m_chi2;
  int m_mcflag;
};

class Pi0UserInfo : public ParticleUserInfo
{
public:
  Pi0UserInfo();
  Pi0UserInfo(const Pi0UserInfo &);
  virtual ~Pi0UserInfo();
  Pi0UserInfo * clone(void) const;
  Pi0UserInfo & operator = (const Pi0UserInfo &);
public:
  void Mass(const double& v) {m_mass = v;}
  double Mass(void) const {return m_mass;}

  void Chi2(const double& v) {m_chi2 = v;}
  double Chi2(void) const {return m_chi2;}

  void FitFlag(const int v) {m_fit_flag = v;}
  int FitFlag(void) const {return m_fit_flag;}

  void Mode(const int v) {m_mode = v;}
  int Mode(void) const {return m_mode;}

  void EGamma(const double& v) {m_gamma_energy = v;}
  double EGamma(void) {return m_gamma_energy;}

  void PtGamma(const double& v) {m_gamma_pt = v;}
  double PtGamma(void) {return m_gamma_pt;}

private:
  double m_mass;
  double m_chi2;
  double m_gamma_energy;
  double m_gamma_pt;
  int m_fit_flag;
  int m_mode;
};

class DUserInfo : public ParticleUserInfo
{
public:
  DUserInfo();
  DUserInfo(const DUserInfo &);
  virtual ~DUserInfo();
  DUserInfo * clone(void) const;
  DUserInfo & operator = (const DUserInfo &);
public:
  void pcm(const double& v) {m_pcm = v;}
  double pcm(void) const {return m_pcm;}

  void m_raw(const double& v) {m_m_raw = v;}
  double m_raw(void) const {return m_m_raw;}

  void FitDone(void) {m_fit_done = true;}
  bool IsFitDone(void) const {return m_fit_done;}

  void FitGood(const bool v) {m_fit_good = v;}
  bool IsFitGood(void) const {return m_fit_good;}

  void MassChi2(const double& v) {m_mass_chi2 = v;}
  double MassChi2(void) const {return m_mass_chi2;}

  void mcflag(const int v) {m_mcflag = v;}
  int mcflag(void) const {return m_mcflag;}

private:
  double m_pcm;
  double m_m_raw;
  double m_mass_chi2;
  bool m_fit_done;
  bool m_fit_good;
  int m_mcflag;
};

class BUserInfo : public ParticleUserInfo
{
public:
  BUserInfo();
  BUserInfo(const BUserInfo &);
  virtual ~BUserInfo();
  BUserInfo * clone(void) const;
  BUserInfo & operator = (const BUserInfo &);
public:
  void pcm(const double& v) {m_pcm = v;}
  double pcm(void) const {return m_pcm;}

  void deltaE(const double& v) {m_de = v;}
  double deltaE(void) const {return m_de;}

  void Mbc(const double& v) {m_mbc = v;}
  double Mbc(void) const {return m_mbc;}

  void CosThetaCMS(const double& v) {m_cos_theta_cms = v;}
  double CosThetaCMS(void) const {return m_cos_theta_cms;}

//////////
//  MC  //
//////////
  void mcflag(const int v) {m_mcflag = v;}
  int mcflag(void) const {return m_mcflag;}

private:
  double m_pcm;
  double m_de;
  double m_mbc;
  double m_cos_theta_cms;
  int m_mcflag;
};

#if defined(BELLE_NAMESPACE)
}
#endif 

#endif
