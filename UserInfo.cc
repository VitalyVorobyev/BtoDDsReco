// Base on example from J. Tanaka

#include "UserInfo.h"
#include "particle/Particle.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif 

//////////////
//  Tracks  //
//////////////
TrkUserInfo::TrkUserInfo() :
  m_rz_svd_hits(0),
  m_rphi_svd_hits(0),
  m_atckpi(0.),
  m_r(0.),
  m_z(0.)
{
}

TrkUserInfo::~TrkUserInfo()
{
}

TrkUserInfo::TrkUserInfo(const TrkUserInfo &x)
  : m_rz_svd_hits(x.m_rz_svd_hits),
    m_rphi_svd_hits(x.m_rphi_svd_hits),
    m_atckpi(x.m_atckpi),
    m_r(x.m_r),
    m_z(x.m_z)
{
}

TrkUserInfo* TrkUserInfo::clone(void) const
{
  TrkUserInfo *x = new TrkUserInfo( *this );
  return x;
}

TrkUserInfo & TrkUserInfo::operator = (const TrkUserInfo &x)
{
  m_rz_svd_hits   = x.m_rz_svd_hits;
  m_rphi_svd_hits = x.m_rphi_svd_hits;
  m_atckpi        = x.m_atckpi;
  m_r             = x.m_r;
  m_z             = x.m_z;
  return *this;
}

//////////////
//    Ks    //
//////////////
KsUserInfo::KsUserInfo() :
  m_dr(0),
  m_dz(0),
  m_dphi(0.),
  m_fl(0.)
{
}

KsUserInfo::~KsUserInfo()
{
}

KsUserInfo::KsUserInfo(const KsUserInfo &x)
  : m_dr(x.m_dr),
    m_dz(x.m_dz),
    m_dphi(x.m_dphi),
    m_fl(x.m_fl)
{
}

KsUserInfo* KsUserInfo::clone(void) const
{
  KsUserInfo *x = new KsUserInfo( *this );
  return x;
}

KsUserInfo & KsUserInfo::operator = (const KsUserInfo &x)
{
  m_dr = x.m_dr;
  m_dz = x.m_dz;
  m_dphi = x.m_dphi;
  m_fl = x.m_fl;
  return *this;
}

///////////////
//    pi0    //
///////////////
Pi0UserInfo::Pi0UserInfo() :
  m_mass(0),
  m_chi2(0),
  m_fit_flag(-1),
  m_gamma_energy(0),
  m_mode(-1),
  m_gamma_pt(0)
{
}

Pi0UserInfo::~Pi0UserInfo()
{
}

Pi0UserInfo::Pi0UserInfo(const Pi0UserInfo &x)
  : m_mass(x.m_mass),
    m_chi2(x.m_chi2),
    m_gamma_energy(x.m_gamma_energy),
    m_gamma_pt(x.m_gamma_pt),
    m_fit_flag(x.m_fit_flag),
    m_mode(x.m_mode)
{
}

Pi0UserInfo* Pi0UserInfo::clone(void) const
{
  Pi0UserInfo *x = new Pi0UserInfo( *this );
  return x;
}

Pi0UserInfo & Pi0UserInfo::operator = (const Pi0UserInfo &x)
{
  m_mass = x.m_mass;
  m_chi2 = x.m_chi2;
  m_gamma_energy = x.m_gamma_energy;
  m_gamma_pt = x.m_gamma_pt;
  m_fit_flag = x.m_fit_flag;
  m_mode = x.m_mode;
  return *this;
}

///////////////
//    phi    //
///////////////
PhiUserInfo::PhiUserInfo() :
  m_mass(0),
  m_chi2(-1),
  m_mcflag(0)
{
}

PhiUserInfo::~PhiUserInfo()
{
}

PhiUserInfo::PhiUserInfo(const PhiUserInfo &x)
  : m_mass(x.m_mass),
    m_chi2(x.m_chi2),
    m_mcflag(x.m_mcflag)
{
}

PhiUserInfo* PhiUserInfo::clone(void) const
{
  PhiUserInfo *x = new PhiUserInfo( *this );
  return x;
}

PhiUserInfo & PhiUserInfo::operator = (const PhiUserInfo &x)
{
  m_mass = x.m_mass;
  m_chi2 = x.m_chi2;
  m_mcflag = x.m_mcflag;
  return *this;
}


////////////////
// D0 and Ds+ //
////////////////
DUserInfo::DUserInfo() :
  m_pcm(0),
  m_mass_chi2(-1.),
  m_m_raw(0.),
  m_fit_done(false),
  m_fit_good(false),
  m_mcflag(0)
{
}

DUserInfo::~DUserInfo()
{
}

DUserInfo::DUserInfo(const DUserInfo &x)
  : m_m_raw(x.m_m_raw),
    m_pcm(x.m_pcm),
    m_fit_done(x.m_fit_done),
    m_fit_good(m_fit_good),
    m_mass_chi2(x.m_mass_chi2),
    m_mcflag(x.m_mcflag)
{
}

DUserInfo* DUserInfo::clone(void) const
{
  DUserInfo *x = new DUserInfo( *this );
  return x;
}

DUserInfo & DUserInfo::operator = (const DUserInfo &x)
{
  m_fit_done = x.m_fit_done;
  m_fit_good = x.m_fit_good;
  m_m_raw = x.m_m_raw;
  m_pcm = x.m_pcm;
  m_mass_chi2 = x.m_mass_chi2;
  m_mcflag = x.m_mcflag;

  return *this;
}

////////////////
//     B+-    //
////////////////
BUserInfo::BUserInfo() :
  m_de(0.),
  m_mbc(0.),
  m_pcm(0.),
  m_mcflag(-2),
  m_cos_theta_cms(-2.)
{
}

BUserInfo::~BUserInfo()
{
}

BUserInfo::BUserInfo(const BUserInfo &x)
  : m_pcm(x.m_pcm),
    m_de(x.m_de),
    m_mbc(x.m_mbc),
    m_mcflag(x.m_mcflag),
    m_cos_theta_cms(x.m_cos_theta_cms)
{
}

BUserInfo* BUserInfo::clone(void) const
{
  BUserInfo *x = new BUserInfo( *this );
  return x;
}

BUserInfo & BUserInfo::operator = (const BUserInfo &x)
{
  m_pcm = x.m_pcm;
  m_de = x.m_de;
  m_mbc = x.m_mbc;
  m_mcflag = x.m_mcflag;
  m_cos_theta_cms = x.m_cos_theta_cms;
  return *this;
}

#if defined(BELLE_NAMESPACE)
}
#endif 
