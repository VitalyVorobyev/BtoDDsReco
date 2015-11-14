#include "particle/Particle.h"
#include "particle/PID.h"
#include "particle/utility.h"
#include "particle/combination.h"
#include "particle/ParticleUserInfo.h"
#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"
#include "eid/eid.h"
#include "kid/atc_pid.h"
#include "panther/panther.h"
#include "mdst/mdst.h"
#include "mdst/Muid_mdst.h"
#include "mdst/findKs.h"
#include "ip/IpProfile.h"
#include MDST_H
#include HEPEVT_H
#include BELLETDF_H
#include EVTCLS_H

#include <stdio.h>
#include <math.h>
#include <strings.h>
#include <string.h>

#include "UserInfo.h"
#include "geninfo.h"
#include "dataStructures.h"

#include "benergy/BeamEnergy.h"

#include "hamlet/Hamlet.h"
#include "hamlet/Fbtag_MultDimLikelihood0.h"
#include "hamlet/Fbtag_NN1.h"
#include "tagv/TagV.h"

#include "toolbox/FuncPtr.h"
#include "toolbox/FoxWolfr.h"
#include "toolbox/Thrust.h"

#include "tuple/BelleTupleManager.h"

#include <vector>
#include <iostream>
#include "belle.h"
#include "../B0toDh0/rooksfw/ksfwmoments.h"
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Vector/ThreeVector.h>
#include "CLHEP/Geometry/Point3D.h"
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/SymMatrix.h>

#include "nisKsFinder/nisKsFinder.h"
#include "exkfitter/ExKFitter.h"

#include "kfitter/kvertexfitter.h"
#include "kfitter/kmassfitter.h"
#include "kfitter/kmakemother.h"

#include "TTree.h"
#include "TFile.h"
#include <TObjArray.h>

using namespace std;

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class b2dds : public Module{
public:
  b2dds(void): clock(0),ntuple_flag(1),n_good_b(0) {
    cout << "Constructor" << endl;
  };
  ~b2dds(void) {};
  void init(int *);
  void term(void);
  void disp_stat(const char*){};
  void hist_def(void) {};
  void begin_run(BelleEvent*, int*);
  void event(BelleEvent*, int*);
  void end_run(BelleEvent*, int*){};

  int m_mode;//0 -> Data
             //1 -> Signal MC
             //2 -> Genegic M
  int ntuple_flag;
  char ofile[1024];

private:
  void make_kpi_list(vector<Particle> &pip_list, vector<Particle> &pim_list,vector<Particle> &kp_list, vector<Particle> &km_list);
  int make_phi_list(vector<Particle> &phi_list,vector<Particle> &kp_list,vector<Particle> &km_list);
  int make_d0_list(vector<Particle> &d0_list,vector<Particle> &pip_list,vector<Particle> &pim_list,vector<Particle> &kp_list,vector<Particle> &km_list);
  int make_dst0_list(vector<Particle> &dst0_list,vector<Particle> &d0_list,vector<Particle> &pi0_list);
  int make_ds_list(vector<Particle> &ds_list,vector<Particle> &phi_list,vector<Particle> &pip_list,vector<Particle> &pim_list);
  int make_bp_list(vector<Particle> &bp_list,vector<Particle> &d0_list,vector<Particle> &ds_list);
  int make_mc_list(void);

  int SetUserInfoTrk(Particle& pi);

  bool CheckSVD(const Particle& pi);
  int TwoBodyMassFit(Particle& p,double& chisq);

  HepPoint3D IP;
  double benergy;

  bool mc_flag;
  TTree* TEvent;
  TFile* tfile;

  bEvnt bEvt;
  MCGenEvnt MCGenEvt;
  int clock;

  int FillBEvent(Particle& bp);
  bool IsDuplicated(const Particle& bp);
  int Mbc_deltaE(const Particle& bp,double& mbc,double& de);
  double CosThr(const Particle& bp,double& thr_sig,double& thr_oth);
  int BoostedList(const Particle& bp,vector<Hep3Vector>&sigtrk,vector<Hep3Vector>&otherb);
  double CosHelAngle(const Particle& p);

  int FillMCGenEvt(Particle& bp);
  vector<Gen_hepevt> mc_list;
  int n_good_b;
};

#if defined(BELLE_NAMESPACE)
}
#endif

