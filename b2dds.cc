#include "b2dds.h"
#include "cuts.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

const bool dump = false;
extern "C" Module_descr *mdcl_b2dds()
{
  b2dds *module = new b2dds;
  Module_descr *dscr = new Module_descr("b2dds",module);
  BeamEnergy::define_global(dscr);

  dscr->define_param("mode","mode",&module->m_mode);
  dscr->define_param("ofile","ofile",1024,module->ofile);
  dscr->define_param("ntuple_flag","ntuple_flag",&module->ntuple_flag);

  return dscr;
}

void b2dds::init(int *){
  std::cout << "init" << endl;
  Ptype dummy("VPHO");
  if(ntuple_flag) Hamlet::init();
  return;
}

void b2dds::begin_run(BelleEvent* evptr, int *status){
  if(dump) cout << "Begin Run" << endl;
  IpProfile::begin_run();
  BeamEnergy::begin_run();
  if(ntuple_flag) Hamlet::begin_run(Hamlet::MULT_DIM_LH);
  Belle_runhead_Manager &rhd_mgr = Belle_runhead_Manager::get_manager();
  std::vector<Belle_runhead>::const_iterator rhd = rhd_mgr.begin();
  if(rhd == rhd_mgr.end()){
    fprintf(stderr,"Constructor: Cannot access to Belle_runhead\n");
    mc_flag = false;
  } else{
    if(rhd->ExpMC() == 1) mc_flag = false;//real data
    else mc_flag = true;
  }
  if(dump) cout << "Begin Run " << ntuple_flag << " " << clock << endl;
  if(ntuple_flag && !clock){
    if(dump) cout << "ntuple_flag!" << endl;
    clock = 1;
    cout << ofile << endl;
    tfile = new TFile(ofile,"RECREATE");
    cout << ofile << endl;

    TEvent = new TTree("TEvent","TEvent");
    TEvent->Branch("exp",&bEvt.exp,"exp/I");
    TEvent->Branch("run",&bEvt.run,"run/I");
    TEvent->Branch("evtn",&bEvt.evtn,"evtn/I");
    TEvent->Branch("charge",&bEvt.charge,"charge/I");
    TEvent->Branch("de",&bEvt.de,"de/D");
    TEvent->Branch("mbc",&bEvt.mbc,"mbc/D");

    TEvent->Branch("costhBcms",&bEvt.costhBcms,"costhBcms/D");

    TEvent->Branch("md0",&bEvt.md0,"md0/D");
    TEvent->Branch("md0_raw",&bEvt.md0_raw,"md0_raw/D");
    TEvent->Branch("mds",&bEvt.mds,"mds/D");
    TEvent->Branch("mds_raw",&bEvt.mds_raw,"mds_raw/D");
    TEvent->Branch("mphi",&bEvt.mphi,"mphi/D");
    TEvent->Branch("mphi_raw",&bEvt.mphi_raw,"mphi_raw/D");

    TEvent->Branch("pcm_d0",&bEvt.pcm_d0,"pcm_d0/D");
    TEvent->Branch("pcm_ds",&bEvt.pcm_ds,"pcm_ds/D");

    TEvent->Branch("cos_thr",&bEvt.cos_thr,"cos_thr/D");
    TEvent->Branch("thr_sig",&bEvt.thr_sig,"thr_sig/D");
    TEvent->Branch("thr_oth",&bEvt.thr_oth,"thr_oth/D");

    TEvent->Branch("cos_hel_d0",&bEvt.cos_hel_d0,"cos_hel_d0/D");
    TEvent->Branch("cos_hel_ds",&bEvt.cos_hel_ds,"cos_hel_ds/D");
    TEvent->Branch("cos_hel_phi",&bEvt.cos_hel_phi,"cos_hel_phi/D");

    TEvent->Branch("px_k_d0",&bEvt.k_d0.px,"px_k_d0/D");
    TEvent->Branch("py_k_d0",&bEvt.k_d0.py,"py_k_d0/D");
    TEvent->Branch("pz_k_d0",&bEvt.k_d0.pz,"pz_k_d0/D");
    TEvent->Branch("r_k_d0",&bEvt.k_d0.r,"r_k_d0/D");
    TEvent->Branch("z_k_d0",&bEvt.k_d0.z,"z_k_d0/D");
    TEvent->Branch("rz_svd_k_d0",&bEvt.k_d0.rz_svd,"rz_svd_k_d0/I");
    TEvent->Branch("rphi_svd_k_d0",&bEvt.k_d0.rphi_svd,"rphi_svd_k_d0/I");
    TEvent->Branch("atckpi_k_d0",&bEvt.k_d0.atckpi,"atckpi_k_d0/D");

    TEvent->Branch("px_pi_d0",&bEvt.pi_d0.px,"px_pi_d0/D");
    TEvent->Branch("py_pi_d0",&bEvt.pi_d0.py,"py_pi_d0/D");
    TEvent->Branch("pz_pi_d0",&bEvt.pi_d0.pz,"pz_pi_d0/D");
    TEvent->Branch("r_pi_d0",&bEvt.pi_d0.r,"r_pi_d0/D");
    TEvent->Branch("z_pi_d0",&bEvt.pi_d0.z,"z_pi_d0/D");
    TEvent->Branch("rz_svd_pi_d0",&bEvt.pi_d0.rz_svd,"rz_svd_pi_d0/I");
    TEvent->Branch("rphi_svd_pi_d0",&bEvt.pi_d0.rphi_svd,"rphi_svd_pi_d0/I");
    TEvent->Branch("atckpi_pi_d0",&bEvt.pi_d0.atckpi,"atckpi_pi_d0/D");

    TEvent->Branch("px_pi_ds",&bEvt.pi_ds.px,"px_pi_ds/D");
    TEvent->Branch("py_pi_ds",&bEvt.pi_ds.py,"py_pi_ds/D");
    TEvent->Branch("pz_pi_ds",&bEvt.pi_ds.pz,"pz_pi_ds/D");
    TEvent->Branch("r_pi_ds",&bEvt.pi_ds.r,"r_pi_ds/D");
    TEvent->Branch("z_pi_ds",&bEvt.pi_ds.r,"z_pi_ds/D");
    TEvent->Branch("rz_svd_pi_ds",&bEvt.pi_ds.rz_svd,"rz_svd_pi_ds/I");
    TEvent->Branch("rphi_svd_pi_ds",&bEvt.pi_ds.rphi_svd,"rphi_svd_pi_ds/I");
    TEvent->Branch("atckpi_pi_ds",&bEvt.pi_ds.atckpi,"atckpi_pi_ds/D");

    TEvent->Branch("px_kp_phi",&bEvt.kp_phi.px,"px_kp_phi/D");
    TEvent->Branch("py_kp_phi",&bEvt.kp_phi.py,"py_kp_phi/D");
    TEvent->Branch("pz_kp_phi",&bEvt.kp_phi.pz,"pz_kp_phi/D");
    TEvent->Branch("r_kp_phi", &bEvt.kp_phi.r,"r_pi_ds/D");
    TEvent->Branch("z_kp_phi", &bEvt.kp_phi.z,"z_kp_phi/D");
    TEvent->Branch("rz_svd_kp_phi",  &bEvt.kp_phi.rz_svd,"rz_svd_kp_phi/I");
    TEvent->Branch("rphi_svd_kp_phi",&bEvt.kp_phi.rphi_svd,"rphi_svd_kp_phi/I");
    TEvent->Branch("atckpi_kp_phi",  &bEvt.kp_phi.atckpi,"atckpi_kp_phi/D");

    TEvent->Branch("px_km_phi",&bEvt.km_phi.px,"px_km_phi/D");
    TEvent->Branch("py_km_phi",&bEvt.km_phi.py,"py_km_phi/D");
    TEvent->Branch("pz_km_phi",&bEvt.km_phi.pz,"pz_km_phi/D");
    TEvent->Branch("r_km_phi", &bEvt.km_phi.r,"r_pi_ds/D");
    TEvent->Branch("z_km_phi", &bEvt.km_phi.z,"z_kp_phi/D");
    TEvent->Branch("rz_svd_km_phi",  &bEvt.km_phi.rz_svd,"rz_svd_km_phi/I");
    TEvent->Branch("rphi_svd_km_phi",&bEvt.km_phi.rphi_svd,"rphi_svd_km_phi/I");
    TEvent->Branch("atckpi_km_phi",  &bEvt.km_phi.atckpi,"atckpi_km_phi/D");

    TEvent->Branch("chi2_fit_d0",&bEvt.chi2_fit_d0,"chi2_fit_d0/D");
    TEvent->Branch("chi2_fit_ds",&bEvt.chi2_fit_ds,"chi2_fit_ds/D");
    TEvent->Branch("chi2_fit_phi",&bEvt.chi2_fit_phi,"chi2_fit_phi/D");

    TEvent->Branch("ipx",&bEvt.ipx,"ipx/D");
    TEvent->Branch("ipy",&bEvt.ipy,"ipy/D");
    TEvent->Branch("ipz",&bEvt.ipz,"ipz/D");

    TEvent->Branch("boostx",&bEvt.boostx,"boostx/D");
    TEvent->Branch("boosty",&bEvt.boosty,"boosty/D");
    TEvent->Branch("boostz",&bEvt.boostz,"boostz/D");

    TEvent->Branch("k0mm2",&bEvt.k0mm2,"k0mm2/D");
    TEvent->Branch("k0et",&bEvt.k0et,"k0et/D");
    TEvent->Branch("k0hso00",&bEvt.k0hso00,"k0hso00/D");
    TEvent->Branch("k0hso01",&bEvt.k0hso01,"k0hso01/D");
    TEvent->Branch("k0hso02",&bEvt.k0hso02,"k0hso02/D");
    TEvent->Branch("k0hso03",&bEvt.k0hso03,"k0hso03/D");
    TEvent->Branch("k0hso04",&bEvt.k0hso04,"k0hso04/D");
    TEvent->Branch("k0hso10",&bEvt.k0hso10,"k0hso10/D");
    TEvent->Branch("k0hso12",&bEvt.k0hso12,"k0hso12/D");
    TEvent->Branch("k0hso14",&bEvt.k0hso14,"k0hso14/D");
    TEvent->Branch("k0hso20",&bEvt.k0hso20,"k0hso20/D");
    TEvent->Branch("k0hso22",&bEvt.k0hso22,"k0hso22/D");
    TEvent->Branch("k0hso24",&bEvt.k0hso24,"k0hso24/D");
    TEvent->Branch("k0hoo0",&bEvt.k0hoo0,"k0hoo0/D");
    TEvent->Branch("k0hoo1",&bEvt.k0hoo1,"k0hoo1/D");
    TEvent->Branch("k0hoo2",&bEvt.k0hoo2,"k0hoo2/D");
    TEvent->Branch("k0hoo3",&bEvt.k0hoo3,"k0hoo3/D");
    TEvent->Branch("k0hoo4",&bEvt.k0hoo4,"k0hoo4/D");

    TEvent->Branch("k1mm2",&bEvt.k1mm2,"k1mm2/D");
    TEvent->Branch("k1et",&bEvt.k1et,"k1et/D");
    TEvent->Branch("k1hso00",&bEvt.k1hso00,"k1hso00/D");
    TEvent->Branch("k1hso01",&bEvt.k1hso01,"k1hso01/D");
    TEvent->Branch("k1hso02",&bEvt.k1hso02,"k1hso02/D");
    TEvent->Branch("k1hso03",&bEvt.k1hso03,"k1hso03/D");
    TEvent->Branch("k1hso04",&bEvt.k1hso04,"k1hso04/D");
    TEvent->Branch("k1hso10",&bEvt.k1hso10,"k1hso10/D");
    TEvent->Branch("k1hso12",&bEvt.k1hso12,"k1hso12/D");
    TEvent->Branch("k1hso14",&bEvt.k1hso14,"k1hso14/D");
    TEvent->Branch("k1hso20",&bEvt.k1hso20,"k1hso20/D");
    TEvent->Branch("k1hso22",&bEvt.k1hso22,"k1hso22/D");
    TEvent->Branch("k1hso24",&bEvt.k1hso24,"k1hso24/D");
    TEvent->Branch("k1hoo0",&bEvt.k1hoo0,"k1hoo0/D");
    TEvent->Branch("k1hoo1",&bEvt.k1hoo1,"k1hoo1/D");
    TEvent->Branch("k1hoo2",&bEvt.k1hoo2,"k1hoo2/D");
    TEvent->Branch("k1hoo3",&bEvt.k1hoo3,"k1hoo3/D");
    TEvent->Branch("k1hoo4",&bEvt.k1hoo4,"k1hoo4/D");

    if(m_mode == 2 || m_mode == 1){
      TEvent->Branch("pi_d0_ch", MCGenEvt.pi_d0_chain,"pi_d0_ch[9]/I");
      TEvent->Branch("pi_ds_ch", MCGenEvt.pi_ds_chain,"pi_ds_ch[9]/I");
      TEvent->Branch("k_d0_ch",  MCGenEvt.k_d0_chain,"k_d0_ch[9]/I");
      TEvent->Branch("phi_ch",   MCGenEvt.phi_chain,"phi_ch[9]/I");
      TEvent->Branch("kp_phi_ch",MCGenEvt.kp_phi_chain,"kp_phi_ch[9]/I");
      TEvent->Branch("km_phi_ch",MCGenEvt.km_phi_chain,"km_phi_ch[9]/I");
      TEvent->Branch("d0_ch",    MCGenEvt.d0_chain,"d0_ch[9]/I");
      TEvent->Branch("ds_ch",    MCGenEvt.ds_chain,"ds_ch[9]/I");

      TEvent->Branch("bpid",&MCGenEvt.bpid,"bpid/I");
      TEvent->Branch("bpf",&MCGenEvt.bpf,"bpf/I");
      TEvent->Branch("d0id",&MCGenEvt.d0id,"d0id/I");
      TEvent->Branch("d0f",&MCGenEvt.d0f,"d0f/I");
      TEvent->Branch("dsid",&MCGenEvt.d0id,"dsid/I");
      TEvent->Branch("dsf",&MCGenEvt.d0f,"dsf/I");
      TEvent->Branch("phiid",&MCGenEvt.phiid,"phiid/I");
      TEvent->Branch("phi0f",&MCGenEvt.phif,"phif/I");

      TEvent->Branch("ngen",&MCGenEvt.ngen,"ngen/I");
      TEvent->Branch("gen_idhep",MCGenEvt.idhep,"gen_idhep[ngen]/I");
      TEvent->Branch("gen_daF",MCGenEvt.daF,"gen_daF[ngen]/I");
      TEvent->Branch("gen_daL",MCGenEvt.daL,"gen_daL[ngen]/I");
      TEvent->Branch("gen_mo",MCGenEvt.mo,"gen_mo[ngen]/I");
    }
  }

  return;
}

void b2dds::term(void){
  cout << "Term " << ofile << endl;
  if(ntuple_flag){
    cout << "Term TEvent: " << TEvent->GetEntries() << endl;
    TEvent->Write();
    tfile->Close();
  }
  cout << "# of Good B0: " << n_good_b << endl;
}

void b2dds::event(BelleEvent* evptr, int *status){
  *status = 0;
  Belle_event_Manager& evtmgr = Belle_event_Manager::get_manager();
  Belle_event& evthead = *evtmgr.begin();
  if(!dump && !((evthead.EvtNo() & 0xfffffff)%100)){
    cout << "run: " << evthead.RunNo() << ", evtn: " << (evthead.EvtNo() & 0xfffffff) << endl;
  }
  if(IpProfile::usable()){
    IP = IpProfile::e_position();
  } else {
    std::cout<<"IP profile is not available!!!\n";
    return;
  }

  vector<Particle> pip_list;
  vector<Particle> pim_list;
  vector<Particle> kp_list;
  vector<Particle> km_list;

  vector<Particle> phi_list;
  vector<Particle> d0_list;
  vector<Particle> ds_list;
  vector<Particle> bp_list;

  make_kpi_list(pip_list,pim_list,kp_list,km_list);
  if(!pip_list.size() || !pim_list.size()) return;
  make_phi_list(phi_list,kp_list,km_list);
  if(!phi_list.size()) return;

  make_ds_list(ds_list,phi_list,pip_list,pim_list);

  make_d0_list(d0_list,pip_list,pim_list,kp_list,km_list);
  if(!d0_list.size()) return;

  const int goodb = make_bp_list(bp_list,d0_list,ds_list);
  n_good_b += goodb;

  if(ntuple_flag){
    if(m_mode) make_mc_list();
    for(int i=0; i<bp_list.size(); i++){
      FillBEvent(bp_list[i]);
      if(m_mode == 2 || m_mode == 1){
        FillMCGenEvt(bp_list[i]);
      }
      TEvent->Fill();
    }
  } else{
    *status = goodb == 0 ? 0 : 1;
  }
  return;
}

int b2dds::make_mc_list(void){
  mc_list.clear();
  Gen_hepevt_Manager &mGHep = Gen_hepevt_Manager::get_manager();
  for(std::vector<Gen_hepevt>::iterator i = mGHep.begin();i != mGHep.end();i++){
    Gen_hepevt& tmp = *i;
//    cout << tmp.idhep() << " " << tmp.daFirst() << " " << tmp.daLast() << endl;
    if(!tmp) continue;
    mc_list.push_back(tmp);
  }
  return mc_list.size();
}

void b2dds::make_kpi_list(vector<Particle> &pip_list, vector<Particle> &pim_list,vector<Particle> &kp_list, vector<Particle> &km_list){
  pip_list.clear(); pim_list.clear(); kp_list.clear(); km_list.clear();

  const Gen_hepevt &null = Gen_hepevt_Manager::get_manager().get_NULL();

  Mdst_charged_Manager &chg_mgr = Mdst_charged_Manager::get_manager();
  for(vector<Mdst_charged>::const_iterator i = chg_mgr.begin(); i != chg_mgr.end(); i++){
    const Mdst_charged& track = *i;
    Particle pi(track,Ptype(track.charge()>0 ? "PI+" : "PI-"));
    Particle K(track,Ptype(track.charge()>0  ? "K+"  : "K-" ));

    if(ntuple_flag){
      pi.userInfo(TrkUserInfo());
      K.userInfo(TrkUserInfo());
    }
    const Gen_hepevt &h = (mc_flag == 1) ? get_hepevt(track) : null;
    if(h){
      pi.relation().genHepevt(h);
      K.relation().genHepevt(h);
    }
    if(ntuple_flag){
      SetUserInfoTrk(pi);
      SetUserInfoTrk(K);
    }
    if(track.charge() > 0){
      pip_list.push_back(pi);
      kp_list.push_back(K);
    } else{
      pim_list.push_back(pi);
      km_list.push_back(K);
    }
  }

  return;
}

int b2dds::make_phi_list(vector<Particle> &phi_list,vector<Particle> &kp_list,vector<Particle> &km_list){
  phi_list.clear();
  combination(phi_list,Ptype("PHI"),kp_list,km_list);
  withMassCut(phi_list,mphi_min,mphi_max);
  if(ntuple_flag){
    for(int i=0; i<phi_list.size(); i++){
      phi_list[i].userInfo(PhiUserInfo());
      dynamic_cast<PhiUserInfo&>(phi_list[i].userInfo()).Mass(phi_list[i].p().m());
      if(mc_flag == 1){
        setMCtruth(phi_list[i]);
        dynamic_cast<PhiUserInfo&>(phi_list[i].userInfo()).mcflag(getMCtruthFlag(phi_list[i]));
      }
      double chisq;
      TwoBodyMassFit(phi_list[i],chisq);
      dynamic_cast<PhiUserInfo&>(phi_list[i].userInfo()).Chi2(chisq);
    }
  }
  return phi_list.size();
}

int b2dds::make_ds_list(vector<Particle> &ds_list,vector<Particle> &phi_list,vector<Particle> &pip_list,vector<Particle> &pim_list){
  ds_list.clear();
  combination(ds_list,Ptype("DS+"),phi_list,pip_list);
  withMassCut(ds_list,mds_min,mds_max);

  vector<Particle> dsm_list;
  combination(dsm_list,Ptype("DS-"),phi_list,pim_list);
  withMassCut(dsm_list,mds_min,mds_max);
  ds_list.insert(ds_list.end(),dsm_list.begin(),dsm_list.end());

  if(ntuple_flag){
    for(int i=0; i<ds_list.size(); i++){
      ds_list[i].userInfo(DUserInfo());
      dynamic_cast<DUserInfo&>(ds_list[i].userInfo()).m_raw(ds_list[i].p().m());
      if(mc_flag == 1){
        setMCtruth(ds_list[i]);
        dynamic_cast<DUserInfo&>(ds_list[i].userInfo()).mcflag(getMCtruthFlag(ds_list[i]));
      }
      double chisq;
      TwoBodyMassFit(ds_list[i],chisq);
      dynamic_cast<DUserInfo&>(ds_list[i].userInfo()).MassChi2(chisq);
    }
  }
  return ds_list.size();
}


bool IsTheSame(const Particle& p1, const Particle& p2){
  const Hep3Vector dp = p1.p().vect() - p2.p().vect();
  return dp.mag() < 0.0001 ? true : false;
}

int b2dds::make_d0_list(vector<Particle> &d0_list,vector<Particle> &pip_list,vector<Particle> &pim_list,vector<Particle> &kp_list,vector<Particle> &km_list){
  d0_list.clear();
  combination(d0_list,Ptype("D0"),km_list,pip_list);
  withMassCut(d0_list,md0_min,md0_max);

  vector<Particle> d0b_list;
  combination(d0b_list,Ptype("D0B"),kp_list,pim_list);
  withMassCut(d0b_list,md0_min,md0_max);
  d0_list.insert(d0_list.end(),d0b_list.begin(),d0b_list.end());

  if(ntuple_flag){
    for(int i=0; i<d0_list.size(); i++){
      d0_list[i].userInfo(DUserInfo());
      dynamic_cast<DUserInfo&>(d0_list[i].userInfo()).m_raw(d0_list[i].p().m());
      if(mc_flag == 1){
        setMCtruth(d0_list[i]);
        dynamic_cast<DUserInfo&>(d0_list[i].userInfo()).mcflag(getMCtruthFlag(d0_list[i]));
      }
      double chisq;
      TwoBodyMassFit(d0_list[i],chisq);
      dynamic_cast<DUserInfo&>(d0_list[i].userInfo()).MassChi2(chisq);
    }
  }
  return d0_list.size();
}

bool b2dds::IsDuplicated(const Particle& b0){
  const Particle& d0 = b0.child(0);
  const Particle& ds = b0.child(1);
  const Particle& phi = ds.child(0);
  const Particle& pi_ds = ds.child(1);
  const Particle& kp_phi = phi.child(0);
  const Particle& km_phi = phi.child(1);
  const Particle& k_d0 = d0.child(0);
  const Particle& pi_d0 = d0.child(1);

  if(IsTheSame(pi_d0,pi_ds)  || IsTheSame(pi_d0,kp_phi) || IsTheSame(pi_d0,km_phi)) return true;
  if(IsTheSame(k_d0,pi_ds)   || IsTheSame(k_d0,kp_phi)  || IsTheSame(k_d0,km_phi))  return true;
  if(IsTheSame(pi_ds,kp_phi) || IsTheSame(pi_ds,km_phi)) return true;

  return false;
}

int b2dds::make_bp_list(vector<Particle> &bp_list,vector<Particle> &d0_list,vector<Particle> &ds_list){
  bp_list.clear();
  double de,mbc;

  combination(bp_list,Ptype("B+"),d0_list,ds_list);
  for(int i=0; i<bp_list.size(); i++){
    Particle& Bp = bp_list[i];
    Mbc_deltaE(Bp,mbc,de);
    if(!ntuple_flag && fabs(de)<0.32 && mbc>5.18 && mbc<5.32){
      return 1;
    } else{
      if(IsDuplicated(Bp) || de>0.3 || de<-0.3 || mbc<5.2 || mbc>5.3){
        bp_list.erase(bp_list.begin()+i); i--;
        continue;
      }
      Bp.userInfo(BUserInfo());
      dynamic_cast<BUserInfo&>(Bp.userInfo()).deltaE(de);
      dynamic_cast<BUserInfo&>(Bp.userInfo()).Mbc(mbc);
      HepLorentzVector b_pcm = Bp.p();
      b_pcm.boost(-BeamEnergy::CMBoost());
      dynamic_cast<BUserInfo&>(Bp.userInfo()).CosThetaCMS(b_pcm.cosTheta());
    }
  }
  return bp_list.size();
}

bool b2dds::CheckSVD(const Particle& pi){
  Mdst_trk_fit &trk_fit = pi.mdstCharged().trk().mhyp(2);//2 -> pi
  if(trk_fit.nhits(3)<1 || trk_fit.nhits(4)<2) return false;// 3 -> rphi, 4 -> z
  return true;
}

int b2dds::SetUserInfoTrk(Particle& trk){
  const int hypo = abs(trk.lund()) == 211 ? 2 : 3;
  Mdst_trk_fit &trk_fit = trk.mdstCharged().trk().mhyp(hypo);//2 -> pi, 3 -> K
  dynamic_cast<TrkUserInfo&>(trk.userInfo()).rphi_svd_hits(trk_fit.nhits(3));
  dynamic_cast<TrkUserInfo&>(trk.userInfo()).rz_svd_hits(trk_fit.nhits(4));

  HepPoint3D pivot(trk_fit.pivot(0),trk_fit.pivot(1),trk_fit.pivot(2));
  HepVector a(5);
  for(int i=0;i<5;i++) a[i] = trk_fit.helix(i);

  Helix helix(pivot,a);
  helix.pivot(IpProfile::e_position());

  dynamic_cast<TrkUserInfo&>(trk.userInfo()).r(helix.dr());
  dynamic_cast<TrkUserInfo&>(trk.userInfo()).z(helix.dz());

  atc_pid selKpi(3,1,5,3,2);
  dynamic_cast<TrkUserInfo&>(trk.userInfo()).atckpi(selKpi.prob(trk.mdstCharged()));

  return 0;
}

/////////////////////////////////////////////////////////////////////
//                     Fill Generic MC Event                       //
/////////////////////////////////////////////////////////////////////
int b2dds::FillMCGenEvt(Particle& bp){
  setMCtruth(bp);

  Particle& d0     = const_cast<Particle&>(bp.child(0));
  Particle& ds     = const_cast<Particle&>(bp.child(1));
  Particle& k_d0   = const_cast<Particle&>(d0.child(0));
  Particle& pi_d0  = const_cast<Particle&>(d0.child(1));
  Particle& phi    = const_cast<Particle&>(ds.child(0));
  Particle& pi_ds  = const_cast<Particle&>(ds.child(1));
  Particle& kp_phi = const_cast<Particle&>(phi.child(0));
  Particle& km_phi = const_cast<Particle&>(phi.child(1));

  genDecayChain(ds,MCGenEvt.ds_chain);
  genDecayChain(d0,MCGenEvt.d0_chain);
  genDecayChain(phi,MCGenEvt.phi_chain);
  genDecayChain(pi_d0,MCGenEvt.pi_d0_chain);
  genDecayChain(pi_ds,MCGenEvt.pi_ds_chain);
  genDecayChain(kp_phi,MCGenEvt.kp_phi_chain);
  genDecayChain(km_phi,MCGenEvt.km_phi_chain);

  MCGenEvt.d0id  = IDhep(d0);
  MCGenEvt.dsf   = getMCtruthFlag(ds);
  MCGenEvt.phiid = IDhep(phi);
  MCGenEvt.phif  = getMCtruthFlag(phi);
  MCGenEvt.bpid  = IDhep(bp);
  MCGenEvt.bpf   = getMCtruthFlag(bp);

  MCGenEvt.ngen = mc_list.size();
  for(int i=0; i<MCGenEvt.ngen; i++){
    MCGenEvt.idhep[i] = mc_list[i].idhep();
    MCGenEvt.daF[i]   = mc_list[i].daFirst();
    MCGenEvt.daL[i]   = mc_list[i].daLast();
    MCGenEvt.mo[i]    = mc_list[i].mother();
  }
  return 0;
}

//////////////////////////////////////////////////////////////////////
//                          Fill B Event                            //
//////////////////////////////////////////////////////////////////////

void SetTrackInfo(TrkInfo& trkinf,Particle& trk){
  TrkUserInfo& info = (TrkUserInfo&)trk.userInfo();
  trkinf.px = trk.px();
  trkinf.py = trk.py();
  trkinf.pz = trk.pz();
  trkinf.r  = info.r();
  trkinf.z  = info.z();
  trkinf.rz_svd   = info.rz_svd_hits();
  trkinf.rphi_svd = info.rphi_svd_hits();
  trkinf.atckpi   = info.atckpi();
}

int b2dds::FillBEvent(Particle& bp){
  Particle& d0     = bp.child(0);
  Particle& ds     = bp.child(1);
  Particle& k_d0   = d0.child(0);
  Particle& pi_d0  = d0.child(1);
  Particle& phi    = ds.child(0);
  Particle& pi_ds  = ds.child(1);
  Particle& kp_phi = phi.child(0);
  Particle& km_phi = phi.child(1);

  BUserInfo& bpinfo    = (BUserInfo&)bp.userInfo();
  DUserInfo& d0info    = (DUserInfo&)d0.userInfo();
  DUserInfo& dsinfo    = (DUserInfo&)ds.userInfo();
  PhiUserInfo& phiinfo = (PhiUserInfo&)phi.userInfo();
  TrkUserInfo& k_d0_info   = (TrkUserInfo&)k_d0.userInfo();
  TrkUserInfo& pi_d0_info  = (TrkUserInfo&)pi_d0.userInfo();
  TrkUserInfo& pi_ds_info  = (TrkUserInfo&)pi_ds.userInfo();
  TrkUserInfo& kp_phi_info = (TrkUserInfo&)kp_phi.userInfo();
  TrkUserInfo& km_phi_info = (TrkUserInfo&)km_phi.userInfo();

  bEvt.mbc = bpinfo.Mbc();
  bEvt.de  = bpinfo.deltaE();
  bEvt.costhBcms = bpinfo.CosThetaCMS();
  bEvt.charge = (int)bp.charge();

  SetTrackInfo(bEvt.k_d0,k_d0);
  SetTrackInfo(bEvt.pi_d0,pi_d0);
  SetTrackInfo(bEvt.pi_ds,pi_ds);
  SetTrackInfo(bEvt.kp_phi,kp_phi);
  SetTrackInfo(bEvt.km_phi,km_phi);

  bEvt.ipx = IP.x();
  bEvt.ipy = IP.y();
  bEvt.ipz = IP.z();

  bEvt.boostx = -BeamEnergy::CMBoost().x();
  bEvt.boosty = -BeamEnergy::CMBoost().y();
  bEvt.boostz = -BeamEnergy::CMBoost().z();

  HepLorentzVector d0_pcm = d0.p();
  d0_pcm.boost(-BeamEnergy::CMBoost());
  bEvt.md0 = d0.mass();
  bEvt.md0_raw = d0info.m_raw();
  bEvt.pcm_d0 = d0_pcm.mag();
  bEvt.chi2_fit_d0 = d0info.MassChi2();
  bEvt.cos_hel_d0 = CosHelAngle(d0);

  HepLorentzVector ds_pcm = ds.p();
  ds_pcm.boost(-BeamEnergy::CMBoost());
  bEvt.mds = ds.mass();
  bEvt.mds_raw = dsinfo.m_raw();
  bEvt.pcm_ds = ds_pcm.mag();
  bEvt.chi2_fit_ds = dsinfo.MassChi2();
  bEvt.cos_hel_ds = CosHelAngle(ds);

  bEvt.mphi = phi.mass();
  bEvt.mphi_raw = phiinfo.Mass();
  bEvt.chi2_fit_phi = phiinfo.Chi2();
  bEvt.cos_hel_phi = CosHelAngle(phi);

  Belle_event_Manager& evtmgr = Belle_event_Manager::get_manager();
  Belle_event& evthead = *evtmgr.begin();
  bEvt.exp = evthead.ExpNo();
  bEvt.run = evthead.RunNo();
  bEvt.evtn = (evthead.EvtNo() & 0xfffffff);

  bEvt.cos_thr = CosThr(bp,bEvt.thr_sig,bEvt.thr_oth);

/////////////////////
//     RooKSFW     //
/////////////////////
  ksfwmoments km(bp,BeamEnergy::E_beam_corr(),-BeamEnergy::CMBoost());//BeamEnergy::Ecm()/2.
  km.usefinal(0);
  bEvt.k0mm2   = km.mm2();
  bEvt.k0et    = km.et();
  bEvt.k0hso00 = km.Hso(0,0);
  bEvt.k0hso01 = km.Hso(0,1);
  bEvt.k0hso02 = km.Hso(0,2);
  bEvt.k0hso03 = km.Hso(0,3);
  bEvt.k0hso04 = km.Hso(0,4);
  bEvt.k0hso10 = km.Hso(1,0);
  bEvt.k0hso12 = km.Hso(1,2);
  bEvt.k0hso14 = km.Hso(1,4);
  bEvt.k0hso20 = km.Hso(2,0);
  bEvt.k0hso22 = km.Hso(2,2);
  bEvt.k0hso24 = km.Hso(2,4);
  bEvt.k0hoo0  = km.Hoo(0);
  bEvt.k0hoo1  = km.Hoo(1);
  bEvt.k0hoo2  = km.Hoo(2);
  bEvt.k0hoo3  = km.Hoo(3);
  bEvt.k0hoo4  = km.Hoo(4);

  km.usefinal(1);
  bEvt.k1mm2   = km.mm2();
  bEvt.k1et    = km.et();
  bEvt.k1hso00 = km.Hso(0,0);
  bEvt.k1hso01 = km.Hso(0,1);
  bEvt.k1hso02 = km.Hso(0,2);
  bEvt.k1hso03 = km.Hso(0,3);
  bEvt.k1hso04 = km.Hso(0,4);
  bEvt.k1hso10 = km.Hso(1,0);
  bEvt.k1hso12 = km.Hso(1,2);
  bEvt.k1hso14 = km.Hso(1,4);
  bEvt.k1hso20 = km.Hso(2,0);
  bEvt.k1hso22 = km.Hso(2,2);
  bEvt.k1hso24 = km.Hso(2,4);
  bEvt.k1hoo0  = km.Hoo(0);
  bEvt.k1hoo1  = km.Hoo(1);
  bEvt.k1hoo2  = km.Hoo(2);
  bEvt.k1hoo3  = km.Hoo(3);
  bEvt.k1hoo4  = km.Hoo(4);

  return 0;
}

double b2dds::CosThr(const Particle& bp,double& thr_sig,double& thr_oth){
  vector<Hep3Vector> sigtrk,otherb;
  BoostedList(bp,sigtrk,otherb);

  Hep3Vector thrSig = thrust(sigtrk.begin(),sigtrk.end(),SelfFunc(Hep3Vector()));
  Hep3Vector thrOth = thrust(otherb.begin(),otherb.end(),SelfFunc(Hep3Vector()));
  thr_sig = thrSig.mag();
  thr_oth = thrOth.mag();

  return thrSig.unit().dot(thrOth.unit());
}

int b2dds::BoostedList(const Particle& bp,vector<Hep3Vector>&sigtrk,vector<Hep3Vector>&otherb){
  // make lists of track/gamma in signal_b (all charged tracks as pi)
  vector<Particle> pi,pi_m,k_p,k_m,allgams;
  makeKPi(k_p,k_m,pi,pi_m);

  Mdst_gamma_Manager &gamma_mag = Mdst_gamma_Manager::get_manager();
  for(std::vector<Mdst_gamma>::iterator i = gamma_mag.begin(); i != gamma_mag.end(); ++i){
    allgams.push_back(Particle(*i));
  }

  vector<Particle*> final = bp.relation().finalStateParticles();
  pi.insert(pi.end(),pi_m.begin(),pi_m.end());
  pi.insert(pi.end(),allgams.begin(),allgams.end());

  for(int i=0; i<bp.relation().nFinalStateParticles(); i++){
    removeParticle(pi,bp.relation().finalStateParticle(i));
  }

  const Hep3Vector boostv = -BeamEnergy::CMBoost();
  sigtrk.clear(); otherb.clear();

  for(int i=0; i<final.size(); i++){
    HepLorentzVector p = final[i]->p();
    p.boost(boostv);
    sigtrk.push_back(p.vect());
  }

  for(int i=0; i<pi.size(); i++){
    HepLorentzVector p = pi[i].p();
    p.boost(boostv);
    otherb.push_back(p.vect());
  }
  return 0;
}

int b2dds::Mbc_deltaE(const Particle& bp,double& mbc,double& de){
  HepLorentzVector lv = bp.p();
  lv.boost(-BeamEnergy::CMBoost());
  const double benergy = BeamEnergy::E_beam_corr();
  de = lv.t() - benergy;
  lv.setT(benergy);
  mbc = lv.mag();
  return 0; 
}

int b2dds::TwoBodyMassFit(Particle& p,double& chisq){
  kmassfitter km;
  km.invariantMass(p.pType().mass());
  addTrack2fit(km,p.child(0));
  addTrack2fit(km,p.child(1));
  const int status = km.fit();

  if(status){
    chisq = -1;
    return status;
  } else{
    chisq = km.chisq()/km.dgf();
  }
  makeMother(km,p);
  return status;
}

double b2dds::CosHelAngle(const Particle& p){
  const Hep3Vector bv = -p.p().boostVector();
  HepLorentzVector ch1_lv = p.child(0).p();
  ch1_lv.boost(bv);
  return p.p().vect().unit().dot(ch1_lv.vect().unit());
}

#if defined(BELLE_NAMESPACE)
}
#endif

