double fold_em_q2_min(const genie::Interaction &state) {
#if XSEC_SCAN_HAS_KPHASESPACECUTS
  return genie::KPhaseSpaceCuts::Instance()->Q2MinCut(&state,0.);
#else
  (void)state;
  return genie::KPhaseSpace::GetQ2MinEM();
#endif
}

// Conditional state normalization shared by the event-matched QE/MEC folds.
// Caches are process-local, immutable after construction, and primed before fork.
const std::vector<WeightedValue> &fold_gauss(int n) {
  static std::map<int,std::vector<WeightedValue>> cache;
  auto found=cache.find(n); if(found!=cache.end()) return found->second;
  auto *t=gsl_integration_glfixed_table_alloc(n);
  std::vector<WeightedValue> nodes;
  for(int j=0;j<n;++j) { double x,w;gsl_integration_glfixed_point(0,1,j,&x,&w,t);nodes.push_back({x,w}); }
  gsl_integration_glfixed_table_free(t);
  return cache.emplace(n,std::move(nodes)).first->second;
}

double fold_cubic(const double *v,double x,int n) {
  int start=std::max(0,std::min(n-4,int(std::floor(x))-1));
  double y=0;
  for(int k=0;k<4;++k) { double w=1;for(int l=0;l<4;++l)if(l!=k)w*=(x-start-l)/double(k-l);y+=w*v[start+k]; }
  return y;
}

std::string fold_norm_key(const genie::XSecAlgorithmI &model,
                         const genie::Interaction &base,const Options &opt) {
  std::ostringstream s;s<<std::setprecision(17)<<&model<<'|'<<model.Id().Key()<<'|'
    <<base.AsString()<<'|'<<base.InitState().ProbeE(genie::kRfLab)<<'|'
    <<fold_em_q2_min(base)<<'|'<<opt.fold_norm_rel_tol;
  return s.str();
}

double qel_free_state_integral(const genie::XSecAlgorithmI &model,
                              genie::Interaction state,int order=32) {
  state.SetBit(genie::kISkipProcessChk);state.SetBit(genie::kISkipKinematicChk);
  state.SetBit(genie::kIAssumeFreeNucleon);
  const auto lim=state.PhaseSpace().Limits(genie::kKVQ2);
  if(!(lim.max>lim.min && lim.min>0)) return 0;
  const double l=std::log(lim.max/lim.min);double sum=0;
  for(auto z:fold_gauss(order)) {double q=lim.min*std::exp(l*z.value);state.KinePtr()->SetQ2(q);sum+=z.weight*q*l*model.XSec(&state,genie::kPSQ2fE);}
  return sum;
}

// Cache the smooth mean density, not a vanishing state integral. Extend the
// mean across an empty phase-space corner only for interpolation; consumers
// always supply the exact positive native width after their phase-space veto.
double qel_state_mean_density(const genie::XSecAlgorithmI &model,
                              genie::Interaction state,int order) {
  auto lim=state.PhaseSpace().Limits(genie::kKVQ2);
  if(lim.max>lim.min && lim.min>0)
    return qel_free_state_integral(model,state,order)/(lim.max-lim.min);
  state.SetBit(genie::kISkipProcessChk);state.SetBit(genie::kISkipKinematicChk);
  state.SetBit(genie::kIAssumeFreeNucleon);
  const auto &init=state.InitState();
  auto uncut=genie::utils::kinematics::electromagnetic::InelQ2Lim_W(
      init.ProbeE(genie::kRfHitNucRest),final_lepton_mass(state),
      init.Tgt().HitNucP4().M(),particle_mass(state.RecoilNucleonPdg()),0.);
  double lo=fold_em_q2_min(state),hi=uncut.max,sum=0;
  if(!(hi>0 && lo>0))throw std::runtime_error("QE normalization extension outside supported phase space");
  for(auto z:fold_gauss(order)) {
    state.KinePtr()->SetQ2(lo+(hi-lo)*z.value);
    double v=model.XSec(&state,genie::kPSQ2fE);
    if(!(v>0) || !std::isfinite(v))throw std::runtime_error("Invalid QE normalization extension");
    sum+=z.weight*v;
  }
  return sum;
}

struct QELStateNormalization {
  int n=0;double pmax=0,acceptance=0,acceptance_error=0,total=0,interpolation_error=0;
  long samples=0;
  std::vector<double> log_sigma;
  double interpolate(const std::vector<double>&table,double p,double c) const {
    if(p<0 || p>pmax*(1+1e-10) || std::abs(c)>1+1e-10)
      throw std::runtime_error("QE state outside normalization table");
    double xp=std::min(double(n-1),p/pmax*(n-1)),xc=(clamp(c,-1.,1.)+1)*.5*(n-1);
    int start=std::max(0,std::min(n-4,int(std::floor(xp))-1));
    double value=0;
    for(int k=0;k<4;++k) {double w=1;for(int l=0;l<4;++l)if(l!=k)w*=(xp-start-l)/double(k-l);value+=w*fold_cubic(&table[(start+k)*n],xc,n);}
    return value;
  }
  double sigma(double p,double c,double width) const { return std::exp(interpolate(log_sigma,p,c))*width; }
};

const QELStateNormalization &qel_state_normalization(
    const genie::XSecAlgorithmI &model,const genie::Interaction &base,
    const QELFoldConfig &cfg,const QELRadialMeasure &m,const Options &opt) {
  static std::map<std::string,QELStateNormalization> cache;
  const auto key=fold_norm_key(model,base,opt);auto it=cache.find(key);if(it!=cache.end())return it->second;
  QELStateNormalization norm;
  for(int j=0;j<m.bins;++j)if(m.suffix[0][j]>0)norm.pmax=(j+1)*m.dp;
  const double removal=qel_fold_removal_energy(base,opt);
  auto state_at=[&](double p,double c) {genie::Interaction state(base);qel_bind_old_fermi_mover(&state,cfg,TVector3(p*std::sqrt(std::max(0.,1-c*c)),0,p*c),removal);return state;};
  for(norm.n=17;norm.n<=257;norm.n=2*norm.n-1) {
    norm.log_sigma.clear();norm.interpolation_error=0;
    for(int i=0;i<norm.n;++i)for(int j=0;j<norm.n;++j) {
      auto state=state_at(norm.pmax*i/(norm.n-1),2.*j/(norm.n-1)-1);
      double mean=qel_state_mean_density(model,state,32);
      if(!(mean>0) || !std::isfinite(mean))throw std::runtime_error("Invalid QE state mean density");
      norm.log_sigma.push_back(std::log(mean));
    }
    for(int i=0;i<norm.n-1;++i)for(int j=0;j<norm.n-1;++j) {
      double p=norm.pmax*(i+.5)/(norm.n-1),c=2.*(j+.5)/(norm.n-1)-1;
      auto state=state_at(p,c);double exact=qel_state_mean_density(model,state,64);
      double error=std::abs(norm.sigma(p,c,1.)/exact-1);
      norm.interpolation_error=std::max(norm.interpolation_error,error);
    }
    std::cerr<<"qe-norm-table hit="<<base.InitState().Tgt().HitNucPdg()<<" nodes="<<norm.n<<" error="<<norm.interpolation_error<<"\n";
    if(norm.interpolation_error<opt.fold_norm_rel_tol*.1)break;
  }
  if(norm.n>257)throw std::runtime_error("QE normalization table did not converge");
  // Sample the native joint radial/momentum measure. Average the Pauli veto
  // analytically over the final lepton's rest-frame azimuth, reducing variance.
  const auto &t=base.InitState().Tgt();const int nr=8192;
  std::vector<double> rcdf(nr+1,0);
  for(int i=1;i<=nr;++i)rcdf[i]=rcdf[i-1]+qel_radial_mass(t.A(),m.rmax*(i-1)/nr,m.rmax*i/nr)/m.norm;
  for(double &v:rcdf)v/=rcdf.back();
  std::vector<std::vector<double>> pcdf(m.pdf.size(),std::vector<double>(m.bins+1,0));
  for(size_t s=0;s<m.pdf.size();++s)for(int j=0;j<m.bins;++j)pcdf[s][j+1]=pcdf[s][j]+m.pdf[s][j]*m.dp;
  auto invert=[](const std::vector<double>&cdf,double u) {size_t j=std::min(cdf.size()-2,size_t(std::upper_bound(cdf.begin(),cdf.end(),u)-cdf.begin()-1));return j+(u-cdf[j])/(cdf[j+1]-cdf[j]);};
  const auto *lfg=dynamic_cast<const genie::LocalFGM*>(cfg.nucl_model);
  const double E=base.InitState().ProbeE(genie::kRfLab),mf=particle_mass(base.RecoilNucleonPdg()),ml=final_lepton_mass(base);
  std::unique_ptr<TLorentzVector> beam(base.InitState().GetProbeP4(genie::kRfLab));
  constexpr int replicas=8;double sums[replicas]={},shift[replicas][4];
  std::uint64_t rng=0x51454e4f524dULL;
  for(auto &row:shift)for(double &v:row){rng+=0x9e3779b97f4a7c15ULL;auto z=rng;z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;z=(z^(z>>27))*0x94d049bb133111ebULL;z^=z>>31;v=(z>>11)*0x1.0p-53;}
  long previous_n=0;double previous=-1;
  for(long n=16384;n<=1048576;n*=2) {
    for(int rep=0;rep<replicas;++rep)for(long i=previous_n;i<n;++i) {
      auto u=[&](int d){double v=qel_lattice_u(i,d)+shift[rep][d];return v-std::floor(v);};
      double r=invert(rcdf,u(0))*m.rmax/nr;
      size_t shell=std::min(m.pdf.size()-1,size_t(std::upper_bound(m.edges.begin(),m.edges.end(),r)-m.edges.begin()-1));
      if(pcdf[shell].back()==0)continue;
      double p=invert(pcdf[shell],u(1))*m.dp,c=2*u(2)-1;
      auto state=state_at(p,c);state.SetBit(genie::kISkipProcessChk);state.SetBit(genie::kISkipKinematicChk);state.SetBit(genie::kIAssumeFreeNucleon);
      auto lim=state.PhaseSpace().Limits(genie::kKVQ2);if(lim.max<=lim.min)continue;
      double l=std::log(lim.max/lim.min),q=lim.min*std::exp(l*u(3));state.KinePtr()->SetQ2(q);
      const auto &hit=state.InitState().Tgt().HitNucP4();TLorentzVector br(*beam);br.Boost(-hit.BoostVector());double er=br.E(),x,y;
      genie::utils::kinematics::WQ2toXY(er,hit.M(),mf,q,x,y);
      double el=(1-y)*er,pll=el-(q+ml*ml)/(2*er),pt2=el*el-pll*pll-ml*ml;if(pt2<0)continue;
      auto beta=hit.BoostVector();double bpar=beta.Dot(br.Vect().Unit()),gamma=hit.E()/hit.M();
      double a=gamma*(el+bpar*pll),b=gamma*std::sqrt(std::max(0.,beta.Mag2()-bpar*bpar))*std::sqrt(pt2);
      double kf=lfg->LocalFermiMomentum(t,base.RecoilNucleonPdg(),r),cut=hit.E()+E-a-std::sqrt(mf*mf+kf*kf);
      double pass=b>1e-15?1-std::acos(clamp(cut/b,-1.,1.))/kPi:(cut>=0?1.:0.);
      sums[rep]+=q*l*model.XSec(&state,genie::kPSQ2fE)/norm.sigma(p,c,lim.max-lim.min)*pass;
    }
    double mean=0;for(double s:sums)mean+=s/n/replicas;double variance=0;for(double s:sums)variance+=std::pow(s/n-mean,2)/(replicas*(replicas-1));
    norm.acceptance=mean;norm.acceptance_error=std::max(3*std::sqrt(variance),previous<0?1.:std::abs(mean-previous));norm.samples=n*replicas;
    std::cerr<<"qe-event-norm-progress target="<<t.Pdg()<<" hit="<<t.HitNucPdg()<<" E="<<E<<" samples="<<norm.samples<<" acceptance="<<mean<<" error="<<norm.acceptance_error<<"\n";
    if(previous>=0 && norm.acceptance_error<opt.fold_norm_rel_tol*mean)break;
    previous=mean;previous_n=n;
  }
  if(!(norm.acceptance>0) || norm.acceptance_error>opt.fold_norm_rel_tol*norm.acceptance)
    throw std::runtime_error("QE event normalization did not converge");
  genie::Interaction selected(base);norm.total=model.Integral(&selected);
  if(!(norm.total>0) || !std::isfinite(norm.total))throw std::runtime_error("Invalid native QE channel integral");
  std::cerr<<"qe-event-norm target="<<t.Pdg()<<" hit="<<t.HitNucPdg()<<" E="<<E<<" table="<<norm.n<<" interpolation_error="<<norm.interpolation_error<<" acceptance="<<norm.acceptance<<" relative_error="<<norm.acceptance_error/norm.acceptance<<" channel_total="<<norm.total<<"\n";
  return cache.emplace(key,std::move(norm)).first->second;
}

// Empirical MEC samples W,Q2 conditionally on the moving on-shell cluster.
// Its Integral() is a separate QE-derived channel rate, not this shape integral.
double mec_state_shape_integral(const genie::XSecAlgorithmI &model,
                               genie::Interaction state,double E,int order=32) {
  const double M=particle_mass(state.InitState().Tgt().HitNucPdg()),ml=final_lepton_mass(state);
  state.InitStatePtr()->TgtPtr()->SetHitNucP4(TLorentzVector(0,0,0,M));
  state.InitStatePtr()->SetProbeP4(TLorentzVector(0,0,E,E));
  auto wl=genie::utils::kinematics::electromagnetic::InelWLim(E,ml,M);
  double lo=std::max(1.88,wl.min),hi=std::min(3.,wl.max),qmin=std::max(.01,fold_em_q2_min(state));
  auto qlim=[&](double w){return genie::utils::kinematics::electromagnetic::InelQ2Lim_W(E,ml,M,w,qmin);};
  if(!(hi>lo) || qlim(lo).max<=qmin)return 0;
  if(qlim(hi).max<=qmin) {double a=lo,b=hi;for(int j=0;j<60;++j){double mid=(a+b)/2;if(qlim(mid).max>qmin)a=mid;else b=mid;}hi=(a+b)/2;}
  double sum=0;
  for(auto w:fold_gauss(order)) {double W=lo+(hi-lo)*w.value;auto lim=qlim(W);double qmax=std::min(8.,lim.max);if(qmax<=qmin)continue;double l=std::log(qmax/qmin);state.KinePtr()->SetW(W);
    for(auto z:fold_gauss(order)){double q=qmin*std::exp(l*z.value);state.KinePtr()->SetQ2(q);sum+=w.weight*(hi-lo)*z.weight*q*l*model.XSec(&state,genie::kPSWQ2fE);}
  }
  return sum;
}
struct MECStateNormalization {
  int n=0;double emin=0,emax=0,total=0,interpolation_error=0;std::vector<double> log_sigma;
  double sigma(double E)const {return std::exp(fold_cubic(log_sigma.data(),(E-emin)/(emax-emin)*(n-1),n));}
};
const MECStateNormalization &mec_state_normalization(
    const genie::XSecAlgorithmI &model,const genie::Interaction &base,const Options &opt) {
  static std::map<std::string,MECStateNormalization> cache;auto key=fold_norm_key(model,base,opt);auto it=cache.find(key);if(it!=cache.end())return it->second;
  if(!base.ProcInfo().IsEM())throw std::runtime_error("Event-normalized empirical MEC currently requires EM");
  double E=base.InitState().ProbeE(genie::kRfLab);MECStateNormalization norm;norm.emin=.65*E;norm.emax=1.55*E;
  for(norm.n=17;norm.n<=513;norm.n=norm.n*2-1) {
    norm.log_sigma.clear();norm.interpolation_error=0;
    for(int j=0;j<norm.n;++j){double energy=norm.emin+(norm.emax-norm.emin)*j/(norm.n-1),s=mec_state_shape_integral(model,base,energy);if(!(s>0))throw std::runtime_error("MEC normalization table has empty phase space");norm.log_sigma.push_back(std::log(s));}
    for(int j=0;j<norm.n-1;++j){double energy=norm.emin+(norm.emax-norm.emin)*(j+.5)/(norm.n-1),exact=mec_state_shape_integral(model,base,energy,64);norm.interpolation_error=std::max(norm.interpolation_error,std::abs(norm.sigma(energy)/exact-1));}
    if(norm.interpolation_error<opt.fold_norm_rel_tol*.1)break;
  }
  if(norm.n>513)throw std::runtime_error("MEC normalization table did not converge");
  genie::Interaction selected(base);norm.total=model.Integral(&selected);
  if(!(norm.total>=0) || !std::isfinite(norm.total))throw std::runtime_error("Invalid native MEC channel integral");
  std::cerr<<"mec-event-norm target="<<base.InitState().Tgt().Pdg()<<" hit="<<base.InitState().Tgt().HitNucPdg()<<" E="<<E<<" table="<<norm.n<<" interpolation_error="<<norm.interpolation_error<<" channel_total="<<norm.total<<"\n";
  return cache.emplace(key,std::move(norm)).first->second;
}
