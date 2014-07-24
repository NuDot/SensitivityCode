// Note on mass:
// Mass enters in 3 places: fRhoIsotope_p_kg, fBG_c_p_keV_kg_y, and
// fDeployMasses_kg. "kg" can refer to either detector mass or to mass of
// isotope. As long as all are three are in the same units the
// calculations will come out correct.

#include <iostream>
#include "Math/ProbFuncMathCore.h"
#include "TAxis.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"

#include <float.h> // this is to include DBL_MAX, DBL_MIN, etc
#include <limits.h> // this is to include LONG_MAX, LONG_MIN, etc
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Units/PhysicalConstants.h"

using namespace std;
using namespace CLHEP;

// Local additional units and constants - add to CLHEP namespace
namespace CLHEP
{
  // new units
  static const double gigahertz = 1.e+9*hertz;
  static const double minute = 60.0*s;
  static const double hour = 60.0*minute;
  static const double day = 24.0*hour;
  static const double year = 365.2425*day;
  static const double tonne = 1000.0*kg;
  static const double meV = 1.e-3*eV;

  // aliases
  static const double us = microsecond;

  static const double Hz = hertz;
  static const double kHz = kilohertz;
  static const double MHz = megahertz;
  static const double GHz = gigahertz;
}


class SensClass {
public:
  enum EStatMethod { kFC90CLSens, k5SigmaDL };
  
  enum EIsotope { kGe76, kTe130, kXe136, kNd150 };
  static double GetNaturalAbundance(EIsotope isotope);
  virtual double GetNaturalAbundance() const { return GetNaturalAbundance(fIsotope); }
  static double GetAtomicMass(EIsotope isotope);
  virtual double GetAtomicMass() const { return GetAtomicMass(fIsotope); }
  
  enum ENMEMethod { kQRPA, kShellModel };
  static double GetNME(EIsotope isotope, ENMEMethod nmeMethod, double ga);
  virtual double GetNME() const { return GetNME(fIsotope, fNMEMethod, fGa); }
  static double GetQRPANME(EIsotope isotope, double ga);
  static double GetShellModelNME(EIsotope isotope, double ga);
  static double GetPhaseSpaceFactor(EIsotope isotope);
  virtual double GetPhaseSpaceFactor() const { return GetPhaseSpaceFactor(fIsotope); }
  static double GetMbb(double tHalf, EIsotope isotope, ENMEMethod nmeMethod, double ga);
  static double GetTHalf(double mbb, EIsotope isotope, ENMEMethod nmeMethod, double ga);
  
  enum ETHalfUnit { kYears, k1e25Years, k1e26Years, k1e27Years };
  static double GetTHalfUnitValue(ETHalfUnit unit);
  static const char* GetTHalfUnitName(ETHalfUnit unit);
  
  enum EExposureUnit { kKgDays, kKgYears, kTonneYears };
  static double GetExposureUnitValue(EExposureUnit unit);
  static const char* GetExposureUnitName(EExposureUnit unit);
  
  enum EMbbUnit { k_eV, kmeV };
  static double GetMbbUnitValue(EMbbUnit unit);
  static const char* GetMbbUnitName(EMbbUnit unit);
  
  enum ETimeUnit { kDays, kCalendarYears };
  static double GetTimeUnitValue(ETimeUnit unit);
  static const char* GetTimeUnitName(ETimeUnit unit);
  
  // FHWM = 2.*sqrt(2.*log(2))*sigma
  static inline double GetFWHMEfficiency() { return TMath::Erf(sqrt(log(2.))); }
  
  // for bg-limited, optimal ROI is ~2.8 x sigma
  static inline double GetROIEfficiency() { return TMath::Erf(2.8/2./sqrt(2.)); }
  
public:
  SensClass(EIsotope isotope = kGe76, double enrFrac = 0., 
	    double efficiency = GetROIEfficiency(), double backgroundRate = 0., 
	    EStatMethod statMethod = kFC90CLSens, ENMEMethod nmeMethod = kQRPA,
	    double ga = 1.27);
  virtual ~SensClass() {}
  
  virtual inline void SetIsotope(EIsotope isotope) { fIsotope = isotope; }
  virtual inline void SetEnrichmentFraction(double enrFrac) { fEnrichmentFraction = enrFrac; }
  virtual inline void SetEfficiency(double efficiency) { fEfficiency = efficiency; }
  virtual inline void SetBackgroundRate(double backgroundRate) { fBackgroundRate = backgroundRate; }
  virtual inline void SetStatMethod(EStatMethod statMethod) { fStatMethod = statMethod; }
  virtual inline void SetNMEMethod(ENMEMethod nmeMethod) { fNMEMethod = nmeMethod; }
  virtual void AddDeployment(double mass, double tStart, double tStop = 0);
  virtual double GetExposure(double time) const;
  virtual inline void CombineWithMe(const SensClass& sensClass) { fExpsToCombine.push_back(&sensClass); }
  
  virtual double GetTHalfVsExposure(double exposure) const;
  double GetTHalfVsExposureForTF1(double* exposure, double* parameters);
  virtual TF1& GetTHalfVsExposureFunction(bool forceSetUnits = true, ETHalfUnit tHalfUnit = k1e26Years, EExposureUnit exposureUnit = kTonneYears);
  
  virtual double GetMbbVsExposure(double exposure) const;
  double GetMbbVsExposureForTF1(double* exposure, double* parameters);
  virtual TF1& GetMbbVsExposureFunction(bool forceSetUnits = true, EMbbUnit mbbUnit = kmeV, EExposureUnit exposureUnit = kTonneYears);

  virtual double GetTHalfVsTime(double time, double liveFraction = 1.0) const;
  double GetTHalfVsTimeForTF1(double* time, double* parameters);
  virtual TF1& GetTHalfVsTimeFunction(double liveFraction = 1.0, bool forceSetUnits = true, ETHalfUnit tHalfUnit = k1e26Years, 
                                        ETimeUnit timeUnit = kCalendarYears, double tOffset = 0);

  virtual double GetMbbVsTime(double time, double liveFraction = 1.0) const;
  double GetMbbVsTimeForTF1(double* time, double* parameters);
  virtual TF1& GetMbbVsTimeFunction(double liveFraction = 1.0, bool forceSetUnits = true, EMbbUnit mbbUnit = kmeV, 
				    ETimeUnit timeUnit = kCalendarYears, double tOffset = 0);

  virtual double GetStatistic(double backgroundCounts) const;
  virtual double GetFC90CLSensitivity(double backgroundCounts) const;
  virtual double Get5SigmaDiscoveryLimit(double backgroundCounts) const;
  virtual inline void SetBeta(double beta);
  
protected:
  EIsotope fIsotope;
  double fEnrichmentFraction;
  double fEfficiency;
  double fBackgroundRate; // counts per unit exposure
  EStatMethod fStatMethod;
  double fBeta; // for error of second kind (see, e.g., 5sig disc. limit)
  ENMEMethod fNMEMethod;
  double fGa;
  vector<double> fDeployMasses;  // the deployed mass
  vector<double> fDeployTimes;   // the times at which masses are deployed
  vector<double> fDeployStopTimes; // the times at which delpoyments are removed
  vector<const SensClass*> fExpsToCombine;
  TF1* fTHalfVsExpFunc;
  TF1* fMbbVsExpFunc;
  TF1* fTHalfVsTimeFunc;
  TF1* fMbbVsTimeFunc;


};

SensClass::SensClass(EIsotope isotope, double enrFrac, double efficiency, double backgroundRate, 
                     EStatMethod statMethod, ENMEMethod nmeMethod, double ga) : 
fIsotope(isotope), fEnrichmentFraction(enrFrac), fEfficiency(efficiency), fBackgroundRate(backgroundRate), 
fStatMethod(statMethod), fBeta(0.9), fNMEMethod(nmeMethod), fGa(ga),
fTHalfVsExpFunc(NULL), fMbbVsExpFunc(NULL), fTHalfVsTimeFunc(NULL), fMbbVsTimeFunc(NULL)
{
  if(fEnrichmentFraction == 0.0) {
    fEnrichmentFraction = GetNaturalAbundance(fIsotope);
  }
}

void SensClass::AddDeployment(double mass, double tStart, double tStop)
{
  fDeployMasses.push_back(mass);
  fDeployTimes.push_back(tStart);
  fDeployStopTimes.push_back(tStop);
}

double SensClass::GetExposure(double time) const
{
  double exposure = 1.e-99;
  for(size_t i=0; i<fDeployMasses.size(); i++) {
    if(time > fDeployTimes[i]) {
      if(fDeployStopTimes[i] > 0 && time >= fDeployStopTimes[i]) {
        exposure += fDeployMasses[i]*(fDeployStopTimes[i]-fDeployTimes[i]);
      }
      else {
        exposure += fDeployMasses[i]*(time-fDeployTimes[i]);
      }
    }
  }
  return exposure;
}

double SensClass::GetTHalfVsExposure(double exposure) const
{
  return TMath::Log(2)*exposure*fEfficiency*fEnrichmentFraction
         / (GetAtomicMass(fIsotope)*GetStatistic(fBackgroundRate*exposure));
}

double SensClass::GetTHalfVsExposureForTF1(double* exposure, double* parameters)
{
  double tHalfUnit = parameters[0];
  double exposureUnit = parameters[1];
  return GetTHalfVsExposure((*exposure)*exposureUnit)/tHalfUnit;
}

TF1& SensClass::GetTHalfVsExposureFunction(bool forceSetUnits, ETHalfUnit tHalfUnit, EExposureUnit exposureUnit)
{
  if(fTHalfVsExpFunc == NULL) {
    int i = 0;
    while(fTHalfVsExpFunc == NULL) {
      TString name = TString::Format("fTHalfVsExpFunc_%d", i++);
      if(gROOT->FindObject(name) == NULL) {
        fTHalfVsExpFunc = new TF1(name, this, &SensClass::GetTHalfVsExposureForTF1, 0, 10., 2, "SensClass", "GetTHalfVsExposureForTF1");
        forceSetUnits = true;
      }
    }
  }
  if(forceSetUnits) {
    fTHalfVsExpFunc->SetParameters(GetTHalfUnitValue(tHalfUnit), GetExposureUnitValue(exposureUnit));
    fTHalfVsExpFunc->GetXaxis()->SetTitle(TString::Format("exposure [%s]", GetExposureUnitName(exposureUnit)));
    fTHalfVsExpFunc->GetYaxis()->SetTitle(TString::Format("T_{1/2} [%s]", GetTHalfUnitName(tHalfUnit)));
  }
  return *fTHalfVsExpFunc;
}

double SensClass::GetMbbVsExposure(double exposure) const
{
  return GetMbb(GetTHalfVsExposure(exposure), fIsotope, fNMEMethod, fGa);
}

double SensClass::GetMbbVsExposureForTF1(double* exposure, double* parameters)
{
  double mbbUnit = parameters[0];
  double exposureUnit = parameters[1];
  return GetMbbVsExposure((*exposure)*exposureUnit)/mbbUnit;
}

TF1& SensClass::GetMbbVsExposureFunction(bool forceSetUnits, EMbbUnit mbbUnit, EExposureUnit exposureUnit)
{
  if(fMbbVsExpFunc == NULL) {
    int i = 0;
    while(fMbbVsExpFunc == NULL) {
      TString name = TString::Format("fMbbVsExpFunc_%d", i++);
      if(gROOT->FindObject(name) == NULL) {
        fMbbVsExpFunc = new TF1(name, this, &SensClass::GetMbbVsExposureForTF1, 0, 10., 2, "SensClass", "GetMbbVsExposureForTF1");
        forceSetUnits = true;
      }
    }
  }
  if(forceSetUnits) {
    fMbbVsExpFunc->SetParameters(GetMbbUnitValue(mbbUnit), GetExposureUnitValue(exposureUnit));
    fMbbVsExpFunc->GetXaxis()->SetTitle(TString::Format("exposure [%s]", GetExposureUnitName(exposureUnit)));
    fMbbVsExpFunc->GetYaxis()->SetTitle(TString::Format("m_{#beta#beta} [%s]", GetMbbUnitName(mbbUnit)));
  }
  return *fMbbVsExpFunc;
}

double SensClass::GetTHalfVsTime(double time, double liveFraction) const
{
  double tHalf = GetTHalfVsExposure(GetExposure(time)*liveFraction);
  if(fExpsToCombine.size() == 0) return tHalf;
  // only valid for FC sens U.L. on the decay RATE so far
  // assumes combining ULs is like adding 0 +/- UL in quadrature
  double weight = tHalf*tHalf;
  for(size_t i=0; i<fExpsToCombine.size(); i++) {
    tHalf = fExpsToCombine[i]->GetTHalfVsTime(time, liveFraction);
    weight += tHalf*tHalf;
  }
  return sqrt(weight);
}

double SensClass::GetTHalfVsTimeForTF1(double* time, double* parameters)
{
  double tHalfUnit = parameters[0];
  double timeUnit = parameters[1];
  double tOffset = parameters[2];
  double liveFraction = parameters[3];
  return GetTHalfVsTime((*time)*timeUnit+tOffset, liveFraction)/tHalfUnit;
}

TF1& SensClass::GetTHalfVsTimeFunction(double liveFraction, bool forceSetUnits, ETHalfUnit tHalfUnit, ETimeUnit timeUnit, double tOffset)
{
  if(fTHalfVsTimeFunc == NULL) {
    int i = 0;
    while(fTHalfVsTimeFunc == NULL) {
      TString name = TString::Format("fTHalfVsTimeFunc_%d", i++);
      if(gROOT->FindObject(name) == NULL) {
        double firstTime = 1.e99;
        for(size_t j=0; j<fDeployTimes.size(); j++) {
          if(fDeployTimes[j] < firstTime) firstTime = fDeployTimes[j];
        }
        if(firstTime == 1.e99) firstTime = tOffset;
        double tUnitVal = GetTimeUnitValue(timeUnit);
        double x0 = (firstTime - tOffset)/tUnitVal;
        fTHalfVsTimeFunc = new TF1(name, this, &SensClass::GetTHalfVsTimeForTF1, x0, x0 + 10., 4, "SensClass", "GetTHalfVsTimeForTF1");
        forceSetUnits = true;
      }
    }
  }
  if(forceSetUnits) {
    fTHalfVsTimeFunc->SetParameters(GetTHalfUnitValue(tHalfUnit), GetTimeUnitValue(timeUnit), tOffset, liveFraction);
    fTHalfVsTimeFunc->GetXaxis()->SetTitle(TString::Format("%s", GetTimeUnitName(timeUnit)));
    fTHalfVsTimeFunc->GetYaxis()->SetTitle(TString::Format("T_{1/2} [%s]", GetTHalfUnitName(tHalfUnit)));
  }
  return *fTHalfVsTimeFunc;
}

double SensClass::GetMbbVsTime(double time, double liveFraction) const
{
  //return GetMbbVsExposure(GetExposure(time));
  // need to call GetTHalfVsTime(time) where experiments can be combined
  return GetMbb(GetTHalfVsTime(time, liveFraction), fIsotope, fNMEMethod, fGa);
}

double SensClass::GetMbbVsTimeForTF1(double* time, double* parameters)
{
  double mbbUnit = parameters[0];
  double timeUnit = parameters[1];
  double tOffset = parameters[2];
  double liveFraction = parameters[3];
  return GetMbbVsTime((*time)*timeUnit+tOffset, liveFraction)/mbbUnit;
}

TF1& SensClass::GetMbbVsTimeFunction(double liveFraction, bool forceSetUnits, EMbbUnit mbbUnit, ETimeUnit timeUnit, double tOffset)
{
  if(fMbbVsTimeFunc == NULL) {
    int i = 0;
    while(fMbbVsTimeFunc == NULL) {
      TString name = TString::Format("fMbbVsTimeFunc_%d", i++);
      if(gROOT->FindObject(name) == NULL) {
        double firstTime = 1.e99;
        for(size_t j=0; j<fDeployTimes.size(); j++) {
          if(fDeployTimes[j] < firstTime) firstTime = fDeployTimes[j];
        }
        if(firstTime == 1.e99) firstTime = tOffset;
        double tUnitVal = GetTimeUnitValue(timeUnit);
        double x0 = (firstTime - tOffset)/tUnitVal;
        fMbbVsTimeFunc = new TF1(name, this, &SensClass::GetMbbVsTimeForTF1, x0, x0 + 10., 4, "SensClass", "GetMbbVsTimeForTF1");
        forceSetUnits = true;
      }
    }
  }
  if(forceSetUnits) {
    fMbbVsTimeFunc->SetParameters(GetMbbUnitValue(mbbUnit), GetTimeUnitValue(timeUnit), tOffset, liveFraction);
    fMbbVsTimeFunc->GetXaxis()->SetTitle(TString::Format("%s", GetTimeUnitName(timeUnit)));
    fMbbVsTimeFunc->GetYaxis()->SetTitle(TString::Format("m_{#beta#beta} [%s]", GetMbbUnitName(mbbUnit)));
    if(fMbbVsTimeFunc->GetMaximum() > 1.*eV/GetMbbUnitValue(mbbUnit)) {
      fMbbVsTimeFunc->SetMaximum(1.*eV/GetMbbUnitValue(mbbUnit));
    }
  }
  return *fMbbVsTimeFunc;
}


double SensClass::GetStatistic(double backgroundCounts) const
{
  if(fStatMethod == kFC90CLSens) return GetFC90CLSensitivity(backgroundCounts);
  if(fStatMethod == k5SigmaDL)   return Get5SigmaDiscoveryLimit(backgroundCounts);
  cout << "Unknown fStatMethod " << fStatMethod << endl;
  return 1.e99;
}

double SensClass::GetFC90CLSensitivity(double backgroundCounts) const
{
  if(backgroundCounts > 300) return 1.21 + 1.73*sqrt(backgroundCounts);
  static TGraph* gFCSensitivity = NULL;
  if(gFCSensitivity == NULL) {
    TFile* interpFile = TFile::Open("FCSensitivityGraph.root");
    gFCSensitivity = (TGraph*) interpFile->Get("gFCSensitivity");
    if(gFCSensitivity == NULL) {
      cout << "Error: coudn't get gFCSensitivity from FCSensitivityGraph.root" << endl;
      return 1e99;
    }
  }
  if(backgroundCounts < 0) backgroundCounts = 0;
  return gFCSensitivity->Eval(backgroundCounts);
}

double SensClass::Get5SigmaDiscoveryLimit(double backgroundCounts) const
{
  static double fiveSigmaCL = TMath::Erfc(5./sqrt(2.));

  static TF1 fUpperGamma("fUpperGamma", "TMath::Gamma(x+1, [0])");
  fUpperGamma.SetParameter(0, backgroundCounts);
  double counts5Sig = fUpperGamma.GetX(fiveSigmaCL, backgroundCounts, TMath::Max(20., backgroundCounts + 10.*sqrt(backgroundCounts)));

  static TF1 fUpperGammaForMu("fUpperGamma", "TMath::Gamma([0], x)");
  fUpperGammaForMu.SetParameter(0, counts5Sig+1.);
  double betaSigma = TMath::ErfInverse(fBeta)*sqrt(2);
  return fUpperGammaForMu.GetX(fBeta, counts5Sig, TMath::Max(20., counts5Sig + 2.*betaSigma*sqrt(counts5Sig))) - backgroundCounts;

/*
  // Old way for integral counts5Sig:
  // find the counts that would represent 5 sigma above bg
  double counts5Sig = TMath::Floor(backgroundCounts + 5.*sqrt(backgroundCounts));
  while(counts5Sig > 0.0 && ROOT::Math::poisson_cdf_c(counts5Sig, backgroundCounts) < fiveSigmaCL) counts5Sig -= 1.0;
  while(ROOT::Math::poisson_cdf_c(counts5Sig, backgroundCounts) > fiveSigmaCL) counts5Sig += 1.0;
  counts5Sig += 1.0;

  // now find the signal such that s+b gives c > c5sig a fraction [beta] of the time
  // s+b - n_beta(sqrt(s+b)) = c: solve quadratic function. s+b => sb, n_beta => b
  // sb = [2c + b^2 +/- sqrt([2c + b^2]^2 - 4c^2)]/2. Want sb > c: choose +
  double nSigBeta = TMath::ErfInverse((fBeta-0.5)*2.)*sqrt(2.);
  double bQuad = counts5Sig + nSigBeta*nSigBeta/2.;
  double sigPlusBG = bQuad + sqrt(bQuad*bQuad - counts5Sig*counts5Sig);

  double delta = sigPlusBG;
  while(TMath::Abs(delta/sigPlusBG) > 1.e-10) {
    double f1 = TMath::PoissonI(counts5Sig, sigPlusBG);
    double f = (ROOT::Math::poisson_cdf(counts5Sig, sigPlusBG) - (1.-fBeta))/f1;
    double f2 = 1.-counts5Sig/sigPlusBG;
    delta = f*(1. + 0.5*f2*f);
    sigPlusBG += delta;
  }
  return sigPlusBG - backgroundCounts;
*/
}

void SensClass::SetBeta(double beta) 
{ 
  if(beta >= 1. || beta <= 0.) {
    cout << "Invalid beta " << beta << " (must satisfy 0 < beta < 1)" << endl;
    return;
  }
  fBeta = beta; 
}

double SensClass::GetTHalfUnitValue(ETHalfUnit unit)
{
  if(unit == kYears)     return       year;
  if(unit == k1e25Years) return 1.e25*year;
  if(unit == k1e26Years) return 1.e26*year;
  if(unit == k1e27Years) return 1.e27*year;
  return 1.;
}

const char* SensClass::GetTHalfUnitName(ETHalfUnit unit)
{
  if(unit == kYears)     return         "years";
  if(unit == k1e25Years) return "10^{25} years";
  if(unit == k1e26Years) return "10^{26} years";
  if(unit == k1e27Years) return "10^{27} years";
  return "arb.";
}

double SensClass::GetExposureUnitValue(EExposureUnit unit)
{
  if(unit == kKgDays)     return kg*day;
  if(unit == kKgYears)    return kg*year;
  if(unit == kTonneYears) return tonne*year;
  return 1.;
}

const char* SensClass::GetExposureUnitName(EExposureUnit unit)
{
  if(unit == kKgDays)     return "kg days";
  if(unit == kKgYears)    return "kg years";
  if(unit == kTonneYears) return "tonne years";
  return "arb.";
}

double SensClass::GetMbbUnitValue(EMbbUnit unit)
{
  if(unit == k_eV) return       eV;
  if(unit == kmeV) return 1.e-3*eV;
  return 1.;
}

const char* SensClass::GetMbbUnitName(EMbbUnit unit)
{
  if(unit == k_eV) return  "eV";
  if(unit == kmeV) return "meV";
  return "arb.";
}

double SensClass::GetTimeUnitValue(ETimeUnit unit)
{
  if(unit == kDays)          return  day;
  if(unit == kCalendarYears) return year;
  return 1.;
}

const char* SensClass::GetTimeUnitName(ETimeUnit unit)
{
  if(unit == kDays)          return "days";
  if(unit == kCalendarYears) return "year";
  return "arb.";
}

double SensClass::GetNaturalAbundance(EIsotope isotope)
{
  if(isotope == kGe76)  return 0.0761;
  if(isotope == kTe130) return 0.3408;
  if(isotope == kXe136) return 0.0887;
  if(isotope == kNd150) return 0.056;
  cout << "Unknown isotope " << isotope << endl;
  return 1.e-99;
}

double SensClass::GetAtomicMass(EIsotope isotope)
{
  if(isotope == kGe76)  return 75.9214016*amu;
  if(isotope == kTe130) return 129.906229*amu;
  if(isotope == kXe136) return 135.907214*amu;
  if(isotope == kNd150) return 149.920887*amu;
  cout << "Unknown isotope " << isotope << endl;
  return 1.e99;
}

double SensClass::GetQRPANME(EIsotope isotope, double ga)
{
  double ga2Fact = (ga == 1.25) ? 1 : ga*ga/(1.25*1.25);
  // Following Steve Elliott's NME_master: The QRPA values from Table III
  // in Ref. [10] using (R)QRPA with CCM SRC. The QRPA values for 150Nd and 160Gd
  // come from Ref. [22], which includes deformation.
  // [10] Simkovic F, Faessler A, uther H, Rodin V and Stauf M 2009 Phys. Rev. C 79 055501 (Preprint arXiv:0902.0331)
  // [22] Fang D L, Faessler A and Rodin V 2011 (Preprint arXiv:1101.2149)
  if(isotope == kGe76)  return (4.07+6.64)/2.*ga2Fact;
  if(isotope == kTe130) return (3.29+5.37)/2.*ga2Fact;
  if(isotope == kXe136) return (2.06+3.23)/2.*ga2Fact;
  if(isotope == kNd150) return (2.55+6.12)/2.*ga2Fact;
  cout << "Unknown isotope " << isotope << endl;
  return 1;
}

double SensClass::GetShellModelNME(EIsotope isotope, double ga)
{
  double ga2Fact = (ga == 1.25) ? 1 : ga*ga/(1.25*1.25);
  // Menendez J, Poves A, Caurier E and Nowacki F 2009 Nucl. Phys. A 818 139
  // Using UCOM-SRC. Xe range from Horoi M and Brown B A 2013 Novel shell-model
  // analysis of the 136xe double beta decay nuclear matrix elements (Preprint arXiv:1301.0256)
  if(isotope == kGe76)  return 2.81*ga2Fact;
  if(isotope == kTe130) return 2.65*ga2Fact;
  if(isotope == kXe136) return (1.46+2.19)/2.*ga2Fact;
  if(isotope == kNd150) {
    static bool first = true;
    if(first) {
      cout << "Sorry, no shell model NME for 150Nd" << endl;
      first = false;
    }
    return 1e-99;
  }
  cout << "Unknown isotope " << isotope << endl;
  return 1;
}

double SensClass::GetPhaseSpaceFactor(EIsotope isotope)
{
  // Kotila J and Iachello F 2012 Phys. Rev. C 85 034316
  // Multiplied by g_A^4 (with g_A = 1.25).
  if(isotope == kGe76)  return 5.77e-15/year/electron_mass_c2/electron_mass_c2;
  if(isotope == kTe130) return 3.47e-14/year/electron_mass_c2/electron_mass_c2;
  if(isotope == kXe136) return 3.56e-14/year/electron_mass_c2/electron_mass_c2;
  if(isotope == kNd150) return 1.54e-13/year/electron_mass_c2/electron_mass_c2;
  cout << "Unknown isotope " << isotope << endl;
  return 1;
}

double SensClass::GetNME(EIsotope isotope, ENMEMethod nmeMethod, double ga)
{
  if(nmeMethod == kQRPA) return GetQRPANME(isotope, ga);
  if(nmeMethod == kShellModel) return GetShellModelNME(isotope, ga);
  cout << "Unknown nmeMethod " << nmeMethod << endl;
  return 1.e-99;
}

double SensClass::GetMbb(double tHalf, EIsotope isotope, ENMEMethod nmeMethod, double ga)
{
  return 1./GetNME(isotope, nmeMethod, ga)/sqrt(GetPhaseSpaceFactor(isotope)*tHalf);
}

double SensClass::GetTHalf(double mbb, EIsotope isotope, ENMEMethod nmeMethod, double ga)
{
  return 1./pow(mbb*GetNME(isotope, nmeMethod, ga), 2)/GetPhaseSpaceFactor(isotope);
}

