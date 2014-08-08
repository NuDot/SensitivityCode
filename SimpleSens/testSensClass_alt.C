//
//
{

  //
  //This is an example macro for using SensClass.C - Author J. Detwiler.
  //Converted to a ROOT Macro by Lindley Winslow, July 24, 2014.
  //


  //Need to load some classes.
  gROOT->ProcessLine(".include \"$CLHEP_INCLUDE_DIR\"");
  gROOT->LoadMacro("SensClass.C+");

  //Make the plots pretty.
  gStyle->SetTitleXOffset(1.25);
  gStyle->SetTitleYOffset(1.25);
 

  //Choose an Isotope and Enrichment.
  Int_t isotope = SensClass::kXe136;
  Float_t enrichment = 0.906; //Phase 1: .906, Phase 2: .9077
  const char *plot_title = "xenon_sensitivity";

  //This is experiment dependent...look to J. Detwiler's original code to lear what he used.
  Float_t efficiency = .998; //Phase 1: .998, Phase 2: .93
  Float_t cPerROIty = 1./CLHEP::tonne/CLHEP::year; //Get Units Right.

  Int_t statModel = SensClass::kFC90CLSens;
  Int_t nmeModel= SensClass::kQRPA;
  
  //What Value for ga?
  //This is almost always 1.25 although there is some 
  Float_t ga=1.25; 
  //Some Ranges.
  Float_t expMin_ty = 1e-3;
  Float_t expMax_ty = 1000;
  Float_t mbbMin_meV = .99999;
  Float_t mbbMax_meV = 1000;
  Float_t tHalfMin_y = 0.9999e24;
  Float_t tHalfMax_y = 1.e30;

/*
  //No Background
  SensClass bgFreeSens(isotope, enrichment, efficiency, 0.*cPerROIty, statModel, nmeModel, ga);
  TF1 bgFreeFunc =  bgFreeSens.GetTHalfVsExposureFunction(kTRUE, SensClass::kYears);
  bgFreeFunc.SetRange(expMin_ty, expMax_ty);

  SensClass bgLowSens(isotope, enrichment, efficiency, 0.1*cPerROIty, statModel, nmeModel, ga);
  TF1 bgLowFunc =  bgLowSens.GetTHalfVsExposureFunction(kTRUE, SensClass::kYears);
  bgLowFunc.SetRange(expMin_ty, expMax_ty);

  SensClass bgMedSens(isotope, enrichment, efficiency, 1.0*cPerROIty, statModel, nmeModel, ga);
  TF1 bgMedFunc =  bgMedSens.GetTHalfVsExposureFunction(kTRUE, SensClass::kYears);
  bgMedFunc.SetRange(expMin_ty, expMax_ty);

  SensClass bgHighSens(isotope, enrichment, efficiency, 10.0*cPerROIty, statModel, nmeModel, ga);
  TF1 bgHighFunc =  bgHighSens.GetTHalfVsExposureFunction(kTRUE, SensClass::kYears);
  bgHighFunc.SetRange(expMin_ty, expMax_ty);
*/

  //Kamland-Zen background count
  SensClass bgKamSens(isotope, enrichment, efficiency, 57*cPerROIty, statModel, nmeModel, ga);
  TF1 bgKamFunc =  bgKamSens.GetTHalfVsExposureFunction(kTRUE, SensClass::kYears, SensClass::kTonneYears);
  bgKamFunc.SetRange(expMin_ty, expMax_ty);
  
  double kam_exposure = .0895; //tonne-years. Phase 1: .0895, Phase 2: .0276
  double kam_limit = bgKamFunc.Eval(kam_exposure);
  cout<<kam_limit<<endl;


  //Set Axes Labels.

  bgKamFunc.GetXaxis()->SetTitle("Exposure [tonne-yrs]");
  bgKamFunc.GetYaxis()->SetTitle("^{136}Xe T_{1/2} 90% C.L. [yrs]");
  //bgFreeFunc.GetXaxis()->SetTitle("Exposure [tonne-yrs]");
  //bgFreeFunc.GetYaxis()->SetTitle("^{136}Xe T_{1/2} 90% C.L. [yrs]");

  //Make Pretty Lines.
  //bgLowFunc.SetLineStyle(9);
  //bgMedFunc.SetLineStyle(7);
  //bgHighFunc.SetLineStyle(10);
  

  //Make a Can/vas and Draw!
  TCanvas c1("c1", "Test Sensitivity", 800, 600);
  c1.SetLogx();
  c1.SetLogy();

  bgKamFunc.Draw();
  //bgFreeFunc.Draw();
  //bgLowFunc.Draw("same");
  //bgMedFunc.Draw("same");
  //bgHighFunc.Draw("same");

  c1.SaveAs(Form("%s.pdf", plot_title));
  c1.Close();



}
