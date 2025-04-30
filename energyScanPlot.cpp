

Color_t darkerShade(Color_t color, float factor = 0.6f) {
    Float_t r, g, b;
    gROOT->GetColor(color)->GetRGB(r, g, b);
    // clamp each component to [0,1]
    if(r * factor < 1.0f) r *= factor; else r = 1.0f;
    if(g * factor < 1.0f) g *= factor; else g = 1.0f;
    if(b * factor < 1.0f) b *= factor; else b = 1.0f;
    return TColor::GetColor(r, g, b);
}


double poly1(Double_t *x, Double_t *par){
    return par[0] + par[1]*x[0];
}


double poly3(Double_t *x, Double_t *par){
    return par[0] + par[1]*x[0] + par[1]*x[0]*x[0] + par[1]*x[0]*x[0]*x[0];
}


double expInv(Double_t *x, Double_t *par){
    // [constant] * gaus([x], [mean], [sigma], norm=false)
    return par[0] * 1/TMath::Exp(x[0]);
}


double gaus(Double_t *x, Double_t *par){
    // [constant] * gaus([x], [mean], [sigma], norm=false)
    return par[0] * TMath::Gaus(x[0], par[1], par[2], false);
}


double crystalball_function(double *x, double *par){
    // [constant] * crystalball_function([x], [alpha], [n], [sigma], [mean])
    return par[0] * ROOT::Math::crystalball_function(x[0], par[1], par[2], par[3], par[4]);
}


double Gaus_ExpInv(Double_t *x, Double_t *par) {
    // Sum of signal and background functions
    return gaus(x, par) + expInv(x, &par[3]);
}


double Gaus_P1(Double_t *x, Double_t *par) {
    // Sum of signal and background functions
    return gaus(x, par) + poly1(x, &par[3]);
}


double Gaus_P2(Double_t *x, Double_t *par) {
    // Sum of signal and background functions
    return gaus(x, par) + poly3(x, &par[3]);
}


TF1* fitDataG_P1(Double_t momentum, TString imagesDir, TH1F* histo,
                 Int_t ElePim, Int_t DepCal, Int_t Containment) {

    Double_t mom = momentum;
    Double_t fitMin = 0.0;
    Double_t fitMax = 1.0;
    Double_t fitConstant = 0.015;  //1200.0
    Double_t fitMean = 1.0;  //1.4
    Double_t fitSigma = 0.5;  //0.5
    Double_t fitBkg0 = 0.01;  //1.0
    Double_t fitBkg1 = -0.01;  //0.0
    Double_t yUp = 0.04;
    Double_t xUp = 4.80;

    bool mom4 = mom==4.0; bool mom4_2 = mom==4.2;
    bool mom6 = mom==6.0; bool mom6_2 = mom==6.2;
    bool mom8 = mom==8.0; bool mom8_2 = mom==8.2;

    bool partEle = ElePim==0; bool partPim = ElePim==1;

    if (mom4 && partPim) {
        fitMin = 0.80;  fitMax = 2.50;  fitMean = 1.5;
        xUp = 4.80; yUp = 0.016; 
    }
    else if (mom6 && partPim) {
        fitMin = 1.00;  fitMax = 3.50;  fitMean = 2.5;
        xUp = 7.20; yUp = 0.012; 
    }
    else if (mom8 && partPim) {
        fitMin = 1.50;  fitMax = 4.50;  fitMean = 3.0;
        xUp = 9.60; yUp = 0.01; 
    }

    TString dataMCString[2] = {"Data", "MC"};
    TString ElePimString[2] = {"Ele", "Pim"};
    TString DepCalString[2] = {"Dep", "Cal"};
    TString ContainmentString[3] = {"_O", "_U", "_u"};
    
    Int_t nPar = 5;
    Double_t par[nPar];

    TH1F *Histo= (TH1F*) (histo->Clone("Histo"));

    TF1 *fitSum = new TF1("fitSum", Gaus_P1, fitMin, fitMax, nPar);
    fitSum->SetParameters(fitConstant, fitMean, fitSigma, fitBkg0, fitBkg1);
    fitSum->SetParNames("constant", "mean", "#sigma", "p0", "p1");
    fitSum->SetParLimits(0, 0.0, 1.0);
    Histo->Fit("fitSum", "R");
    fitSum->SetRange(fitMin, fitMax);
    fitSum->GetParameters(par);
    fitSum->SetParameters(par);
    Histo->Fit("fitSum", "R");
    fitSum->GetParameters(par);
    fitSum->SetLineColor(kOrange + 2);
    fitSum->SetLineWidth(3);

    TF1 *fitGaus = new TF1("fitGaus", gaus, fitMin, fitMax, 3);
    fitGaus->SetParameters(par);
    fitGaus->SetLineColor(kYellow + 2);
    fitGaus->SetLineWidth(3);

    TF1 *fitP1 = new TF1("fitP1", poly1, fitMin, fitMax, 2);
    fitP1->SetParameters(&par[3]);
    fitP1->SetLineColor(kBlack);
    fitP1->SetLineWidth(3);

    TCanvas *cFit = new TCanvas("cFit","cFit", 1200, 800);

    gPad->DrawFrame(0.0, 0.0, xUp, yUp, "; Reconstructed Energy (GeV); Normalized Counts");
    Histo->Draw("same");
    fitSum->Draw("same");
    fitGaus->Draw("same");
    fitP1->Draw("same");

    TLegend *legFit = new TLegend(0.66, 0.66, 0.89, 0.89);
    legFit->SetBorderSize(0); legFit->SetFillStyle(0); legFit->SetTextSize(0.03);
    legFit->AddEntry(Histo, Form("Data - %.1f GeV", mom), "ep");
    legFit->AddEntry(fitSum, "Gaus + poly1", "l");
    legFit->AddEntry(fitGaus, "Gaus", "l");
    legFit->AddEntry(fitP1, "poly1", "l");
    legFit->Draw();

    TLatex *l = new TLatex();
    l->SetNDC();
    l->SetTextSize(0.030);
    l->SetTextColor(kOrange+2);
    Double_t chi2NDF = (fitSum->GetChisquare()) / (fitSum->GetNDF());
    l->DrawLatex(0.66, 0.61, Form("Chi2/NDF = %f", chi2NDF));
    l->DrawLatex(0.66, 0.56, Form("MaximumX = %f", fitSum->GetMaximumX()));
    l->Draw();

    TLatex *name = new TLatex();
    name->SetNDC();
    name->SetTextSize(0.030);
    name->SetTextColor(kBlack);
    name->DrawLatex(0.21, 0.85, Form("Data_%s%s_%.1fGeV%s",
                                ElePimString[ElePim].Data(),
                                DepCalString[DepCal].Data(),
                                mom,
                                ContainmentString[Containment].Data()));
    name->Draw();

    cFit->Print(imagesDir + Form("fit_Data_GP1_%s%s_%.1fGeV%s.png",
                                ElePimString[ElePim].Data(),
                                DepCalString[DepCal].Data(),
                                mom,
                                ContainmentString[Containment].Data()));

    return fitSum;
}


TF1* fitMCG_P1(Double_t momentum, TString imagesDir, TH1F* histo,
             Int_t ElePim, Int_t DepCal, Int_t kB, Int_t pList, Int_t Containment) {
             //Double_t min, Double_t max,
             //Double_t constant, Double_t mean, Double_t sigma,
             //Double_t p0, Double_t p1) {

    Double_t mom = momentum;
    Double_t fitMin = 0.0;
    Double_t fitMax = 1.0;
    Double_t fitConstant = 1200.0;  //1200.0
    Double_t fitMean = 1.0;  //1.4
    Double_t fitSigma = 0.5;  //0.5
    Double_t fitBkg0 = 0.01;  //1.0
    Double_t fitBkg1 = -0.01;  //0.0
    Double_t xUp = 4.80;
    Double_t yUp = 600.0;

    bool mom4 = mom==4.0; bool mom4_2 = mom==4.2; bool mom4_4 = mom==4.4; bool mom4_5 = mom==4.5; bool mom4_6 = mom==4.6; bool mom4_8 = mom==4.8;
    bool mom5 = mom==5.0;
    bool mom6 = mom==6.0; bool mom6_2 = mom==6.2;
    bool mom8 = mom==8.0; bool mom8_2 = mom==8.2;

    bool partEle = ElePim==0; bool partPim = ElePim==1;

    if (mom4 && partPim) {
        fitMin = 0.55;  fitMax = 3.00;  fitMean = 1.5;  xUp = 4.80;
    }
    else if (mom4_2 && partPim) {
        fitMin = 0.55;  fitMax = 3.25;  fitMean = 1.5;  xUp = 4.80;
    }
    else if (mom4_4 && partPim) {
        fitMin = 0.55;  fitMax = 3.25;  fitMean = 1.5;  xUp = 4.80;
    }
    else if (mom4_5 && partPim) {
        fitMin = 0.65;  fitMax = 3.50;  fitMean = 1.5;  xUp = 4.80;
    }
    else if (mom4_6 && partPim) {
        fitMin = 0.65;  fitMax = 3.00;  fitMean = 1.5;  xUp = 4.80;
    }
    else if (mom4_8 && partPim) {
        fitMin = 0.65;  fitMax = 3.00;  fitMean = 1.5;  xUp = 4.80;
    }
    else if (mom5 && partPim) {
        fitMin = 0.65;  fitMax = 3.00;  fitMean = 1.5;  xUp = 4.80;
    }
    else if (mom6 && partPim) {
        fitMin = 0.75;  fitMax = 4.80;  fitMean = 2.5;  xUp = 7.20;
    }
    else if (mom6_2 && partPim) {
        fitMin = 1.00;  fitMax = 5.00;  fitMean = 2.5;  xUp = 7.20;
    }
    else if (mom8 && partPim) {
        fitMin = 1.25;  fitMax = 6.00;  fitMean = 3.0;  xUp = 9.60;
    }
    else if (mom8_2 && partPim) {
        fitMin = 1.25;  fitMax = 6.00;  fitMean = 3.0;  xUp = 9.60;
    }

    TString dataMCString[2] = {"Data", "MC"};
    TString ElePimString[2] = {"Ele", "Pim"};
    TString DepCalString[2] = {"Dep", "Cal"};
    TString kBString[4] = {"_kB0132", "_kB033", "_kB066", "_kB132"};
    TString pListString[5] = {"_QGSP_BERT", "_FTFP_BERT", "_QGSP_BIC", "_QGSP_BERT_HP", "_FTFP_BERT_HP"};
    TString ContainmentString[3] = {"_O", "_U", "_u"};
    
    Int_t nPar = 5;
    Double_t par[nPar];

    TH1F *Histo= (TH1F*) (histo->Clone("Histo"));

    TF1 *fitSum = new TF1("fitSum", Gaus_P1, fitMin, fitMax, nPar);
    fitSum->SetParameters(fitConstant, fitMean, fitSigma, fitBkg0, fitBkg1);
    fitSum->SetParNames("constant", "mean", "#sigma", "p0", "p1");
    fitSum->SetParLimits(0, 0.0, 3000.0);
    Histo->Fit("fitSum", "R");
    fitSum->SetRange(fitMin, fitMax);
    fitSum->GetParameters(par);
    fitSum->SetParameters(par);
    Histo->Fit("fitSum", "R");
    fitSum->GetParameters(par);
    fitSum->SetLineColor(kOrange + 2);
    fitSum->SetLineWidth(3);

    TF1 *fitGaus = new TF1("fitGaus", gaus, fitMin, fitMax, 3);
    fitGaus->SetParameters(par);
    fitGaus->SetLineColor(kYellow + 2);
    fitGaus->SetLineWidth(3);

    TF1 *fitP1 = new TF1("fitP1", poly1, fitMin, fitMax, 2);
    fitP1->SetParameters(&par[3]);
    fitP1->SetLineColor(kBlack);
    fitP1->SetLineWidth(3);

    TCanvas *cFit = new TCanvas("cFit","cFit", 1200, 800);

    gPad->DrawFrame(0.0, 0.0, xUp, yUp, "; Reconstructed Energy (GeV); Normalized Counts");
    //Histo->Draw();
    Histo->Draw("same");
    fitSum->Draw("same");
    fitGaus->Draw("same");
    fitP1->Draw("same");

    TLegend *legFit = new TLegend(0.66, 0.66, 0.89, 0.89);
    legFit->SetBorderSize(0); legFit->SetFillStyle(0); legFit->SetTextSize(0.03);
    legFit->AddEntry(Histo, Form("MC - %.1f GeV", mom), "ep");
    legFit->AddEntry(fitSum, "Gaus + poly1", "l");
    legFit->AddEntry(fitGaus, "Gaus", "l");
    legFit->AddEntry(fitP1, "poly1", "l");
    legFit->Draw();

    TLatex *l = new TLatex();
    l->SetNDC();
    l->SetTextSize(0.030);
    l->SetTextColor(kOrange+2);
    Double_t chi2NDF = (fitSum->GetChisquare()) / (fitSum->GetNDF());
    l->DrawLatex(0.66, 0.61, Form("Chi2/NDF = %f", chi2NDF));
    l->DrawLatex(0.66, 0.56, Form("MaximumX = %f", fitSum->GetMaximumX()));
    l->Draw();

    TLatex *name = new TLatex();
    name->SetNDC();
    name->SetTextSize(0.030);
    name->SetTextColor(kBlack);
    name->DrawLatex(0.21, 0.85, Form("MC_%s%s_%.1fGeV%s%s%s",
                                ElePimString[ElePim].Data(),
                                DepCalString[DepCal].Data(),
                                mom,
                                kBString[kB].Data(),
                                pListString[pList].Data(),
                                ContainmentString[Containment].Data()));
    name->Draw();

    cFit->Print(imagesDir + Form("fit_MC_GP1_%s%s_%.1fGeV%s%s%s.png",
                                ElePimString[ElePim].Data(),
                                DepCalString[DepCal].Data(),
                                mom,
                                kBString[kB].Data(),
                                pListString[pList].Data(),
                                ContainmentString[Containment].Data()));

    return fitSum;
}


TH1F* scaledHistogram(Double_t momentum, Int_t kB, TString imagesDir, TH1F* MC, Double_t peakData,
                      Int_t ElePim, Int_t DepCal, Int_t pList, Int_t Containment) {

    TH1F *hMC= (TH1F*) (MC->Clone("hMC"));

    if (momentum==4.0) {hMC->Rebin(10);}
    else if (momentum==4.2) {hMC->Rebin(10);}
    else if (momentum==4.4) {hMC->Rebin(10);}
    else if (momentum==4.5) {hMC->Rebin(10);}
    else if (momentum==4.6) {hMC->Rebin(10);}
    else if (momentum==4.8) {hMC->Rebin(10);}
    else if (momentum==5.0) {hMC->Rebin(10);}
    else if (momentum==6.0) {hMC->Rebin(10);}
    else if (momentum==6.2) {hMC->Rebin(10);}
    else if (momentum==8.0) {hMC->Rebin(10);}
    else if (momentum==8.2) {hMC->Rebin(10);}
    

    TF1 *fitMC = fitMCG_P1(momentum, imagesDir, hMC, ElePim, DepCal, kB, pList, Containment); 
    Double_t peakMC = fitMC->GetMaximum();

    hMC->Scale(peakData/peakMC);
    hMC->SetLineStyle(1);
    hMC->SetLineWidth(2);

    return hMC;
}


void energyScanPlot () {

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(111);

    int col4 = TColor::GetColor("#E69F00");
    int dark4 = darkerShade(col4);
    int darker4 = darkerShade(dark4);
    int col6 = TColor::GetColor("#56B4E9");
    int dark6 = darkerShade(col6);
    int darker6 = darkerShade(dark6);
    int col8 = TColor::GetColor("#009E73");
    int dark8 = darkerShade(col8);
    int darker8 = darkerShade(dark8);
    int col10 = TColor::GetColor("#D55E00");
    int dark10 = darkerShade(col10);
    int darker10 = darkerShade(dark10);

    Double_t low4 = 1.00;
    Double_t high4 = 3.00;
    Double_t low6 = 1.25;
    Double_t high6 = 4.50;
    Double_t low8 = 1.40;
    Double_t high8 = 6.00;
    Double_t low10 = 2.00;
    Double_t high10 = 6.00;

    //TString imagesDir = "~/Projects/ANL/bbcal/images_newGeo/pim/v2.1/paper/";
    TString imagesDir = "~/Projects/ANL/bbcal/images_newGeo/pim/v2.1/paper/energyScan/updated/";

    TFile dataFile("~/Projects/ANL/bbcal/FTBF_data/Pion_Analysis_Histograms_U.root", "READ");  // e calib "U"
    //TFile dataFile("~/Projects/ANL/bbcal/FTBF_data/Pion_Analysis_Histograms.root", "READ");  // e calib "O"
    //TFile dataFile("~/Projects/ANL/bbcal/FTBF_data/Pion_Analysis_Histograms_NN.root", "READ");  // e calib "O" and nearest neighbor
    TFile mcFile("~/Projects/ANL/bbcal/histos_newGeo/v2.1/paper/paperHistosNoGeoCut.root", "READ");
    //TFile mcFile("~/Projects/ANL/bbcal/histos_newGeo/v2.1/paper/paperHistosAll4.root", "READ");
    //TFile mcFile("~/Projects/ANL/bbcal/histos_newGeo/v2.1/paper/paperHistos2.root", "READ");

    //--------------Histos from Henry--------------------------------------
    TH1F *hData_4 = (TH1F*) dataFile.Get("h4_Electron");
    TH1F *hData_6 = (TH1F*) dataFile.Get("h6_Electron");
    TH1F *hData_8 = (TH1F*) dataFile.Get("h8_Electron");

    hData_4->Rebin(1);
    hData_6->Rebin(1);
    hData_8->Rebin(1);

    ////--------------Histos from MonteCarlo sims----------------------------
    TH1F *hPion_mc_kB132_FTFP_BERT_4_0_U = (TH1F*) mcFile.Get("hPimCal_4.0GeV_APC_TC_kB132_Adjusted_FTFP_BERT_50.3_37.3_U");
    TH1F *hPion_mc_kB132_FTFP_BERT_4_2_U = (TH1F*) mcFile.Get("hPimCal_4.2GeV_APC_TC_kB132_Adjusted_FTFP_BERT_50.3_37.3_U");
    TH1F *hPion_mc_kB132_FTFP_BERT_4_4_U = (TH1F*) mcFile.Get("hPimCal_4.4GeV_APC_TC_kB132_Adjusted_FTFP_BERT_50.3_37.3_U");
    TH1F *hPion_mc_kB132_FTFP_BERT_4_5_U = (TH1F*) mcFile.Get("hPimCal_4.5GeV_APC_TC_kB132_Adjusted_FTFP_BERT_50.3_37.3_U");
    //TH1F *hPion_mc_kB132_FTFP_BERT_4_6_U = (TH1F*) mcFile.Get("hPimCal_4.6GeV_APC_TC_kB132_Adjusted_FTFP_BERT_U");
    //TH1F *hPion_mc_kB132_FTFP_BERT_4_8_U = (TH1F*) mcFile.Get("hPimCal_4.8GeV_APC_TC_kB132_Adjusted_FTFP_BERT_U");
    //TH1F *hPion_mc_kB132_FTFP_BERT_5_0_U = (TH1F*) mcFile.Get("hPimCal_5.0GeV_APC_TC_kB132_Adjusted_FTFP_BERT_U");

    TH1F *hPion_mc_kB132_FTFP_BERT_6_0_U = (TH1F*) mcFile.Get("hPimCal_6.0GeV_APC_TC_kB132_Adjusted_FTFP_BERT_29.4_14.3_U");
    TH1F *hPion_mc_kB132_FTFP_BERT_6_2_U = (TH1F*) mcFile.Get("hPimCal_6.2GeV_APC_TC_kB132_Adjusted_FTFP_BERT_29.4_14.3_U");

    TH1F *hPion_mc_kB132_FTFP_BERT_8_0_U = (TH1F*) mcFile.Get("hPimCal_8.0GeV_APC_TC_kB132_Adjusted_FTFP_BERT_20.6_7.3_U");
    TH1F *hPion_mc_kB132_FTFP_BERT_8_2_U = (TH1F*) mcFile.Get("hPimCal_8.2GeV_APC_TC_kB132_Adjusted_FTFP_BERT_20.6_7.3_U");

    //----------------Fit Data Histos-------------------------------------------
    Int_t ElePim = 1; Int_t DepCal = 1; Int_t O_U_u = 1;

    TF1* fitData4 = fitDataG_P1(4.0, imagesDir, hData_4, ElePim, DepCal, O_U_u);
    Double_t peakData4 = fitData4->GetMaximum();
    hData_4->Rebin(3);
    std::cout << Form("4GeV, TF1::Maximum() =  %f", peakData4) << std::endl;

    TF1* fitData6 = fitDataG_P1(6.0, imagesDir, hData_6, ElePim, DepCal, O_U_u);
    Double_t peakData6 = fitData6->GetMaximum();
    hData_6->Rebin(4);
    std::cout << Form("6GeV, TF1::Maximum() =  %f", peakData6) << std::endl;

    TF1* fitData8 = fitDataG_P1(8.0, imagesDir, hData_8, ElePim, DepCal, O_U_u);
    Double_t peakData8 = fitData8->GetMaximum();
    hData_8->Rebin(5);
    std::cout << Form("8GeV, TF1::Maximum() =  %f", peakData8) << std::endl;

    //----------------Scale MC Histos-------------------------------------------
    Int_t kB = 3; Int_t pList = 1;
    TH1F *scaledMC_4_0 = scaledHistogram(4.0, kB, imagesDir, hPion_mc_kB132_FTFP_BERT_4_0_U, peakData4, ElePim, DepCal, pList, O_U_u);
    scaledMC_4_0->Rebin(3);
    TH1F *scaledMC_4_2 = scaledHistogram(4.2, kB, imagesDir, hPion_mc_kB132_FTFP_BERT_4_2_U, peakData4, ElePim, DepCal, pList, O_U_u);
    scaledMC_4_2->Rebin(3);
    TH1F *scaledMC_4_4 = scaledHistogram(4.4, kB, imagesDir, hPion_mc_kB132_FTFP_BERT_4_4_U, peakData4, ElePim, DepCal, pList, O_U_u);
    scaledMC_4_4->Rebin(3);
    TH1F *scaledMC_4_5 = scaledHistogram(4.5, kB, imagesDir, hPion_mc_kB132_FTFP_BERT_4_5_U, peakData4, ElePim, DepCal, pList, O_U_u);
    scaledMC_4_5->Rebin(3);
    //TH1F *scaledMC_4_6 = scaledHistogram(4.6, kB, imagesDir, hPion_mc_kB132_FTFP_BERT_4_6_U, peakData4, ElePim, DepCal, pList, O_U_u);
    //TH1F *scaledMC_4_8 = scaledHistogram(4.8, kB, imagesDir, hPion_mc_kB132_FTFP_BERT_4_8_U, peakData4, ElePim, DepCal, pList, O_U_u);
    //TH1F *scaledMC_5_0 = scaledHistogram(5.0, kB, imagesDir, hPion_mc_kB132_FTFP_BERT_5_0_U, peakData4, ElePim, DepCal, pList, O_U_u);

    TH1F *scaledMC_6_0 = scaledHistogram(6.0, kB, imagesDir, hPion_mc_kB132_FTFP_BERT_6_0_U, peakData6, ElePim, DepCal, pList, O_U_u);
    scaledMC_6_0->Rebin(4);
    TH1F *scaledMC_6_2 = scaledHistogram(6.2, kB, imagesDir, hPion_mc_kB132_FTFP_BERT_6_2_U, peakData6, ElePim, DepCal, pList, O_U_u);
    scaledMC_6_2->Rebin(4);

    TH1F *scaledMC_8_0 = scaledHistogram(8.0, kB, imagesDir, hPion_mc_kB132_FTFP_BERT_8_0_U, peakData8, ElePim, DepCal, pList, O_U_u);
    scaledMC_8_0->Rebin(5);
    TH1F *scaledMC_8_2 = scaledHistogram(8.2, kB, imagesDir, hPion_mc_kB132_FTFP_BERT_8_2_U, peakData8, ElePim, DepCal, pList, O_U_u);
    scaledMC_8_2->Rebin(5);

    //---------------Plot--------------------------------------------------
    TCanvas *c4energy = new TCanvas("c4energy","c4energy", 1200, 800);

    // Set pad margins before drawing
    gPad->SetLeftMargin(0.13);   // Increase left margin for Y-axis title
    gPad->SetBottomMargin(0.12); // Increase bottom margin for X-axis title

    Double_t xUp = 4.80;
    Double_t yUp = 0.045;
    // Create a dummy histogram to define the axes
    TH1F *hDummy4 = new TH1F("hDummy4", ";Reconstructed Energy (GeV);Normalized Counts", 100, 0.0, xUp);
    hDummy4->SetMaximum(yUp);
    hDummy4->SetMinimum(0.0);
    hDummy4->Draw("AXIS"); // Draw only the axes

    // Adjust axis label sizes
    hDummy4->GetXaxis()->SetTitleSize(0.05); // X-axis label size
    hDummy4->GetYaxis()->SetTitleSize(0.05); // Y-axis label size
    //gPad->DrawFrame(0.0, 0.0, 4.80, yUp, "; Reconstructed Energy (GeV); Normalized Counts");

    scaledMC_4_2->SetLineColor(col4);
    scaledMC_4_2->SetLineWidth(1);
    scaledMC_4_2->Draw("histe same");

    scaledMC_4_4->SetLineColor(col10);
    scaledMC_4_4->SetLineWidth(1);
    scaledMC_4_4->Draw("histe same");

    scaledMC_4_5->SetLineColor(darker10);
    scaledMC_4_5->SetLineWidth(1);
    scaledMC_4_5->Draw("histe same");

    scaledMC_4_0->SetLineColor(dark6);
    scaledMC_4_0->SetLineWidth(3);
    scaledMC_4_0->Draw("histe same");

    hData_4->SetMarkerColor(dark10);
    hData_4->SetLineColor(dark10);
    hData_4->SetLineWidth(2);
    hData_4->SetMarkerStyle(20);
    hData_4->SetMarkerSize(2);
    hData_4->Draw("p same");

    auto leg4energy = new TLegend(0.58, 0.55, 0.89, 0.89);
    leg4energy->SetBorderSize(0); leg4energy->SetFillStyle(0);leg4energy->SetTextSize(0.06);
    leg4energy->AddEntry(hData_4, "Data", "ep");
    leg4energy->AddEntry(scaledMC_4_0, "MC", "el");
    leg4energy->AddEntry(scaledMC_4_2, "MC +200 MeV", "el");
    leg4energy->AddEntry(scaledMC_4_4, "MC +400 MeV", "el");
    leg4energy->AddEntry(scaledMC_4_5, "MC +500 MeV", "el");
    leg4energy->Draw();

    TLatex *l4 = new TLatex();
    l4->SetNDC();
    l4->SetTextSize(0.100);
    l4->DrawLatex(0.30, 0.80, "4 GeV");
    l4->Draw();

    TLatex *lPrelim = new TLatex();
    lPrelim->SetNDC();
    lPrelim->SetTextSize(0.080);
    lPrelim->SetTextColorAlpha(kBlack, 0.35);
    lPrelim->DrawLatex(0.15, 0.40, "Preliminary Results, Baby BCAL");
    lPrelim->Draw();

    c4energy->Print(imagesDir + "pion4GeVenergy.pdf");


    TCanvas *c6energy = new TCanvas("c6energy","c6energy", 1200, 800);

    // Set pad margins before drawing
    gPad->SetLeftMargin(0.13);   // Increase left margin for Y-axis title
    gPad->SetBottomMargin(0.12); // Increase bottom margin for X-axis title

    xUp = 7.20;
    yUp = 0.040;
    // Create a dummy histogram to define the axes
    TH1F *hDummy6 = new TH1F("hDummy6", ";Reconstructed Energy (GeV);Normalized Counts", 100, 0.0, xUp);
    hDummy6->SetMaximum(yUp);
    hDummy6->SetMinimum(0.0);
    hDummy6->Draw("AXIS"); // Draw only the axes

    // Adjust axis label sizes
    hDummy6->GetXaxis()->SetTitleSize(0.05); // X-axis label size
    hDummy6->GetYaxis()->SetTitleSize(0.05); // Y-axis label size

    scaledMC_6_2->SetLineColor(col4);
    scaledMC_6_2->SetLineWidth(1);
    scaledMC_6_2->Draw("histe same");

    scaledMC_6_0->SetLineColor(dark6);
    scaledMC_6_0->SetLineWidth(3);
    scaledMC_6_0->Draw("histe same");

    hData_6->SetMarkerColor(dark10);
    hData_6->SetLineColor(dark10);
    hData_6->SetLineWidth(2);
    hData_6->SetMarkerStyle(20);
    hData_6->SetMarkerSize(2);
    hData_6->Draw("p same");

    auto leg6energy = new TLegend(0.58, 0.66, 0.89, 0.89);
    leg6energy->SetBorderSize(0); leg6energy->SetFillStyle(0);leg6energy->SetTextSize(0.06);
    leg6energy->AddEntry(hData_6, "Data", "ep");
    leg6energy->AddEntry(scaledMC_6_0, "MC", "el");
    leg6energy->AddEntry(scaledMC_6_2, "MC +200 MeV", "el");
    leg6energy->Draw();

    TLatex *l6 = new TLatex();
    l6->SetNDC();
    l6->SetTextSize(0.100);
    l6->DrawLatex(0.30, 0.80, "6 GeV");
    l6->Draw();

    lPrelim->DrawLatex(0.15, 0.40, "Preliminary Results, Baby BCAL");
    lPrelim->Draw();

    c6energy->Print(imagesDir + "pion6GeVenergy.pdf");

    TCanvas *c8energy = new TCanvas("c8energy","c8energy", 1200, 800);

    // Set pad margins before drawing
    gPad->SetLeftMargin(0.13);   // Increase left margin for Y-axis title
    gPad->SetBottomMargin(0.12); // Increase bottom margin for X-axis title

    xUp = 9.60;
    yUp = 0.035;
    // Create a dummy histogram to define the axes
    TH1F *hDummy8 = new TH1F("hDummy8", ";Reconstructed Energy (GeV);Normalized Counts", 100, 0.0, xUp);
    hDummy8->SetMaximum(yUp);
    hDummy8->SetMinimum(0.0);
    hDummy8->Draw("AXIS"); // Draw only the axes

    // Adjust axis label sizes
    hDummy8->GetXaxis()->SetTitleSize(0.05); // X-axis label size
    hDummy8->GetYaxis()->SetTitleSize(0.05); // Y-axis label size

    scaledMC_8_2->SetLineColor(col4);
    scaledMC_8_2->SetLineWidth(1);
    scaledMC_8_2->Draw("histe same");

    scaledMC_8_0->SetLineColor(dark6);
    scaledMC_8_0->SetLineWidth(3);
    scaledMC_8_0->Draw("histe same");

    hData_8->SetMarkerColor(dark10);
    hData_8->SetLineColor(dark10);
    hData_8->SetLineWidth(2);
    hData_8->SetMarkerStyle(20);
    hData_8->SetMarkerSize(2);
    hData_8->Draw("p same");

    auto leg8energy = new TLegend(0.58, 0.66, 0.89, 0.89);
    leg8energy->SetBorderSize(0); leg8energy->SetFillStyle(0);leg8energy->SetTextSize(0.06);
    leg8energy->AddEntry(hData_8, "Data", "ep");
    leg8energy->AddEntry(scaledMC_8_0, "MC", "el");
    leg8energy->AddEntry(scaledMC_8_2, "MC +200 MeV", "el");
    leg8energy->Draw();

    TLatex *l8 = new TLatex();
    l8->SetNDC();
    l8->SetTextSize(0.100);
    l8->DrawLatex(0.30, 0.80, "8 GeV");
    l8->Draw();

    lPrelim->DrawLatex(0.15, 0.40, "Preliminary Results, Baby BCAL");
    lPrelim->Draw();

    c8energy->Print(imagesDir + "pion8GeVenergy.pdf");

}

