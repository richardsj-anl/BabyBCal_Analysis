

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


double gaus(Double_t *x, Double_t *par){
    // [constant] * gaus([x], [mean], [sigma], norm=false)
    return par[0] * TMath::Gaus(x[0], par[1], par[2], false);
}


double Gaus_P1(Double_t *x, Double_t *par) {
    // Sum of signal and background functions
    return gaus(x, par) + poly1(x, &par[3]);
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
    //Histo->Draw();
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

    cFit->Print(imagesDir + Form("fit_Data_%s%s_%.1fGeV%s.pdf",
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

    bool mom4 = mom==4.0; bool mom4_2 = mom==4.2;
    bool mom6 = mom==6.0; bool mom6_2 = mom==6.2;
    bool mom8 = mom==8.0; bool mom8_2 = mom==8.2;

    bool partEle = ElePim==0; bool partPim = ElePim==1;

    if (mom4 && partPim) {
        fitMin = 0.55;  fitMax = 3.30;  fitMean = 1.5;  xUp = 4.80;
    }
    else if (mom6 && partPim) {
        fitMin = 1.00;  fitMax = 4.80;  fitMean = 2.5;  xUp = 7.20;
    }
    else if (mom8 && partPim) {
        fitMin = 1.25;  fitMax = 6.00;  fitMean = 3.0;  xUp = 9.60;
    }

    //Double_t low4 = 0.85;
    //Double_t high4 = 3.00;
    //Double_t low6 = 1.25;
    //Double_t high6 = 4.50;
    //Double_t low8 = 1.40;
    //Double_t high8 = 6.00;
    //Double_t low10 = 2.00;
    //Double_t high10 = 6.00;
    
    TString dataMCString[2] = {"Data", "MC"};
    TString ElePimString[2] = {"Ele", "Pim"};
    TString DepCalString[2] = {"Dep", "Cal"};
    TString kBString[6] = {"_kB0132", "_kB033", "_kB0462", "_kB0528", "_kB066", "_kB132"};
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

    cFit->Print(imagesDir + Form("fit_MC_%s%s_%.1fGeV%s%s%s.pdf",
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
    else if (momentum==6.0) {hMC->Rebin(10);}
    else if (momentum==8.0) {hMC->Rebin(10);}
    

    TF1 *fitMC = fitMCG_P1(momentum, imagesDir, hMC, ElePim, DepCal, kB, pList, Containment); 
    Double_t peakMC = fitMC->GetMaximum();

    hMC->Scale(peakData/peakMC);
    hMC->SetLineStyle(1);
    hMC->SetLineWidth(1);

    return hMC;
}


void beamProfileScanPlot () {

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

    TString imagesDir = "beamProfileScan/";

    // Check if imagesDir exists, create it if it doesn't
    if (gSystem->AccessPathName(imagesDir)) {
        std::cout << "Directory " << imagesDir << " does not exist. Creating it..." << std::endl;
        gSystem->mkdir(imagesDir, kTRUE); // kTRUE for recursive creation
    }

    TFile dataFile("FTBF_data/Pion_Analysis_Histograms_U.root", "READ");  // e calib "U"
    TFile mcFile("histograms/Histos.root", "READ");

    //--------------Histos from Henry--------------------------------------
    TH1F *hData_4 = (TH1F*) dataFile.Get("h4_Electron");  //Electron Calibration
    TH1F *hData_6 = (TH1F*) dataFile.Get("h6_Electron");  //Electron Calibration
    TH1F *hData_8 = (TH1F*) dataFile.Get("h8_Electron");  //Electron Calibration

    hData_4->Rebin(1);
    hData_6->Rebin(1);
    hData_8->Rebin(1);

    ////--------------Histos from MonteCarlo sims----------------------------
    TH1F *hPion_mc_kB132_FTFP_BERT_4_Low_U = (TH1F*) mcFile.Get("hPimCal_4.0GeV_APC_TC_kB132_Adjusted_FTFP_BERT_9.4_37.3_U");
    TH1F *hPion_mc_kB132_FTFP_BERT_4_Nom_U = (TH1F*) mcFile.Get("hPimCal_4.0GeV_APC_TC_kB132_Adjusted_FTFP_BERT_50.3_37.3_U");
    TH1F *hPion_mc_kB132_FTFP_BERT_4_High_U = (TH1F*) mcFile.Get("hPimCal_4.0GeV_APC_TC_kB132_Adjusted_FTFP_BERT_91.2_37.3_U");

    TH1F *hPion_mc_kB132_FTFP_BERT_6_Low_U = (TH1F*) mcFile.Get("hPimCal_6.0GeV_APC_TC_kB132_Adjusted_FTFP_BERT_22.0_14.3_U");
    TH1F *hPion_mc_kB132_FTFP_BERT_6_Nom_U = (TH1F*) mcFile.Get("hPimCal_6.0GeV_APC_TC_kB132_Adjusted_FTFP_BERT_29.4_14.3_U");
    TH1F *hPion_mc_kB132_FTFP_BERT_6_High_U = (TH1F*) mcFile.Get("hPimCal_6.0GeV_APC_TC_kB132_Adjusted_FTFP_BERT_36.8_14.3_U");

    TH1F *hPion_mc_kB132_FTFP_BERT_8_Low_U = (TH1F*) mcFile.Get("hPimCal_8.0GeV_APC_TC_kB132_Adjusted_FTFP_BERT_18.7_7.3_U");
    TH1F *hPion_mc_kB132_FTFP_BERT_8_Nom_U = (TH1F*) mcFile.Get("hPimCal_8.0GeV_APC_TC_kB132_Adjusted_FTFP_BERT_20.6_7.3_U");
    TH1F *hPion_mc_kB132_FTFP_BERT_8_High_U = (TH1F*) mcFile.Get("hPimCal_8.0GeV_APC_TC_kB132_Adjusted_FTFP_BERT_22.5_7.3_U");

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
    TH1F *scaledMC_4_kB132_FTFP_BERT_Low = scaledHistogram(4.0, 5, imagesDir, hPion_mc_kB132_FTFP_BERT_4_Low_U, peakData4, ElePim, DepCal, 1, O_U_u);
    scaledMC_4_kB132_FTFP_BERT_Low->Rebin(3);
    TH1F *scaledMC_4_kB132_FTFP_BERT_Nom = scaledHistogram(4.0, 5, imagesDir, hPion_mc_kB132_FTFP_BERT_4_Nom_U, peakData4, ElePim, DepCal, 1, O_U_u);
    scaledMC_4_kB132_FTFP_BERT_Nom->Rebin(3);
    TH1F *scaledMC_4_kB132_FTFP_BERT_High = scaledHistogram(4.0, 5, imagesDir, hPion_mc_kB132_FTFP_BERT_4_High_U, peakData4, ElePim, DepCal, 1, O_U_u);
    scaledMC_4_kB132_FTFP_BERT_High->Rebin(3);

    TH1F *scaledMC_6_kB132_FTFP_BERT_Low = scaledHistogram(6.0, 5, imagesDir, hPion_mc_kB132_FTFP_BERT_6_Low_U, peakData6, ElePim, DepCal, 1, O_U_u);
    scaledMC_6_kB132_FTFP_BERT_Low->Rebin(4);
    TH1F *scaledMC_6_kB132_FTFP_BERT_Nom = scaledHistogram(6.0, 5, imagesDir, hPion_mc_kB132_FTFP_BERT_6_Nom_U, peakData6, ElePim, DepCal, 1, O_U_u);
    scaledMC_6_kB132_FTFP_BERT_Nom->Rebin(4);
    TH1F *scaledMC_6_kB132_FTFP_BERT_High = scaledHistogram(6.0, 5, imagesDir, hPion_mc_kB132_FTFP_BERT_6_High_U, peakData6, ElePim, DepCal, 1, O_U_u);
    scaledMC_6_kB132_FTFP_BERT_High->Rebin(4);

    TH1F *scaledMC_8_kB132_FTFP_BERT_Low = scaledHistogram(8.0, 5, imagesDir, hPion_mc_kB132_FTFP_BERT_8_Low_U, peakData8, ElePim, DepCal, 1, O_U_u);
    scaledMC_8_kB132_FTFP_BERT_Low->Rebin(5);
    TH1F *scaledMC_8_kB132_FTFP_BERT_Nom = scaledHistogram(8.0, 5, imagesDir, hPion_mc_kB132_FTFP_BERT_8_Nom_U, peakData8, ElePim, DepCal, 1, O_U_u);
    scaledMC_8_kB132_FTFP_BERT_Nom->Rebin(5);
    TH1F *scaledMC_8_kB132_FTFP_BERT_High = scaledHistogram(8.0, 5, imagesDir, hPion_mc_kB132_FTFP_BERT_8_High_U, peakData8, ElePim, DepCal, 1, O_U_u);
    scaledMC_8_kB132_FTFP_BERT_High->Rebin(5);


    //---------------Plot--------------------------------------------------
    TCanvas *c4= new TCanvas("c4beamProfile","c4beamProfile", 1200, 800);

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

    scaledMC_4_kB132_FTFP_BERT_Low->SetLineColor(col4);
    scaledMC_4_kB132_FTFP_BERT_Low->Draw("histe same");

    scaledMC_4_kB132_FTFP_BERT_High->SetLineColor(darker10);
    scaledMC_4_kB132_FTFP_BERT_High->Draw("histe same");

    scaledMC_4_kB132_FTFP_BERT_Nom->SetLineColor(dark6);
    scaledMC_4_kB132_FTFP_BERT_Nom->SetLineWidth(3);
    scaledMC_4_kB132_FTFP_BERT_Nom->Draw("histe same");

    hData_4->SetMarkerColor(dark10);
    hData_4->SetLineColor(dark10);
    hData_4->SetLineWidth(2);
    hData_4->SetMarkerStyle(20);
    hData_4->SetMarkerSize(2);
    hData_4->Draw("p same");

    auto leg4= new TLegend(0.53, 0.61, 0.89, 0.89);
    leg4->SetBorderSize(0); leg4->SetFillStyle(0);leg4->SetTextSize(0.045);
    leg4->AddEntry(hData_4, "Data", "ep");
    leg4->AddEntry(scaledMC_4_kB132_FTFP_BERT_Low, "MC: #sigma = 9.4 mm", "el");
    leg4->AddEntry(scaledMC_4_kB132_FTFP_BERT_Nom, "MC: #sigma = 50.3 mm", "el");
    leg4->AddEntry(scaledMC_4_kB132_FTFP_BERT_High, "MC: #sigma = 91.2 mm", "el");
    leg4->Draw();

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

    c4->Print(imagesDir + "pion4GeVbeamProfile.pdf");


    TCanvas *c6= new TCanvas("c6beamProfile","c6beamProfile", 1200, 800);

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

    scaledMC_6_kB132_FTFP_BERT_Low->SetLineColor(col4);
    scaledMC_6_kB132_FTFP_BERT_Low->Draw("histe same");

    scaledMC_6_kB132_FTFP_BERT_High->SetLineColor(darker10);
    scaledMC_6_kB132_FTFP_BERT_High->Draw("histe same");

    scaledMC_6_kB132_FTFP_BERT_Nom->SetLineColor(dark6);
    scaledMC_6_kB132_FTFP_BERT_Nom->SetLineWidth(3);
    scaledMC_6_kB132_FTFP_BERT_Nom->Draw("histe same");

    hData_6->SetMarkerColor(dark10);
    hData_6->SetLineColor(dark10);
    hData_6->SetLineWidth(2);
    hData_6->SetMarkerStyle(20);
    hData_6->SetMarkerSize(2);
    hData_6->Draw("p same");

    auto leg6= new TLegend(0.53, 0.61, 0.89, 0.89);
    leg6->SetBorderSize(0); leg6->SetFillStyle(0); leg6->SetTextSize(0.045);
    leg6->AddEntry(hData_6, "Data", "ep");
    leg6->AddEntry(scaledMC_6_kB132_FTFP_BERT_Low, "MC: #sigma = 22.0 mm", "el");
    leg6->AddEntry(scaledMC_6_kB132_FTFP_BERT_Nom, "MC: #sigma = 29.4 mm", "el");
    leg6->AddEntry(scaledMC_6_kB132_FTFP_BERT_High, "MC: #sigma = 36.8 mm", "el");
    leg6->Draw();

    TLatex *l6 = new TLatex();
    l6->SetNDC();
    l6->SetTextSize(0.100);
    l6->DrawLatex(0.30, 0.80, "6 GeV");
    l6->Draw();

    lPrelim->DrawLatex(0.15, 0.40, "Preliminary Results, Baby BCAL");
    lPrelim->Draw();

    c6->Print(imagesDir + "pion6GeVbeamProfile.pdf");


    TCanvas *c8= new TCanvas("c8beamProfile","c8beamProfile", 1200, 800);

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

    scaledMC_8_kB132_FTFP_BERT_Low->SetLineColor(col4);
    scaledMC_8_kB132_FTFP_BERT_Low->Draw("histe same");

    scaledMC_8_kB132_FTFP_BERT_High->SetLineColor(darker10);
    scaledMC_8_kB132_FTFP_BERT_High->Draw("histe same");

    scaledMC_8_kB132_FTFP_BERT_Nom->SetLineColor(dark6);
    scaledMC_8_kB132_FTFP_BERT_Nom->SetLineWidth(3);
    scaledMC_8_kB132_FTFP_BERT_Nom->Draw("histe same");

    hData_8->SetMarkerColor(dark10);
    hData_8->SetLineColor(dark10);
    hData_8->SetLineWidth(2);
    hData_8->SetMarkerStyle(20);
    hData_8->SetMarkerSize(2);
    hData_8->Draw("p same");

    auto leg8= new TLegend(0.53, 0.61, 0.89, 0.89);
    leg8->SetBorderSize(0); leg8->SetFillStyle(0);leg8->SetTextSize(0.045);
    leg8->AddEntry(hData_8, "Data", "ep");
    leg8->AddEntry(scaledMC_8_kB132_FTFP_BERT_Low, "MC: #sigma = 18.7 mm", "el");
    leg8->AddEntry(scaledMC_8_kB132_FTFP_BERT_Nom, "MC: #sigma = 20.6 mm", "el");
    leg8->AddEntry(scaledMC_8_kB132_FTFP_BERT_High, "MC: #sigma = 22.6 mm", "el");
    leg8->Draw();

    TLatex *l8 = new TLatex();
    l8->SetNDC();
    l8->SetTextSize(0.100);
    l8->DrawLatex(0.30, 0.80, "8 GeV");
    l8->Draw();

    lPrelim->DrawLatex(0.15, 0.40, "Preliminary Results, Baby BCAL");
    lPrelim->Draw();

    c8->Print(imagesDir + "pion8GeVbeamProfile.pdf");

}

