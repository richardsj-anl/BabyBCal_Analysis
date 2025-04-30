#include <TF1.h>
#include <TH1.h>
#include <TMath.h>
#include <TColor.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>

Color_t darkerShade(Color_t color, float factor = 0.6f) {
    Float_t r, g, b;
    TColor* col = gROOT->GetColor(color);
    if (!col) {
        std::cerr << "Invalid color index: " << color << std::endl;
        return color;
    }
    col->GetRGB(r, g, b);
    if (r * factor < 1.0f) r *= factor; else r = 1.0f;
    if (g * factor < 1.0f) g *= factor; else g = 1.0f;
    if (b * factor < 1.0f) b *= factor; else b = 1.0f;
    return TColor::GetColor(r, g, b);
}

Double_t crystalball_function(Double_t *x, Double_t *par) {
    return par[0] * ROOT::Math::crystalball_function(x[0], par[1], par[2], par[3], par[4]);
}

Double_t res(Double_t *x, Double_t *par) {
    auto rm = TMath::Sqrt(x[0]);
    return TMath::Sqrt((par[0]/rm * par[0]/rm) + (par[1] * par[1]));
}

void eResCB(TF1& fit, Double_t& ERes, Double_t& eERes) {
    Double_t Sig = fit.GetParameter(3);
    Double_t eSig = fit.GetParError(3);
    Double_t E = fit.GetParameter(4);
    Double_t eE = fit.GetParError(4);
    if (E == 0 || Sig == 0) {
        ERes = 0;
        eERes = 0;
        std::cerr << "Warning: Division by zero avoided in eResCB" << std::endl;
        return;
    }
    ERes = Sig / E;
    eERes = ERes * TMath::Sqrt((eSig/Sig)*(eSig/Sig) + (eE/E)*(eE/E));
}

TF1* fitCB(Int_t mom, TH1F* hDep, Double_t* par, Double_t& ERes, Double_t& eERes) {
    TString fitName = Form("fit_%dGeV", mom);
    TF1* DepFitCB = new TF1(fitName, crystalball_function, 0.060 * mom, 0.12 * mom, 5);
    if (mom == 4) {
        DepFitCB->SetParameters(4000.0, 1.00, 1.0, 0.01, 0.4);
    } else if (mom == 6) {
        DepFitCB->SetParameters(4000.0, 1.00, 1.0, 0.01, 0.6);
    } else if (mom == 8) {
        DepFitCB->SetParameters(4000.0, 1.00, 1.0, 0.02, 0.8);
    } else if (mom == 10) {
        DepFitCB->SetParameters(4000.0, 1.00, 1.0, 0.02, 0.9);
    }
    DepFitCB->SetParNames("const", "#alpha", "n", "#sigma", "mean");
    //DepFitCB->SetParLimits(0, 3000, 6000);
    DepFitCB->SetParLimits(4, 0.01 * mom, 0.12 * mom);
    hDep->Fit(DepFitCB, "R");
    DepFitCB->GetParameters(par);
    DepFitCB->SetParLimits(4, 0.01 * mom, 0.12 * mom);
    DepFitCB->SetParameters(par);
    hDep->Fit(DepFitCB, "R");
    eResCB(*DepFitCB, ERes, eERes);
    std::cout << "Momentum " << mom << " GeV - ERes: " << ERes << ", eERes: " << eERes << std::endl;
    return DepFitCB;
}

void eResPlot() {
    gStyle->SetOptFit(111);
    // Define colors
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

    // File and directory
    TString imagesDir = "ERes/";

    // Check if imagesDir exists, create it if it doesn't
    if (gSystem->AccessPathName(imagesDir)) {
        std::cout << "Directory " << imagesDir << " does not exist. Creating it..." << std::endl;
        gSystem->mkdir(imagesDir, kTRUE); // kTRUE for recursive creation
    }

    TFile mcFile("histograms/Histos.root", "READ");
    if (!mcFile.IsOpen()) {
        std::cerr << "Error: Could not open Histo file!" << std::endl;
        return;
    }

    // Load histograms
    TH1F *hEleDep_4_sigE0 = (TH1F*)mcFile.Get("hEleDep_4.0GeV_APC_TC_kB132_sigE0_FTFP_BERT_50.3_U");
    TH1F *hEleDep_6_sigE0 = (TH1F*)mcFile.Get("hEleDep_6.0GeV_APC_TC_kB132_sigE0_FTFP_BERT_29.4_U");
    TH1F *hEleDep_8_sigE0 = (TH1F*)mcFile.Get("hEleDep_8.0GeV_APC_TC_kB132_sigE0_FTFP_BERT_20.6_U");
    TH1F *hEleDep_10_sigE0 = (TH1F*)mcFile.Get("hEleDep_10.0GeV_APC_TC_kB132_sigE0_FTFP_BERT_17.9_U");

    TH1F *hEleDep_4_Cited = (TH1F*)mcFile.Get("hEleDep_4.0GeV_APC_TC_kB132_Cited_FTFP_BERT_50.3_U");
    TH1F *hEleDep_6_Cited = (TH1F*)mcFile.Get("hEleDep_6.0GeV_APC_TC_kB132_Cited_FTFP_BERT_29.4_U");
    TH1F *hEleDep_8_Cited = (TH1F*)mcFile.Get("hEleDep_8.0GeV_APC_TC_kB132_Cited_FTFP_BERT_20.6_U");
    TH1F *hEleDep_10_Cited = (TH1F*)mcFile.Get("hEleDep_10.0GeV_APC_TC_kB132_Cited_FTFP_BERT_17.9_U");

    TH1F *hEleDep_4_Adjusted = (TH1F*)mcFile.Get("hEleDep_4.0GeV_APC_TC_kB132_Adjusted_FTFP_BERT_50.3_U");
    TH1F *hEleDep_6_Adjusted = (TH1F*)mcFile.Get("hEleDep_6.0GeV_APC_TC_kB132_Adjusted_FTFP_BERT_29.4_U");
    TH1F *hEleDep_8_Adjusted = (TH1F*)mcFile.Get("hEleDep_8.0GeV_APC_TC_kB132_Adjusted_FTFP_BERT_20.6_U");  // Fix me when sims are complete
    TH1F *hEleDep_10_Adjusted = (TH1F*)mcFile.Get("hEleDep_10.0GeV_APC_TC_kB132_Adjusted_FTFP_BERT_17.9_U");  // Fix me when sims are complete

    // Check histograms
    if (!hEleDep_4_sigE0 || !hEleDep_6_sigE0 || !hEleDep_8_sigE0 || !hEleDep_10_sigE0 ||
        !hEleDep_4_Cited || !hEleDep_6_Cited || !hEleDep_8_Cited || !hEleDep_10_Cited ||
        !hEleDep_4_Adjusted || !hEleDep_6_Adjusted || !hEleDep_8_Adjusted || !hEleDep_10_Adjusted) {
        std::cerr << "Error: One or more histograms not found!" << std::endl;
        mcFile.Close();
        return;
    }

    // Arrays for plotting
    const Int_t nMom = 4;
    Double_t momValues[nMom] = {4, 6, 8, 10};
    Double_t eMomValues[nMom] = {0, 0, 0, 0};
    Double_t EResSigE0[nMom], eEResSigE0[nMom];
    Double_t EResCited[nMom], eEResCited[nMom];
    Double_t EResAdjusted[nMom], eEResAdjusted[nMom];
    TH1F* hDepsSigE0[nMom] = {hEleDep_4_sigE0, hEleDep_6_sigE0, hEleDep_8_sigE0, hEleDep_10_sigE0};
    TH1F* hDepsCited[nMom] = {hEleDep_4_Cited, hEleDep_6_Cited, hEleDep_8_Cited, hEleDep_10_Cited};
    TH1F* hDepsAdjusted[nMom] = {hEleDep_4_Adjusted, hEleDep_6_Adjusted, hEleDep_8_Adjusted, hEleDep_10_Adjusted};
    TF1* fitsSigE0[nMom];
    TF1* fitsCited[nMom];
    TF1* fitsAdjusted[nMom];

    // Rebin histograms (rebinFactor = 1 means no rebinning)
    const Int_t rebinFactor = 1;
    for (Int_t i = 0; i < nMom; i++) {
        hDepsSigE0[i]->Rebin(rebinFactor);
        hDepsCited[i]->Rebin(rebinFactor);
        hDepsAdjusted[i]->Rebin(rebinFactor);
    }

    // Fit and calculate resolutions, store fits
    Double_t par[5];
    for (Int_t i = 0; i < nMom; i++) {
        fitsSigE0[i] = fitCB(momValues[i], hDepsSigE0[i], par, EResSigE0[i], eEResSigE0[i]);
        fitsCited[i] = fitCB(momValues[i], hDepsCited[i], par, EResCited[i], eEResCited[i]);
        fitsAdjusted[i] = fitCB(momValues[i], hDepsAdjusted[i], par, EResAdjusted[i], eEResAdjusted[i]);
    }

    // Plot individual fits with legends
    for (Int_t i = 0; i < nMom; i++) {
        // sigE0 plot
        TCanvas* cSigE0 = new TCanvas(Form("cSigE0_%.0f", momValues[i]), Form("Fit sigE0 %.0f GeV", momValues[i]), 800, 600);
        hDepsSigE0[i]->SetTitle(Form("Energy Deposition at %.0f GeV (sigE0);Energy (GeV);Counts", momValues[i]));
        hDepsSigE0[i]->GetXaxis()->SetRangeUser(0, 0.2 * momValues[i]);
        hDepsSigE0[i]->SetLineColor(kBlack);
        hDepsSigE0[i]->Draw("HIST");
        fitsSigE0[i]->SetLineColor(col6);
        fitsSigE0[i]->SetLineWidth(2);
        fitsSigE0[i]->Draw("SAME");

        TLegend* legSigE0 = new TLegend(0.1, 0.7, 0.3, 0.9);
        legSigE0->SetBorderSize(0);
        legSigE0->SetFillStyle(0);
        legSigE0->AddEntry(hDepsSigE0[i], "Data", "l");
        legSigE0->AddEntry(fitsSigE0[i], "Crystal Ball Fit", "l");
        legSigE0->Draw();

        cSigE0->SaveAs(imagesDir + Form("fit_sigE0_%.0f.png", momValues[i]));
        delete legSigE0;
        delete cSigE0;

        // Cited plot
        TCanvas* cCited = new TCanvas(Form("cCited_%.0f", momValues[i]), Form("Fit Cited %.0f GeV", momValues[i]), 800, 600);
        gPad->SetLogy();
        hDepsCited[i]->SetTitle(Form("Energy Deposition at %.0f GeV (Cited);Energy (GeV);Counts", momValues[i]));
        hDepsCited[i]->GetXaxis()->SetRangeUser(0, 0.2 * momValues[i]);
        hDepsCited[i]->SetLineColor(kBlack);
        hDepsCited[i]->Draw("HIST");
        fitsCited[i]->SetLineColor(col10);
        fitsCited[i]->SetLineWidth(2);
        fitsCited[i]->Draw("SAME");

        TLegend* legCited = new TLegend(0.1, 0.7, 0.3, 0.9);
        legCited->SetBorderSize(0);
        legCited->SetFillStyle(0);
        legCited->AddEntry(hDepsCited[i], "Data", "l");
        legCited->AddEntry(fitsCited[i], "Crystal Ball Fit", "l");
        legCited->Draw();

        cCited->SaveAs(imagesDir + Form("fit_Cited_%.0f.png", momValues[i]));
        delete legCited;
        delete cCited;

        // Adjusted plot
        TCanvas* cAdjusted = new TCanvas(Form("cAdjusted_%.0f", momValues[i]), Form("Fit Adjusted %.0f GeV", momValues[i]), 800, 600);
        hDepsAdjusted[i]->SetTitle(Form("Energy Deposition at %.0f GeV (Adjusted);Energy (GeV);Counts", momValues[i]));
        hDepsAdjusted[i]->GetXaxis()->SetRangeUser(0, 0.2 * momValues[i]);
        hDepsAdjusted[i]->SetLineColor(kBlack);
        hDepsAdjusted[i]->Draw("HIST");
        fitsAdjusted[i]->SetLineColor(col8);  // Changed to col8 for distinction
        fitsAdjusted[i]->SetLineWidth(2);
        fitsAdjusted[i]->Draw("SAME");

        TLegend* legAdjusted = new TLegend(0.1, 0.7, 0.3, 0.9);
        legAdjusted->SetBorderSize(0);
        legAdjusted->SetFillStyle(0);
        legAdjusted->AddEntry(hDepsAdjusted[i], "Data", "l");
        legAdjusted->AddEntry(fitsAdjusted[i], "Crystal Ball Fit", "l");
        legAdjusted->Draw();

        cAdjusted->SaveAs(imagesDir + Form("fit_Adjusted_%.0f.png", momValues[i]));
        delete legAdjusted;
        delete cAdjusted;
    }

    // Beam test data
    //Double_t EBT[nMom] = {0.048, 0.043, 0.039, 0.036};
    Double_t EBT[nMom] = {0.0475, 0.0417, 0.0375, 0.035};
    Double_t eEBT[nMom] = {0.0, 0.0, 0.0, 0.0};
    Double_t eData[nMom] = {0.2, 0.2, 0.2, 0.2};

    // Calculate sigE_total to match EBT using Cited simulations
    std::cout << "\nCalculating sigE_total to match beam test data (EBT) using Cited simulations:\n";
    Double_t sigE_Cited[nMom] = {0.027, 0.025, 0.023, 0.021}; // Cited sigE values
    for (Int_t i = 0; i < nMom; i++) {
        Double_t ERes_intrinsic_sq = EResCited[i] * EResCited[i] - sigE_Cited[i] * sigE_Cited[i];
        Double_t sigE_total = TMath::Sqrt(EBT[i] * EBT[i] - ERes_intrinsic_sq);
        if (ERes_intrinsic_sq >= 0 && sigE_total >= 0 && !std::isnan(sigE_total) && !std::isinf(sigE_total)) {
            std::cout << "Momentum " << momValues[i] << " GeV: Required sigE_total = " << sigE_total
                      << " (Cited ERes = " << EResCited[i] << " at sigE_Cited = " << sigE_Cited[i] << ") to match EBT = " << EBT[i] << std::endl;
        } else {
            std::cout << "Momentum " << momValues[i] << " GeV: No valid sigE_total solution (EResCited^2 < sigE_Cited^2 or EBT^2 < intrinsic^2)" << std::endl;
        }
    }

    // Graphs for resolution
    auto gEResBT = new TGraphErrors(nMom, momValues, EBT, eData, eEBT);
    gEResBT->SetTitle("Beam Test");
    gEResBT->SetMarkerColor(col4);
    gEResBT->SetMarkerStyle(21);
    gEResBT->SetMarkerSize(3);

    auto gERes0 = new TGraphErrors(nMom, momValues, EResSigE0, eMomValues, eEResSigE0);
    gERes0->SetTitle("Simulation (sigE0)");
    gERes0->SetMarkerColor(col8);
    gERes0->SetMarkerStyle(23);
    gERes0->SetMarkerSize(2);

    auto gEResCited = new TGraphErrors(nMom, momValues, EResCited, eMomValues, eEResCited);
    gEResCited->SetTitle("Simulation (Cited)");
    //gEResCited->RemovePoint(1);
    //gEResCited->RemovePoint(2);
    //gEResCited->RemovePoint(3);
    gEResCited->SetMarkerColor(col10);
    gEResCited->SetMarkerStyle(22);
    gEResCited->SetMarkerSize(2);

    auto gEResAdjusted = new TGraphErrors(nMom, momValues, EResAdjusted, eMomValues, eEResAdjusted);
    gEResAdjusted->SetTitle("Simulation (Adjusted)");
    gEResAdjusted->SetMarkerColor(col6);  // Changed to col8 for distinction
    gEResAdjusted->SetMarkerStyle(20);    // Changed to triangle for distinction
    gEResAdjusted->SetMarkerSize(2);

    // Requirement region
    auto fReqMax = new TF1("ReqMax", res, 2.0, 12.0, 2);
    fReqMax->SetParameters(0.10, 0.03);

    auto fReqMin = new TF1("ReqMin", res, 2.0, 12.0, 2);
    fReqMin->SetParameters(0.10, 0.02);

    const Int_t nX = 110;
    Double_t fX[nX], fMax[nX], fMin[nX];
    for (int i = 0; i < nX; i++) {
        fX[i] = 2.0 + i * 0.1;
        fMax[i] = fReqMax->Eval(fX[i]);
        fMin[i] = fReqMin->Eval(fX[i]);
    }

    TGraph *grMax = new TGraph(nX, fX, fMax);
    TGraph *grMin = new TGraph(nX, fX, fMin);
    TGraph *grShade = new TGraph(2 * nX);
    for (int i = 0; i < nX; i++) {
        grShade->SetPoint(i, fX[i], fMax[i]);
        grShade->SetPoint(nX + i, fX[nX - i - 1], fMin[nX - i - 1]);
    }
    grShade->SetFillStyle(3013);
    grShade->SetFillColor(16);

    // Resolution plot
    TCanvas* c1 = new TCanvas("c1", "Energy Resolution", 1200, 800);
    gEResBT->SetTitle("Energy Resolution vs Momentum;Beam Energy (GeV);#sigma_{E}/E");
    gEResBT->GetYaxis()->SetRangeUser(0, 0.09);
    gEResBT->Draw("APE");
    grShade->Draw("F SAME");
    grMax->Draw("L SAME");
    grMin->Draw("L SAME");
    gERes0->Draw("PE SAME");
    gEResCited->Draw("PE SAME");
    gEResAdjusted->Draw("PE SAME");

    // Legend
    TLegend* leg = new TLegend(0.35, 0.60, 0.89, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(grShade, "EIC Requirement: 10%/#sqrt{E} #oplus (2-3)%", "f");
    leg->AddEntry(gEResBT, "Test Beam Data Including FTBF #frac{dp}{p}", "ep");
    leg->AddEntry(gEResAdjusted, "MC with Adjusted Beam #frac{dp}{p}", "ep");
    leg->AddEntry(gEResCited, "MC with Cited Beam #frac{dp}{p}", "ep");
    leg->AddEntry(gERes0, "MC with No Beam #frac{dp}{p}", "ep");
    leg->Draw();

    TLatex *lPrelim = new TLatex();
    lPrelim->SetNDC();
    lPrelim->SetTextSize(0.080);
    lPrelim->SetTextColorAlpha(kBlack, 0.35);
    lPrelim->DrawLatex(0.15, 0.40, "Preliminary Results, Baby BCAL");
    lPrelim->Draw();

    // Save resolution plot
    c1->SaveAs(imagesDir + "eResPlot.pdf");
    mcFile.Close();

    // Clean up dynamically allocated fits
    for (Int_t i = 0; i < nMom; i++) {
        delete fitsSigE0[i];
        delete fitsCited[i];
        delete fitsAdjusted[i];
    }
}

