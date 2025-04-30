

TH1F* hHits(Double_t Energy, Int_t nbins, Int_t ElePim, Int_t DepCal, Int_t cut, Int_t kB, Int_t sigE, Int_t pList, Int_t sigHV, Int_t sigMean, Int_t Containment){

    Double_t xLow;  Double_t xHigh;
    if (DepCal == 0) {
        xLow = -0.01*Energy;
        xHigh = 0.02*Energy;
    } else {
        xLow = 0;
        xHigh = 16;
    }

    TString ElePimString[2] = {"Ele", "Pim"};
    TString DepCalString[2] = {"Dep", "Cal"};
    TString cutString[4] = {"", "_APC", "_APC_T", "_APC_TC"};
    TString kBString[6] = {"_kB0132", "_kB033", "_kB0462", "_kB0528", "_kB066", "_kB132"};
    TString pListString[5] = {"_QGSP_BERT", "_FTFP_BERT", "_QGSP_BIC", "_QGSP_BERT_HP", "_FTFP_BERT_HP"};
    TString sigEString[6] = {"_sigE0", "_Cited", "_Adjusted"};
    TString sigHVString[10] = {"_9.4", "_50.3", "_91.2", "_22.0", "_29.4", "_36.8", "_18.7", "_20.6", "_22.5", "_17.9"};
    TString sigMeanString[10] = {"_33.7", "_37.3", "_40.9", "_2.6", "_14.3", "_26.0", "_4.4", "_7.3", "_10.2", "_7.8"};
    TString ContainmentString[3] = {"_O", "_U", "_u"};
    TH1F* h = new TH1F(Form("h%s%s_%.1fGeV%s%s%s%s%s%s%s",
                            ElePimString[ElePim].Data(),
                            DepCalString[DepCal].Data(),
                            Energy,
                            cutString[cut].Data(),
                            kBString[kB].Data(),
                            sigEString[sigE].Data(),
                            pListString[pList].Data(),
                            sigHVString[sigHV].Data(),
                            sigMeanString[sigMean].Data(),
                            ContainmentString[Containment].Data()),
                       Form("pim %.1f GeV, %s; E_{%s} [GeV]; bin count",
                            Energy, cutString[cut].Data(), DepCalString[DepCal].Data()),
                       nbins, xLow, xHigh);
    return h;
}


double SylvestersFormula (double *x, double *parSyl) {
    // [I0] * ([alpha] * Exp(-[x]/[lambda1]) + (1-[alpha]) * Exp(-[x]/[lambda2]))
    return parSyl[0] * (parSyl[1] * TMath::Exp(-x[0]/parSyl[2]) + (1-parSyl[1]) * TMath::Exp(-x[0]/parSyl[3]));
}


TH1F* Histos_DepCal_APC_TC (TChain* chain, Double_t mom, Int_t part, Int_t contain, Int_t pList, Int_t sigE, Int_t Birks, Int_t sigHV, Int_t sigMean) {

    Int_t i = 0;
    if (mom == 4) {i = 0;}
    else if (mom == 4.2) {i = 0;}
    else if (mom == 4.4) {i = 0;}
    else if (mom == 4.5) {i = 0;}
    else if (mom == 4.6) {i = 0;}
    else if (mom == 4.8) {i = 0;}
    else if (mom == 5) {i = 0;}
    else if (mom == 6) {i = 1;}
    else if (mom == 6.2) {i = 1;}
    else if (mom == 8) {i = 2;}
    else if (mom == 8.2) {i = 2;}
    else if (mom == 10) {i = 3;}
    else if (mom == 10.2) {i = 3;}

    const int nE= 4;
    int Ene[nE] = {4, 6, 8, 10};
    auto nbins = 3000;
    int r = 10000;  // [mm]
    float threshold = 0.0010;  // [GeV]

    Double_t ESF[4] = {0.0982306, 0.0978238, 0.0974763, 0.0971651};
    Double_t ECC[4] = {0.0957395, 0.0958678, 0.0958741, 0.0956800};
    Double_t lightGuideEffFactor[10] = {1.151630, 1.117192,
                                        1.084063, 1.053218,
                                        1.021398, 0.992446,
                                        0.801124, 0.775977,
                                        0.752811, 0.736212};

    //--------------------------------------------------------------------------
    // ScFi Attenuation
    //--------------------------------------------------------------------------
    // Random number genertor
    TRandom* rand = new TRandom();

    int npePreAtten = 1100;
    int length = 58;  // [cm]
    int nParSyl = 4;
    double parSyl[4] = {1.37130409e+02, 1.81201485e-01,
                        6.08666936e+01, 4.18209843e+02};
     
    // Current Function
    TF1 *yValuesCurrent = new TF1("yValuesCurrent", SylvestersFormula,
                                  0.0, length, nParSyl);
    yValuesCurrent->SetParameters(parSyl);
    yValuesCurrent->SetParNames("I0", "#alpha", "#lambda_{1}", "#lambda_{2}");
    double currentToNpeFactor = npePreAtten / yValuesCurrent->Eval(0);

    // Number of photoelectrons per GeV Function
    TF1 *npePerGeV = new TF1("Number of PhotoElectrons per GeV", SylvestersFormula,
                             0.0, length, nParSyl);
    npePerGeV->SetParameters(parSyl);
    npePerGeV->SetParameter(0, parSyl[0] * currentToNpeFactor);
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    

    //---------------------------------------------------------------------
    // Create Histograms
    //---------------------------------------------------------------------
    TH1F *hDepAPC_T = hHits(mom, nbins, part, 0, 2, Birks, sigE, pList, sigHV, sigMean, contain);
    TH1F *hDepAPC_TC = hHits(mom, nbins, part, 0, 3, Birks, sigE, pList, sigHV, sigMean, contain);
    TH1F *hCalAPC_T = hHits(mom, nbins, part, 1, 2, Birks, sigE, pList, sigHV, sigMean, contain);
    TH1F *hCalAPC_TC = hHits(mom, nbins, part, 1, 3, Birks, sigE, pList, sigHV, sigMean, contain);
    
    //---------------------------------------------------------------------
    // Loop through data-members and fill histograms
    //---------------------------------------------------------------------
    // Initialize reader
    TTreeReader tree_reader(chain);

    // Access whatever data-members you need
    TTreeReaderArray<float> array_e(tree_reader, "EcalBabyBarrelScFiHits.energy");
    TTreeReaderArray<float> array_hitX(tree_reader, "EcalBabyBarrelScFiHits.position.x");
    TTreeReaderArray<uint64_t> array_c(tree_reader, "EcalBabyBarrelScFiHits.cellID");
    TTreeReaderArray<int> array_gs(tree_reader, "MCParticles.generatorStatus");
    TTreeReaderArray<double> array_x(tree_reader, "MCParticles.endpoint.x");
    TTreeReaderArray<double> array_y(tree_reader, "MCParticles.endpoint.y");
    TTreeReaderArray<double> array_z(tree_reader, "MCParticles.endpoint.z");

    // Loop through data-members and fill histograms
    while(tree_reader.Next()) { // Loop over events
        // Project Endpoints to leading surface of BBCal
        auto endX = 0.0;
        auto endY = 0.0;
        auto endZ = 0.0;
        for (int k=0; k<array_gs.GetSize(); k++) {
            if (array_gs[k]==1) {
                endX = array_x[k];  //depth
                endY = array_y[k];  //height
                endZ = array_z[k];  //width
            }
        }
        auto projY = -r*endY/endX;  //negative to cancel -x
        auto projZ = -r*endZ/endX;  //negative to cancel -x
        
    // Layer and Grid with              LAYER                    
    // Array Indices        1---2------3--------4----------5------6
    //            \   G  1 |-| 0,0 |  0,1  |   0,2   |    0,3    |-| Row0
    //----BEAM-----\  R  2 |-| 1,0 |  1,1  |   1,2   |    1,3    |-| Row1
    //-------------/  I  3 |-| 2,0 |  2,1  |   2,2   |    2,3    |-| Row2
    //            /   D  4 |-| 3,0 |  3,1  |   3,2   |    3,3    |-| Row3
    //                      ---------------------------------------

        Double_t Esum = 0.0;
        Double_t SipmSum[4][10] = {0};
        Double_t EsumA_N = 0.0;
        Double_t EsumAP_N = 0.0;
        Double_t EsumAPC_N = 0.0;
        Double_t EsumA_S = 0.0;
        Double_t EsumAP_S = 0.0;
        Double_t EsumAPC_S = 0.0;
        Double_t EsumAPC_T = 0.0;
        Double_t EsumAPC_TC = 0.0;
        Double_t SumNpeA_N = 0.0;
        Double_t SumNpeAP_N = 0.0;
        Double_t SumNpeA_S = 0.0;
        Double_t SumNpeAP_S = 0.0;
        Double_t SipmSumA[4][10] = {0};
        Double_t SipmSumAP[4][10] = {0};
        Double_t SipmSumAPC[4][10] = {0};
        Double_t SipmSumA_N[4][10] = {0};
        Double_t SipmSumAP_N[4][10] = {0};
        Double_t SipmSumAPC_N[4][10] = {0};
        Double_t SipmSumA_S[4][10] = {0};
        Double_t SipmSumAP_S[4][10] = {0};
        Double_t SipmSumAPC_S[4][10] = {0};
        Double_t Sipm2ChannelSumA[4][4] = {0};
        Double_t Sipm2ChannelSumAP[4][4] = {0};
        Double_t Sipm2ChannelSumAPC[4][4] = {0};
        Double_t Sipm2ChannelSumAPC_S[4][4] = {0};
        Double_t Sipm2ChannelSumAPC_N[4][4] = {0};
        Double_t ChannelSumAPC_T[4][4] = {0};
        Double_t ChannelSumAPC_T_S[4][4] = {0};
        Double_t ChannelSumAPC_T_N[4][4] = {0};
        Double_t ChannelSumAPC_TC[4][4] = {0};
        Double_t ChannelSumAPC_TC_S[4][4] = {0};
        Double_t ChannelSumAPC_TC_N[4][4] = {0};

        // Sort ScFiHits.energy by Cell ID
        for (Int_t ii=0; ii<array_e.GetSize(); ii++){
            // Bit Masking to retrieve cellID
            uint64_t cellID = array_c[ii];
            uint64_t gridAND = 1023;  // 2e10 - 1
            uint64_t layerAND = 63;  // 2e6 - 1
            auto Grid = (cellID >> 24) & gridAND;
            auto Layer = (cellID >> 14) & layerAND;
            Esum += array_e[ii]; 
            Double_t e = array_e[ii];

            bool L2 = Layer == 2;
            bool L3 = Layer == 3;
            bool L4 = Layer == 4;
            bool L5 = Layer == 5;
            bool L6 = Layer == 6;
            bool L7 = Layer == 7;
            bool L8 = Layer == 8;
            bool L9 = Layer == 9;
            bool L10 = Layer == 10;
            bool L11 = Layer == 11;

            // Add deposited energy to appropriate cell and sipm
            for (int g=0; g<4; g++) {
                for (int l=0; l<10; l++) {
                    int grid = g+1; 
                    int layer = l+2; 
                    //if ((Grid==grid) && (Layer==layer) && (Grid!=0) && (Layer!=0)) {
                    if ((Grid==grid) && (Layer==layer)) {
                        if (L2) {
                            SipmSum[g][0] += e;
                        } else if (L3) {
                            SipmSum[g][1] += e;
                        } else if (L4) {
                            SipmSum[g][2] += e;
                        } else if (L5) {
                            SipmSum[g][3] += e;
                        } else if (L6) {
                            SipmSum[g][4] += e;
                        } else if (L7) {
                            SipmSum[g][5] += e;
                        } else if (L8) {
                            SipmSum[g][6] += e;
                        } else if (L9) {
                            SipmSum[g][7] += e;
                        } else if (L10) {
                            SipmSum[g][8] += e;
                        } else if (L11) {
                            SipmSum[g][9] += e;
                        }
                    } 
                }
            }
        }

        // Attenuation North
        double dist_N;
        double distParam_N;
        double attenuationFactor_N;
        double npeA_N;
        double eneA_N;
        double npeAP_N;
        double eneAP_N;
        double eneAC_N;
        double eneAPC_N;

        for (int g=0; g<4; g++) {
            for (int l=0; l<10; l++) {
                npePerGeV->SetParameter(0, parSyl[0] * currentToNpeFactor * lightGuideEffFactor[l]);

                dist_N = (length/2.0) - (projZ/10.0);  // [cm]
                distParam_N = npePerGeV->Eval(dist_N);  // npe per GeV at dist
                attenuationFactor_N = distParam_N / npePerGeV->Eval(0);

                npeA_N = std::round(distParam_N * SipmSum[g][l] / ESF[i]);  // npe at dist
                eneA_N = npeA_N / npePerGeV->Eval(0) * ESF[i];  // Energy at dist
                npeAP_N = rand->Poisson(npeA_N);  // npe at dist with Poisson smearing
                eneAP_N = npeAP_N / npePerGeV->Eval(0) * ESF[i];  // Energy at dist with Poisson smearing
                eneAC_N = eneA_N / attenuationFactor_N;  // Energy corrected for attenuation
                eneAPC_N = eneAP_N / attenuationFactor_N;  // Poisson smeared energy corrected for attenuation

                EsumA_N += eneA_N;
                EsumAP_N += eneAP_N;
                EsumAPC_N += eneAPC_N;
                SumNpeA_N += npeA_N;
                SumNpeAP_N += npeAP_N;
                SipmSumA_N[g][l] += eneA_N;
                SipmSumAP_N[g][l] += eneAP_N;
                SipmSumAPC_N[g][l] += eneAPC_N;
            }
        }
        
        // Attenuation South
        double dist_S;
        double distParam_S;
        double attenuationFactor_S;
        double npeA_S;
        double eneA_S;
        double npeAP_S;
        double eneAP_S;
        double eneAC_S;
        double eneAPC_S;

        for (int g=0; g<4; g++) {
            for (int l=0; l<10; l++) {
                npePerGeV->SetParameter(0, parSyl[0] * currentToNpeFactor * lightGuideEffFactor[l]);

                dist_S = (length/2.0) + (projZ/10.0);  // [cm]
                distParam_S = npePerGeV->Eval(dist_S);  // npe per GeV at dist
                attenuationFactor_S = distParam_S / npePerGeV->Eval(0);

                npeA_S = std::round(distParam_S * SipmSum[g][l] / ESF[i]);  // npe at dist
                eneA_S = npeA_S / npePerGeV->Eval(0) * ESF[i];  // Energy at dist
                npeAP_S = rand->Poisson(npeA_S);  // npe at dist with Poisson smearing
                eneAP_S = npeAP_S / npePerGeV->Eval(0) * ESF[i];  // Energy at dist with Poisson smearing
                eneAC_S = eneA_S / attenuationFactor_S;  // Energy corrected for attenuation
                eneAPC_S = eneAP_S / attenuationFactor_S;  // Poisson smeared energy corrected for attenuation

                EsumA_S += eneA_S;
                EsumAP_S += eneAP_S;
                EsumAPC_S += eneAPC_S;
                SumNpeA_S += npeA_S;
                SumNpeAP_S += npeAP_S;
                SipmSumA_S[g][l] += eneA_S;
                SipmSumAP_S[g][l] += eneAP_S;
                SipmSumAPC_S[g][l] += eneAPC_S;
            }
        }
    
        // Combine North and South sides
        for (int g=0; g<4; g++) {
            for (int l=0; l<10; l++) {
                SipmSumA[g][l] = TMath::Sqrt(SipmSumA_N[g][l]*SipmSumA_S[g][l]);
                SipmSumAP[g][l] = TMath::Sqrt(SipmSumAP_N[g][l]*SipmSumAP_S[g][l]);
                SipmSumAPC[g][l] = TMath::Sqrt(SipmSumAPC_N[g][l]*SipmSumAPC_S[g][l]);
            }
        }

        // Sum Sipms to Readout Channels
        for (int g=0; g<4; g++) {
            for (int l=0; l<10; l++) {
                if (l==0) {
                    Sipm2ChannelSumA[g][0] += SipmSumA[g][l];
                    Sipm2ChannelSumAP[g][0] += SipmSumAP[g][l];
                    Sipm2ChannelSumAPC[g][0] += SipmSumAPC[g][l];
                    Sipm2ChannelSumAPC_S[g][0] += SipmSumAPC_S[g][l];
                    Sipm2ChannelSumAPC_N[g][0] += SipmSumAPC_N[g][l];
                } else if (l==1) {
                    Sipm2ChannelSumA[g][1] += SipmSumA[g][l];
                    Sipm2ChannelSumAP[g][1] += SipmSumAP[g][l];
                    Sipm2ChannelSumAPC[g][1] += SipmSumAPC[g][l];
                    Sipm2ChannelSumAPC_S[g][1] += SipmSumAPC_S[g][l];
                    Sipm2ChannelSumAPC_N[g][1] += SipmSumAPC_N[g][l];
                } else if (l==2) {
                    Sipm2ChannelSumA[g][1] += SipmSumA[g][l];
                    Sipm2ChannelSumAP[g][1] += SipmSumAP[g][l];
                    Sipm2ChannelSumAPC[g][1] += SipmSumAPC[g][l];
                    Sipm2ChannelSumAPC_S[g][1] += SipmSumAPC_S[g][l];
                    Sipm2ChannelSumAPC_N[g][1] += SipmSumAPC_N[g][l];
                } else if (l==3) {
                    Sipm2ChannelSumA[g][2] += SipmSumA[g][l];
                    Sipm2ChannelSumAP[g][2] += SipmSumAP[g][l];
                    Sipm2ChannelSumAPC[g][2] += SipmSumAPC[g][l];
                    Sipm2ChannelSumAPC_S[g][2] += SipmSumAPC_S[g][l];
                    Sipm2ChannelSumAPC_N[g][2] += SipmSumAPC_N[g][l];
                } else if (l==4) {
                    Sipm2ChannelSumA[g][2] += SipmSumA[g][l];
                    Sipm2ChannelSumAP[g][2] += SipmSumAP[g][l];
                    Sipm2ChannelSumAPC[g][2] += SipmSumAPC[g][l];
                    Sipm2ChannelSumAPC_S[g][2] += SipmSumAPC_S[g][l];
                    Sipm2ChannelSumAPC_N[g][2] += SipmSumAPC_N[g][l];
                } else if (l==5) {
                    Sipm2ChannelSumA[g][2] += SipmSumA[g][l];
                    Sipm2ChannelSumAP[g][2] += SipmSumAP[g][l];
                    Sipm2ChannelSumAPC[g][2] += SipmSumAPC[g][l];
                    Sipm2ChannelSumAPC_S[g][2] += SipmSumAPC_S[g][l];
                    Sipm2ChannelSumAPC_N[g][2] += SipmSumAPC_N[g][l];
                } else if (l==6) {
                    Sipm2ChannelSumA[g][3] += SipmSumA[g][l];
                    Sipm2ChannelSumAP[g][3] += SipmSumAP[g][l];
                    Sipm2ChannelSumAPC[g][3] += SipmSumAPC[g][l];
                    Sipm2ChannelSumAPC_S[g][3] += SipmSumAPC_S[g][l];
                    Sipm2ChannelSumAPC_N[g][3] += SipmSumAPC_N[g][l];
                } else if (l==7) {
                    Sipm2ChannelSumA[g][3] += SipmSumA[g][l];
                    Sipm2ChannelSumAP[g][3] += SipmSumAP[g][l];
                    Sipm2ChannelSumAPC[g][3] += SipmSumAPC[g][l];
                    Sipm2ChannelSumAPC_S[g][3] += SipmSumAPC_S[g][l];
                    Sipm2ChannelSumAPC_N[g][3] += SipmSumAPC_N[g][l];
                } else if (l==8) {
                    Sipm2ChannelSumA[g][3] += SipmSumA[g][l];
                    Sipm2ChannelSumAP[g][3] += SipmSumAP[g][l];
                    Sipm2ChannelSumAPC[g][3] += SipmSumAPC[g][l];
                    Sipm2ChannelSumAPC_S[g][3] += SipmSumAPC_S[g][l];
                    Sipm2ChannelSumAPC_N[g][3] += SipmSumAPC_N[g][l];
                } else if (l==9) {
                    Sipm2ChannelSumA[g][3] += SipmSumA[g][l];
                    Sipm2ChannelSumAP[g][3] += SipmSumAP[g][l];
                    Sipm2ChannelSumAPC[g][3] += SipmSumAPC[g][l];
                    Sipm2ChannelSumAPC_S[g][3] += SipmSumAPC_S[g][l];
                    Sipm2ChannelSumAPC_N[g][3] += SipmSumAPC_N[g][l];
                }
            }
        }

        // Threshold Cut
        for (int g=0; g<4; g++) {
            for (int l=0; l<4; l++) {
                if (Sipm2ChannelSumAP[g][l] > threshold) {
                    ChannelSumAPC_T[g][l] = Sipm2ChannelSumAPC[g][l];
                    ChannelSumAPC_T_S[g][l] = Sipm2ChannelSumAPC_S[g][l];
                    ChannelSumAPC_T_N[g][l] = Sipm2ChannelSumAPC_N[g][l];
                    EsumAPC_T += ChannelSumAPC_T[g][l];
                }
            }
        }
        if (EsumAPC_T > 0.0) {
            hDepAPC_T->Fill(EsumAPC_T);
            hCalAPC_T->Fill(EsumAPC_T/ECC[i]);
        }
        
        // Containment Cut
        Double_t MaxCell = 0.0;
        int maxG = 0;
        int maxL = 0;
        for (int g=0; g<4; g++) {
            for (int l=0; l<4; l++) {
                if (ChannelSumAPC_T[g][l]>MaxCell) {
                    MaxCell = ChannelSumAPC_T[g][l];
                    maxG = g;
                    maxL = l;
                }
            }
        }

        if (contain == 0) {
            if ((maxG==1 || maxG==2) && (maxL==1 || maxL==2)) { // "O"
                for (int g=0; g<4; g++) {
                    for (int l=0; l<4; l++) {
                        ChannelSumAPC_TC[g][l] = ChannelSumAPC_T[g][l];
                        ChannelSumAPC_TC_S[g][l] = ChannelSumAPC_T_S[g][l];
                        ChannelSumAPC_TC_N[g][l] = ChannelSumAPC_T_N[g][l];
                        EsumAPC_TC += ChannelSumAPC_TC[g][l];
                    }
                }
                hDepAPC_TC->Fill(EsumAPC_TC);
                hCalAPC_TC->Fill(EsumAPC_TC/ECC[i]);
            }
        } else if (contain == 1) {
            if ((maxG==1 || maxG==2) && (maxL==1 || maxL==2 || maxL==3)) { // "U"
                for (int g=0; g<4; g++) {
                    for (int l=0; l<4; l++) {
                        ChannelSumAPC_TC[g][l] = ChannelSumAPC_T[g][l];
                        ChannelSumAPC_TC_S[g][l] = ChannelSumAPC_T_S[g][l];
                        ChannelSumAPC_TC_N[g][l] = ChannelSumAPC_T_N[g][l];
                        EsumAPC_TC += ChannelSumAPC_TC[g][l];
                    }
                }
                hDepAPC_TC->Fill(EsumAPC_TC);
                hCalAPC_TC->Fill(EsumAPC_TC/ECC[i]);
            }
        } else if (contain == 2) {
            if (((maxG==1 || maxG==2) && (maxL==1 || maxL==2)) || maxL==3) { // "u"
                for (int g=0; g<4; g++) {
                    for (int l=0; l<4; l++) {
                        ChannelSumAPC_TC[g][l] = ChannelSumAPC_T[g][l];
                        ChannelSumAPC_TC_S[g][l] = ChannelSumAPC_T_S[g][l];
                        ChannelSumAPC_TC_N[g][l] = ChannelSumAPC_T_N[g][l];
                        EsumAPC_TC += ChannelSumAPC_TC[g][l];
                    }
                }
                hDepAPC_TC->Fill(EsumAPC_TC);
                hCalAPC_TC->Fill(EsumAPC_TC/ECC[i]);
            }
        }

    }

    TString ElePimString[2] = {"Ele", "Pim"};
    std::cout << Form("Histogram Filled: %.1fGeV %s", mom, ElePimString[part].Data()) << std::endl;


    // Write Histos to file
    hDepAPC_T->Write();
    hDepAPC_TC->Write();
    hCalAPC_T->Write();
    hCalAPC_TC->Write();

    return 0;
}


void createHistos(TString simDir, TString simFileName, TString histoDir, TString histoFileName, Int_t part, Double_t mom, Int_t contain, Int_t sigE, Int_t pList, Int_t Birks, Int_t sigHV, Int_t sigMean){

    TFile outFile(histoDir + histoFileName, "UPDATE");

    auto treename = "events";
    TChain *chain = new TChain(treename);
    chain->Add(simDir + simFileName);
    Histos_DepCal_APC_TC (chain, mom, part, contain, pList, sigE, Birks, sigHV, sigMean);

    // Delete the pointer to chain
    delete chain;
}


void histos() {
    //std::cout << "paperHistos2.cpp(Int_t particle, Int_t momentum, Int_t centrality, Int_t sigmaE, Int_t physicsList, Int_t BirksConstant, Int_t sigHV, Int_t sigMean)" << std::endl;
    //std::cout << "Int_t particle = {0=ele, 1=pim}" << std::endl;
    //std::cout << "Int_t momentum = {4, 4.2, 4.4, 4.6, 4.8, 5, 6, 6.2, 8, 8.2, 10, 10.2}" << std::endl;
    //std::cout << "Int_t centrality = {0=O, 1=U, 2=u}" << std::endl;
    //std::cout << "Int_t sigmaE = {0=O, 1=Cited, 2=Adjusted}" << std::endl;
    //std::cout << "Int_t physicsList = {0=QGSP_BERT, 1=FTFP_BERT, 2=QGSP_BIC, 3=QGSP_BERT_HP, 4=FTFP_BERT_HP}" << std::endl;
    //std::cout << "Int_t BirksConstant = {0=kB0132, 1=kB033, 2=kB0462, 3=kB0528, 4=kB066, 5=kB132}" << std::endl;
    //std::cout << "Int_t sigHV = {0=9.4, 1=50.3, 2=91.2, 3=22.0, 4=29.4, 5=36.8, 6=18.7, 7=20.6, 8=22.5, 9=17.9};
    //std::cout << "Int_t sigMean = {0=33.7, 1=37.3, 2=40.9, 3=2.6, 4=14.3, 5=26.0, 6=4.4, 7=7.3, 8=10.2, 9=7.8};

    TString simDir = "sims/";
    TString simFileName = "MCSimulation.edm4hep.root";

    TString histoDir = "histograms/";
    TString histoFileName = "Histos.root";
    
    // Check if histoDir exists, create it if it doesn't
    if (gSystem->AccessPathName(histoDir)) {
        std::cout << "Directory " << histoDir << " does not exist. Creating it..." << std::endl;
        gSystem->mkdir(histoDir, kTRUE); // kTRUE for recursive creation
    }

    createHistos(simDir, simFileName, histoDir, histoFileName, 0, 4, 1, 2, 1, 5, 1, 0);
    
}

