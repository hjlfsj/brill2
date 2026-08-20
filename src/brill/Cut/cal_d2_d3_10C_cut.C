{
//========= Macro generated from object: cal_d2_d3_10C_cut/Graph
//========= by ROOT version6.38.04
   
   std::vector<Double_t> cutg_vect0{ 128.2433581978345, 131.4813299606135, 179.6911317619895, 238.6941727726287, 271.4336650407273, 267.116369357022, 182.2095542441509, 133.2802031621574, 128.2433581978345, 128.2433581978345 };
   std::vector<Double_t> cutg_vect1{ 112.0663908947444, 96.29703769827121, 106.3320806414814, 106.8099398292533, 115.4114052091478, 137.3929278466559, 130.7028992178491, 121.6235746501827, 112.0663908947444, 112.0663908947444 };
   TCutG *cutg = new TCutG("cal_d2_d3_10C_cut", 10, cutg_vect0.data(), cutg_vect1.data());
   cutg->SetVarX("e3");
   cutg->SetVarY("e2");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetLineColor(2);
   cutg->SetLineWidth(2);
   cutg->Draw();
}
