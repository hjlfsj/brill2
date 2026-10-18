{
//========= Macro generated from object: cal_d2_d3_tran_He4_cut/Graph
//========= by ROOT version6.38.04
   
   std::vector<Double_t> cutg_vect0{ 6.275323910527367, 11.15844086640837, 21.6519475162803, 34.63896069681488, 50.22337651345636, 53.23636357134038, 46.89870113923951, 28.61298658104683, 12.19740192085113, 7.002596648637306, 6.275323910527367, 6.275323910527367 };
   std::vector<Double_t> cutg_vect1{ 10.77678153492016, 7.297728164253414, 13.26181965682498, 18.14906129657113, 21.62811466723787, 24.36165660133317, 26.43252170292052, 23.53331056069823, 17.48638446406317, 13.0133158446345, 10.61111232679318, 10.77678153492016 };
   TCutG *cutg = new TCutG("cal_d2_d3_tran_He4_cut", 12, cutg_vect0.data(), cutg_vect1.data());
   cutg->SetVarX("PID D2-D3 (pre-cal hit0)");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetLineColor(2);
   cutg->SetLineWidth(2);
   cutg->Draw();
}
