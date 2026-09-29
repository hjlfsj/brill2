{
//========= Macro generated from object: cal_d2_d3_14O_beam_cut/Graph
//========= by ROOT version6.38.04
   
   std::vector<Double_t> cutg_vect0{ 153.1250000465661, 173.076923420796, 239.6634628745513, 250.4807707280494, 173.076923420796, 153.3653846655328, 153.1250000465661 };
   std::vector<Double_t> cutg_vect1{ 303.1159264008477, 271.7751708515494, 221.3758477384887, 270.0810759569927, 318.3627804518576, 304.8100212954043, 303.1159264008477 };
   TCutG *cutg = new TCutG("cal_d2_d3_14O_beam_cut", 7, cutg_vect0.data(), cutg_vect1.data());
   cutg->SetVarX("PID D2-D3 (hit0)");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetLineColor(2);
   cutg->SetLineWidth(2);
   cutg->Draw();
}
