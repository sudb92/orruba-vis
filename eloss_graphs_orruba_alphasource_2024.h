//Radii in mm, from GODDESS_NIM
const float mm = 1e-3; //in m

const float dEradius = 96.2*mm;
const float Eradius = 100.2*mm;

const float sx3_center_to_target_center = 40.8*mm;
const float sx3_strip_length = 75.0*mm;
const float target_to_sx3_0 = sx3_center_to_target_center - sx3_strip_length/2.0;
const float sx3width = 40.0*mm;

TRandom3 rnd(0); //for rnd.Rndm()

//Should stay untouched unless using a different kind of barrel detector/alpha source
//const float deadlayer_vert_eloss = 0.133;//deadlayer vertical energy loss in MeV
const float deadlayer_vert_eloss = 0.020;//deadlayer vertical energy loss in MeV
const float source_vert_eloss = 0.040;//source vertical energy loss in MeV

//Coordinate system assumed - centered at orruba geometric center, x-axis on beam right when flying down with the beam, y-axis straight upwards, +z-axis 'upstream' from target, i.e. beam moves along -z
//awkward assumption, apologies -  upstream barrel z=0 being 'close to target' and z=75mm being away from target was the goal

//Things to change/modify
float rotnewy = 2.5; //fine correction in degrees to target rotation (due to slack etc)
int boltPosition = 0; //Source's orientation direction as indicated by ladder bolt positions 0,1,2,3,4,5 for 0,60,120,180,240,300deg in phi
std::string orruba_layer = "Eu"; //Eu : E layer upstream, dEu: dE layer upstream, Eds: E layer downstream for FRIB2024 experiment

TVector3 geom_orruba(int clkpos=4, int stripnum=2, std::string layer="dEu", float z=37.5*mm) {
    assert(clkpos<12 && stripnum <4 && "clock positions 0-11 and strip numbers 0,1,2,3 only");
    assert((layer == "dEu" || layer == "Eu") && "layer 'dEu' or 'Eu' only for now");

    TVector3 target_position(0,0,0);
    TVector3 sx3_position;

    float vradius = layer=="dEu"?dEradius:Eradius;

    //int clkpos=4; // at 4pm position
    //int stripnum=2; // 0,1,2,3

    float yoffset, phi, radius;
    if(stripnum==1 || stripnum ==2)
        yoffset = sx3width*1.0/8.0;
    if(stripnum==0 || stripnum ==3)
        yoffset = sx3width*3.0/8.0;

    //SX3 fronts facing up(true for dEs) and pointing to you, read 0123 right to left, also the direction phi curls in if +z is along beam.
    //SX3 fronts facing down(true for Es) and pointing to you, read 3210 right to left

    //dE has backs facing down, E's have fronts facing down.

    yoffset = std::pow(-1,stripnum/2)*yoffset;     //0,1 +ve for E, 2,3 -ve for E
    yoffset = layer=="dEu"?(-1)*yoffset:yoffset;
    //reverse for dE, 0,1, are -ve, 2,3 +ve.

    z= 1.0*(target_to_sx3_0+z);//+ve z upstream
    phi = (TMath::Pi()/6.)*(3-clkpos); //awkward 90deg bend to line things up right

    float deltaphi =  std::atan(yoffset/vradius);
    phi+=deltaphi;
    radius = std::sqrt(vradius*vradius + yoffset*yoffset);

    //sx3_position.SetPerp(radius);
    //sx3_position.SetPhi(phi);
    sx3_position.SetXYZ(radius*cos(phi), radius*sin(phi), z);

    //std::cout << radius << " " << phi << " " << z << std::endl;
    //return sx3_position.Theta();
    return sx3_position;
}

void eloss_graphs_orruba_alphasource_2024() {
        TCanvas c("c1","c1",2000,0,800,600);
        c.cd();
        float ledge,redge;
        int clk,strip;

        TGraph g[12][4];
        TGraph g2[12][4];
        TGraph g3[12][4];

        float tht=27.0*TMath::Pi()/180.; //Fixed ladder angle, 27 deg
        float ph=(boltPosition*60+rotnewy)*TMath::Pi()/180.; //how much 'twist' the ladder goes through in the current bolt position, in deg

        TVector3 source_dir(0.0,TMath::Sin(tht), TMath::Cos(tht)); //area vector out of source is oriented along the normal to phi=0, theta=27.0 originally.
        TVector3 lad_dir(0.0, TMath::Cos(tht), -TMath::Sin(tht));
        source_dir.Rotate(ph,lad_dir); //rotate 3d vector 'source_dir' by angle 'ph' about the vector 'lad_dir'

        for(int clock=0; clock<12; clock++) {
          for(int strip=0; strip<4; strip++) {
            for(float z=0; z<sx3_strip_length; z+=0.2*mm) {

                //std::cout << z << "\t" << geom_orruba(clock,strip,"Eu",z)*180./TMath::Pi() << std::endl;
                //g[clock][strip].AddPoint(z,geom_orruba(clock,strip,"dEu",z)*180./TMath::Pi());

                //Get ORRUBA position vector as a TVector3, for the given (clock, strip, z) values. Assumes geometric center of the barrel.
                TVector3 vec = geom_orruba(clock,strip,orruba_layer,z); //"dEu" for dE layer upstream, "Eu" for E layer upstream
                if(z<0.2*mm) std::cout << strip << " " << vec.X() << " " << vec.Y() << " " << vec.Z() << std::endl;
                float theta = vec.Theta();

                float cos_sourceangle = vec.Dot(source_dir)/(vec.Mag()*source_dir.Mag());
                float sourceangle =  (180.0/TMath::Pi())*std::acos(cos_sourceangle); //angle between the source's normal vector and the barrel position vector
                float sourcethick = 1.0/std::abs(cos_sourceangle); //distance penetrated through the alpha source as a fraction of its normal thickness

                float thickness_seen = 1.0/std::sin(theta); //in um 
                //^^ 133keV eloss empirically observed by SD Pain by flying alphas vertically into backs/fronts and calculating difference. SRIM chalks this to 1um roughly of Si at 5.485 MeV

                float eloss = 5.486-deadlayer_vert_eloss*thickness_seen; //only deadlayer loss
                float eloss2 = 5.486-deadlayer_vert_eloss*thickness_seen-source_vert_eloss*sourcethick; //both source and deadlayer loss
                float eloss3 = 5.486-source_vert_eloss*sourcethick; //only sourceloss
                if(sourceangle>76.0) {
                    //ignore points for which the source's emission cone larger than 76deg. no evidence alphas leave this cone
                    eloss2=0.0;
                    eloss3=0.0;
                }
                g[clock][strip].AddPoint(z,eloss);
                g2[clock][strip].AddPoint(z,eloss2);
                g3[clock][strip].AddPoint(z,eloss3);
                //std::cout << cos_sourceangle << " " << eloss <<" "  << eloss2 << std::endl;

            } //z for
            g[clock][strip].SetTitle(Form("clock:%d, strip:%d, only dead layer loss",clock,strip));
            g[clock][strip].SetMarkerColor(kOrange-3);
            g2[clock][strip].SetMarkerColor(kRed-3);
            g3[clock][strip].SetMarkerColor(kGreen-3);
            g[clock][strip].SetLineColor(kOrange-3);
            g2[clock][strip].SetLineColor(kRed-3);
            g3[clock][strip].SetLineColor(kGreen-3);
            //g[clock][strip].SetName("Only dead-layer loss");
            g2[clock][strip].SetName("Source+dead-layer loss");
            g3[clock][strip].SetName("Only Source loss");

            g[clock][strip].SetMarkerStyle(kFullCircle);             g[clock][strip].SetMarkerSize(0.4);
            g2[clock][strip].SetMarkerStyle(kFullCircle);            g2[clock][strip].SetMarkerSize(0.4);
            g3[clock][strip].SetMarkerStyle(kFullCircle);            g3[clock][strip].SetMarkerSize(0.4);
            g[clock][strip].GetYaxis()->SetRangeUser(5.0,5.8);
            TLine L1(0.0,5.486,76.0*mm,5.486); L1.SetLineColor(kBlue-3); 
            g[clock][strip].Draw("AP");
            g2[clock][strip].Draw("P");
            g3[clock][strip].Draw("P");
            c.BuildLegend();
            L1.Draw("SAME");
            c.Update();
            while(c.WaitPrimitive());
         }//strip for
      }//clock for
}








////
////    Miscellaneous unused functions usable elsewhere
////
////



TVector3 geom_orruba_eds(int clkpos=4, int stripnum=2, std::string layer="Eds", float z=37.5*mm) {
    assert( (clkpos==5 || clkpos==17) && stripnum <4 && "clock positions 5 or 17 and strip numbers 0,1,2,3 only");
    assert((layer == "Eds") && "layer 'Eds' only for downstream");

    TVector3 target_position(0,0,0);
    TVector3 sx3_position;

    if(clkpos==5) clkpos = 4;
    if(clkpos==17) clkpos = 6;

    float vradius = Eradius;

    //int clkpos=4; // at 4pm position
    //int stripnum=2; // 0,1,2,3

    float yoffset, phi, radius;
    if(stripnum==1 || stripnum ==2)
        yoffset = sx3width*1.0/8.0;
    if(stripnum==0 || stripnum ==3)
        yoffset = sx3width*3.0/8.0;

    //SX3 fronts facing up(true for dEs) and pointing to you, read 0123 right to left, also the direction phi curls in if +z is along beam.
    //SX3 fronts facing down(true for Es) and pointing to you, read 3210 right to left
    //SX3 fronts downstream face down(true for Eds) and pointing towards target

    //dE has backs facing down, E's have fronts facing down.

    yoffset = std::pow(-1,stripnum/2)*yoffset;
    //yoffset = (-1)*yoffset;
    //offset correction for Eds, 0,1, are -ve, 2,3 +ve.

    z= -1.0*(target_to_sx3_0+z); //+ve sign because upstream sx3s are in +ve z axis 
    std::cout <<" --> " << clkpos << " " << stripnum << " " << z << std::endl;
    phi = (2*TMath::Pi()/12.)*((clkpos+4)%12);
    float deltaphi =  std::atan(yoffset/vradius);
    phi+=deltaphi;
    radius = std::sqrt(vradius*vradius + yoffset*yoffset);

    //sx3_position.SetPerp(radius);
    //sx3_position.SetPhi(phi);
    sx3_position.SetXYZ(radius*cos(phi), radius*sin(phi), z);

    //std::cout << radius << " " << phi << " " << z << std::endl;
    //return sx3_position.Theta();
    return sx3_position;
}

//Function to get 8 SX3 vertices given a clock position and a layer
std::array<float,24> vertices(int clkpos, std::string layer="E", float* opt=nullptr) {
    std::array<TVector3,8> result;
    float center_phi = 30*(3-clkpos)*TMath::Pi()/180.; //in rad
    float center_z = sx3_strip_length+target_to_sx3_0;
    float center_z_near = target_to_sx3_0;
    float center_r;
    if(layer=="dE")
        center_r = dEradius;
    else if (layer=="E")
        center_r = Eradius;
    else {//Eds
        center_r = Eradius;
        center_z *= -1.0;
        center_z_near*=-1.0;
    }
    TVector3 nn(TMath::Cos(center_phi),TMath::Sin(center_phi),0.); //'normal': radial vector out of cylinder
    TVector3 tt(-TMath::Sin(center_phi),TMath::Cos(center_phi),0.); //'tangent', will be along the detector's short edge
    TVector3 c(center_r*TMath::Cos(center_phi), center_r*TMath::Sin(center_phi), center_z); //center of edge farthest from source
    TVector3 cnear(center_r*TMath::Cos(center_phi), center_r*TMath::Sin(center_phi),center_z_near); //center of edge nearest to source

    result.at(0) = c+nn*0.5*mm-tt*0.5*sx3width;
    result.at(1) = c+nn*0.5*mm+tt*0.5*sx3width;
    result.at(2) = c-nn*0.5*mm+tt*0.5*sx3width;
    result.at(3) = c-nn*0.5*mm-tt*0.5*sx3width;

    result.at(4) = cnear+nn*0.5*mm-tt*0.5*sx3width;
    result.at(5) = cnear+nn*0.5*mm+tt*0.5*sx3width;
    result.at(6) = cnear-nn*0.5*mm+tt*0.5*sx3width;
    result.at(7) = cnear-nn*0.5*mm-tt*0.5*sx3width;

    std::array<float,24> result_c_arr;
    for(int i=0; i<8; i++) {
            result_c_arr.at(i*3+0) = result.at(i).x();
            result_c_arr.at(i*3+1) = result.at(i).y();
            result_c_arr.at(i*3+2) = result.at(i).z();
            if(opt) {
                opt[i*3+0] = result.at(i).x();
                opt[i*3+1] = result.at(i).y();
                opt[i*3+2] = result.at(i).z();
                //std::cout << "(" << opt[i+0] << " " << opt[i+1] << " " << opt[i+2] <<")" << std::endl;
            }
    }
    return result_c_arr;
}
