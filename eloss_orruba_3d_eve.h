//Radii in mm, from GODDESS_NIM
const float mm = 1e-3; //in m

const float dEradius = 96.2*mm;
const float Eradius = 100.2*mm;

const float sx3_center_to_target_center = 40.8*mm;
const float sx3_strip_length = 75.0*mm;
const float target_to_sx3_0 = sx3_center_to_target_center - sx3_strip_length/2.0;
const float sx3width = 40.0*mm;

TRandom3 rnd(0); //for rnd.Rndm()

TVector3 geom_orruba(int clkpos=4, int stripnum=2, std::string layer="dEu", float z=37.5*mm) {
    assert(clkpos<12 && stripnum <4 && "clock positions 0-11 and strip numbers 0,1,2,3 only");
    assert((layer == "dEu" || layer == "Eu") && "layer 'dEu' or 'Eu' only for now");

    TVector3 target_position(0,0,0);
    TVector3 sx3_position;

    float vradius = layer=="dEu"?dEradius:Eradius;

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

    z= 1.0*(target_to_sx3_0+z); //+ve sign because upstream sx3s are in +ve z axis
    //phi = (2*TMath::Pi()/12.)*((clkpos+3)%12); //
    phi = (TMath::Pi()/6.)*(3-clkpos); //awkward 90deg bend to line things up right

    float deltaphi =  std::atan(yoffset/vradius);
    phi+=deltaphi;
    radius = std::sqrt(vradius*vradius + yoffset*yoffset);
    sx3_position.SetXYZ(radius*cos(phi), radius*sin(phi), z);

    return sx3_position;
}

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

    z= -1.0*(target_to_sx3_0+z); //-ve sign because upstream sx3s are in +ve z axis
    //std::cout <<" --> " << clkpos << " " << stripnum << " " << z << std::endl;
//    phi = (2*TMath::Pi()/12.)*((clkpos+3)%12);
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

//Function to get 8 SX3 vertices given a clock position and a layer
std::array<float,24> vertices(int clkpos, std::string layer="E", float* opt=nullptr) {
    std::array<TVector3,8> result;
    //float center_phi = 30*(3-clkpos)*TMath::Pi()/180.; //in rad
    float center_phi = (TMath::Pi()/6.)*(3-clkpos); //awkward 90deg bend to line things up right

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

void eloss_orruba_3d_eve() {
    gSystem->IgnoreSignal(kSigSegmentationViolation, true);
    TEveManager::Create();

    std::cout << gEve->GetWindowManager()->IsA()->GetName() << std::endl;
    TEveViewer *ev = gEve->GetDefaultViewer();
    TEveParamList *x = new TEveParamList("Click to rotate");
    gEve->AddToListTree(x,0);
    x->AddParameter(TEveParamList::BoolConfig_t("Cick to rotate",0));

    auto ccc = new TCanvas("c1","c1",0,0,80,80);
    ccc->SetFrameFillColor(kBlue-5);
    ccc->Update();
    while(gPad->WaitPrimitive());

    TGLViewer *gv = ev->GetGLViewer();
    gv->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);
    //gv->SetCurrentCamera(TGLViewer::kCameraPerspYOZ);

    auto source0 = new TGeoTube(0, 5*mm, 0.5*mm); //0 inner radius, 5 outer radius, 0.5 thickness
    TEveGeoShape *source1 = new TEveGeoShape("source"); //Convert to Eve friendly datatype
    source1->SetShape(source0);

    //Target ladder is angled at 27 degrees from +Z direction. The below gives (theta,phi)_i for i=x,y,z in order
    TGeoRotation mat("rotate_to_ladder",90.,0.,117,90.,27.0,90.0);
    source1->SetTransMatrix(mat);

    //Draw ORRUBA
    auto orruba = new TEveBoxSet("orruba","ORRUBA");
    orruba->UseSingleColor();
    orruba->SetMainColor(kCyan-2);
    orruba->SetMainTransparency(20);
    orruba->Reset(TEveBoxSet::kBT_FreeBox, kFALSE, 1);

    for(int i=0; i<12; i++) {
        if(i==6) continue;
        if(i==8) continue;
        if(i==3) continue;
        //Upstream barrel
        std::array<float,24> vertexlist_E = vertices(i,"E");
        std::array<float,24> vertexlist_dE = vertices(i,"dE");
        orruba->AddBox(vertexlist_E.data());
        orruba->AddBox(vertexlist_dE.data());
    }
    //Downstream detectors
    std::array<float,24> vertexlist_Eds = vertices(5,"Eds");
    orruba->AddBox(vertexlist_Eds.data());

    vertexlist_Eds = vertices(7,"Eds");
    orruba->AddBox(vertexlist_Eds.data());
    /////////


    //Source vector, and its normal along ladder direction
    TVector3 sourcevec;
    sourcevec.SetMagThetaPhi(500.0*mm,0.0,0.0);
    TVector3 ladderdir;
    ladderdir.SetXYZ(0,TMath::Sin(63*TMath::Pi()/180.),-TMath::Cos(63*TMath::Pi()/180.));
    ladderdir*=200.0*mm;//long long ladder

    //Cone to guesstimate the alpha illumination pattern
    auto alphacone = new TEveBoxSet("acone","a-cone");
    alphacone->UseSingleColor();
    alphacone->SetMainColor(kYellow);
    alphacone->SetMainTransparency(0);
//    alphacone->Reset(TEveBoxSet::kBT_EllipticCone,0,1);
    alphacone->Reset(TEveBoxSet::kBT_Cone,0,1);

    //The above are TVectors, we need Eve things to initizliae a cone
    TEveVector origin;
    origin.Set(TVector3(0,0,0));
    TEveVector sourcevec_eve;
    sourcevec_eve.Set(sourcevec);
    sourcevec_eve*=100.0;
    sourcevec_eve*=mm;
    alphacone->SetTransMatrix(mat); //transform the cone to point in the right direction
//    alphacone->AddEllipticCone(origin, sourcevec_eve, 150.0*mm, 150.0*mm, 160.0*TMath::Pi()/180.); //origin, direction, elliptical radii 1 and 2, angle
    alphacone->AddCone(origin, sourcevec_eve, 100.0*mm*TMath::Tan(63.0*TMath::Pi()/180.)); //origin, direction+magnitude of axis, radius

    //Arrow along source normal
    TEveArrow *a1 = new TEveArrow(sourcevec.x(), sourcevec.y(), sourcevec.z(), 0., 0., 0.); //vector x,y,z; then origin-point-of-vector x,y,z (0,0,0)
    a1->SetMainColor(kOrange-3);
    a1->SetTubeR(1*mm);
    a1->SetPickable(kTRUE);
    a1->SetTransMatrix(mat);

    //Arrow along target ladder, the axis of rotation
    auto a2 = new TEveArrow(ladderdir.x(), ladderdir.y(), ladderdir.z(), 0.,0.,0.);
    a2->SetMainColor(kCyan);
    a2->SetTubeR(1*mm);
    a2->SetPickable(kTRUE);

    //A few test arrows
    TVector3 vec0 = geom_orruba_eds(5,0,"Eds",75.*mm), vec1, vec2=geom_orruba_eds(5,0,"Eds",0.);
    std::cout << "Vec0 theta (deg): " << 180-vec0.Theta()*180./M_PI << " Vec0 phi (deg):" <<  vec0.Phi()*180./M_PI << std::endl;
    std::cout << "Vec2 theta (deg): " << 180-vec2.Theta()*180./M_PI << " Vec2 phi (deg):" <<  vec2.Phi()*180./M_PI << std::endl;
    vec1.SetMagThetaPhi(0.2,130.*M_PI/180.,240.*M_PI/180.);
    //vec0.SetXYZ(TMath::Sin(51*TMath::Pi()/180.)*TMath::Cos(120*M_PI/180.),TMath::Sin(51.*TMath::Pi()/180.)*TMath::Sin(120*M_PI/180.),-TMath::Cos(51*TMath::Pi()/180.));
    auto downstreamArrow_0 = new TEveArrow(vec0.x(),vec0.y(),vec0.z(),0.,0.,0.);
    downstreamArrow_0->SetMainColor(kMagenta-3);
    downstreamArrow_0->SetTubeR(1*mm);

    auto downstreamArrow_1 = new TEveArrow(vec1.x(),vec1.y(),vec1.z(),0.,0.,0.);
    downstreamArrow_1->SetMainColor(kMagenta-3);
    downstreamArrow_1->SetTubeR(1*mm);

    auto downstreamArrow_2 = new TEveArrow(vec2.x(),vec2.y(),vec2.z(),0.,0.,0.);
    downstreamArrow_2->SetMainColor(kMagenta-3);
    downstreamArrow_2->SetTubeR(1*mm);

    auto zeropos_arrow = new TEveArrow();
    zeropos_arrow->SetMainColor(kMagenta);
    zeropos_arrow->SetTubeR(2*mm);

    int ctr2=0;
    float rotnewy=2.5;

    //Add everything to viewer, render
    gEve->AddElement(a1);
    gEve->AddElement(a2);
//    gEve->AddElement(zeropos_arrow);
    gEve->AddElement(orruba);
    gEve->AddElement(alphacone);
    gEve->AddElement(source1);
    gEve->AddElement(downstreamArrow_0);
    gEve->AddElement(downstreamArrow_1);
    gEve->AddElement(downstreamArrow_2);
    gEve->FullRedraw3D(kTRUE,0);

    int count=0;
    while(1) {
        TGeoRotation mat3;
        mat3.RotateY(count*60+rotnewy); //new Y after 'mat' is along ladder
        source1->SetTransMatrix(mat*mat3);
        a1->SetTransMatrix(mat*mat3);
        alphacone->SetTransMatrix(mat*mat3);
        gEve->FullRedraw3D(kFALSE);

        //usleep(1e5);
        //gSystem->ProcessEvents();
        while(gPad->WaitPrimitive());
        count+=1;
        if(count==6*5) count=0;
    }
}
