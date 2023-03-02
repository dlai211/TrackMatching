
void extract_align_fix3layers_combine(int iter=0){
  TString outputpath="/eos/user/k/keli/Faser/alignment/data/8023_8025_8115_newfieldmap/";
  TString outputname="data_iter";
  outputname+=iter;
  //do not align the first and last layer
  double align_sta[4][6]={0.};
  double align_lay[12][6]={0.};
  double align_mod[12][8][6]={0.};
  int ntrk_mod[12][8]={0};
  int ntrk_mod1[12][8]={0};
  int ntrk_sta[4]={0};
  int ntrk_lay[12]={0};
  //define matrix
  int Ndim=6;
  std::vector<TMatrixD> denom_lay;
  std::vector<TMatrixD> num_lay;
  std::vector<std::vector<TMatrixD>> denom_mod;
  std::vector<std::vector<TMatrixD>> num_mod;
  std::vector<std::vector<TMatrixD>> denom_mod1;
  std::vector<std::vector<TMatrixD>> num_mod1;
  //initialize to 0
  auto num_matrix = [](int n) {
    TMatrixD num_t(n,1);
    for(int i=0;i<n;i++){
      num_t[i][0]=0.;
    }
    return num_t;
  };
  auto denom_matrix = [](int n) {
    TMatrixD denom_t(n,n);
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
	denom_t[i][j]=0.;
      }
    }
    return denom_t;
  };
  std::vector<TMatrixD> denom_sta;
  std::vector<TMatrixD> num_sta;
  for(int i=0;i<4;i++){
    denom_sta.push_back(denom_matrix(Ndim));
    num_sta.push_back(num_matrix(Ndim));
  }
  for(int i=0;i<12;i++){
    denom_lay.push_back(denom_matrix(Ndim));
    num_lay.push_back(num_matrix(Ndim));
    std::vector<TMatrixD> denom_l;
    std::vector<TMatrixD> num_l;
    std::vector<TMatrixD> denom_l1;
    std::vector<TMatrixD> num_l1;
    for(int j=0;j<8;j++){
      denom_l.push_back(denom_matrix(Ndim));
      num_l.push_back(num_matrix(Ndim));
      denom_l1.push_back(denom_matrix(Ndim));
      num_l1.push_back(num_matrix(Ndim));
    }
    denom_mod.push_back(denom_l1);
    denom_mod1.push_back(denom_l1);
    num_mod.push_back(num_l);
    num_mod1.push_back(num_l1);
  }
  double delta_sta[6]={0.};
  std::ifstream input;
  int id,inum,ntrk;
  double m0,m1,m2,m3,m4,m5;
  std::cout<<"start reading input matrix"<<std::endl;
  for(int i=0;i<94;i++){
    input.open(Form(outputpath+"%d/alignment_matrix.txt",i));
    std::cout<<"open "<<Form(outputpath+"%d/alignment_matrix.txt",i)<<std::endl;
    bool iflay=false;
    bool ifmod=false;
    while(true){
      input>>id>>inum>>m0>>m1>>m2>>m3>>m4>>m5>>ntrk;
      if(input.eof())break;
      if(!iflay){
	ntrk_sta[id]+=ntrk;
	if(inum<6){
	  denom_sta[id][inum][0]+=m0;
	  denom_sta[id][inum][1]+=m1;
	  denom_sta[id][inum][2]+=m2;
	  denom_sta[id][inum][3]+=m3;
	  denom_sta[id][inum][4]+=m4;
	  denom_sta[id][inum][5]+=m5;
	}
	else{
	  std::cout<<i<<" "<<id<<" "<<inum<<" "<<m0<<" "<<num_sta[id][0][0]<<std::endl;
	  num_sta[id][0][0]+=m0;
	  num_sta[id][1][0]+=m1;
	  num_sta[id][2][0]+=m2;
	  num_sta[id][3][0]+=m3;
	  num_sta[id][4][0]+=m4;
	  num_sta[id][5][0]+=m5;
	}
      }
      if(iflay&&!ifmod){
	int ilay=(id/10)*3+id%10;
	ntrk_lay[ilay]+=ntrk;
	if(inum<6){
	  denom_lay[ilay][inum][0]+=m0;
	  denom_lay[ilay][inum][1]+=m1;
	  denom_lay[ilay][inum][2]+=m2;
	  denom_lay[ilay][inum][3]+=m3;
	  denom_lay[ilay][inum][4]+=m4;
	  denom_lay[ilay][inum][5]+=m5;
	}
	else{
	  num_lay[ilay][0][0]+=m0;
	  num_lay[ilay][1][0]+=m1;
	  num_lay[ilay][2][0]+=m2;
	  num_lay[ilay][3][0]+=m3;
	  num_lay[ilay][4][0]+=m4;
	  num_lay[ilay][5][0]+=m5;
	}
      }
      if(ifmod){
	int ista=id/100;
	int ilay=id/10%10+3*ista;
	int imod=id%10;
	ntrk_mod[ilay][imod]+=ntrk;
	if(inum<6){
	  denom_mod[ilay][imod][inum][0]+=m0;
	  denom_mod[ilay][imod][inum][1]+=m1;
	  denom_mod[ilay][imod][inum][2]+=m2;
	  denom_mod[ilay][imod][inum][3]+=m3;
	  denom_mod[ilay][imod][inum][4]+=m4;
	  denom_mod[ilay][imod][inum][5]+=m5;
	}
	else{
	  num_mod[ilay][imod][0][0]+=m0;
	  num_mod[ilay][imod][1][0]+=m1;
	  num_mod[ilay][imod][2][0]+=m2;
	  num_mod[ilay][imod][3][0]+=m3;
	  num_mod[ilay][imod][4][0]+=m4;
	  num_mod[ilay][imod][5][0]+=m5;
	}
      }
      if(id==3&&inum==6)iflay=true;
      if(id==32&&inum==6)ifmod=true;

    }
    input.close();
  }

  std::cout<<"Get the nominal transform"<<std::endl;
  //get the alignment constants
  std::ofstream align_output("all_alignment_"+outputname+".txt",ios::out);
  std::ofstream alignparam_output("all_alignment_input_"+outputname+".txt",ios::out);
  double inte=1.;
  for(int i=0;i<4;i++){
    if(ntrk_sta[i]>1000){
      inte=1;
      TMatrixD inv_sta=denom_sta[i].Invert(&inte);
      TMatrixD alignsta=inv_sta*num_sta[i];
      align_sta[i][0]=alignsta[0][0];
      align_sta[i][1]=alignsta[1][0];
      align_sta[i][2]=alignsta[2][0];
      align_sta[i][3]=alignsta[3][0];
      align_sta[i][4]=alignsta[4][0];
      align_sta[i][5]=alignsta[5][0];
      align_output<<i<<" "<<align_sta[i][0]<<"+-"<<sqrt(inv_sta[0][0])<<" "<<align_sta[i][1]<<"+-"<<sqrt(inv_sta[1][1])<<" "<<align_sta[i][2]<<"+-"<<sqrt(inv_sta[2][2])<<" "<<align_sta[i][3]<<"+-"<<sqrt(inv_sta[3][3])<<" "<<align_sta[i][4]<<"+-"<<sqrt(inv_sta[4][4])<<" "<<align_sta[i][5]<<"+-"<<sqrt(inv_sta[5][5])<<" "<<ntrk_sta[i]<<std::endl;
      alignparam_output<<i<<" "<<align_sta[i][0]<<" "<<align_sta[i][1]<<" "<<align_sta[i][2]<<" "<<align_sta[i][3]<<" "<<align_sta[i][4]<<" "<<align_sta[i][5]<<std::endl;
      if(i==0){
	delta_sta[0]=align_sta[i][0];
	delta_sta[1]=align_sta[i][1];
	delta_sta[2]=align_sta[i][2];
	delta_sta[3]=align_sta[i][3];
	delta_sta[4]=align_sta[i][4];
	delta_sta[5]=align_sta[i][5];
      }
    }
    else{
      alignparam_output<<i<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<std::endl;
      //	ROOT::Math::Translation3D translation(0,0,0);
      //	ROOT::Math::Transform3D trans(translation);
      //	delta_trans_sta[i]=trans;
    }
  }
  //layers
  for(int i=0;i<denom_lay.size();i++){
    std::cout<<"like get layer level "<<i<<std::endl;
    //      if(i==0||i==11)continue;
    inte=1.;
    if(ntrk_lay[i]<500){
      alignparam_output<<(i)/3<<(i)%3<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<std::endl;
    }
    else{
      TMatrixD inv_lay=denom_lay[i].Invert(&inte);
      TMatrixD alignlay=inv_lay*num_lay[i];
      align_lay[i][0]=alignlay[0][0];
      align_lay[i][1]=alignlay[1][0];
      align_lay[i][2]=alignlay[2][0];
      align_lay[i][3]=alignlay[3][0];
      align_lay[i][4]=alignlay[4][0];
      align_lay[i][5]=alignlay[5][0];
      align_output<<(i)/3<<(i)%3<<" "<<align_lay[i][0]<<"+-"<<sqrt(inv_lay[0][0])<<" "<<align_lay[i][1]<<"+-"<<sqrt(inv_lay[1][1])<<" "<<align_lay[i][2]<<"+-"<<sqrt(inv_lay[2][2])<<" "<<align_lay[i][3]<<"+-"<<sqrt(inv_lay[3][3])<<" "<<align_lay[i][4]<<"+-"<<sqrt(inv_lay[4][4])<<" "<<align_lay[i][5]<<"+-"<<sqrt(inv_lay[5][5])<<" "<<ntrk_lay[i]<<std::endl;
      if(i/3>0)
	alignparam_output<<(i)/3<<(i)%3<<" "<<align_lay[i][0]<<" "<<align_lay[i][1]<<" "<< align_lay[i][2]<<" "<< align_lay[i][3]<<" "<< align_lay[i][4]<<" "<<align_lay[i][5]<<std::endl;
      else
	alignparam_output<<(i)/3<<(i)%3<<" "<<align_lay[i][0]-delta_sta[0]<<" "<<align_lay[i][1]-delta_sta[1]<<" "<<align_lay[i][2]-delta_sta[2]<<" "<<align_lay[i][3]-delta_sta[3]<<" "<<align_lay[i][4]-delta_sta[4]<<" "<<align_lay[i][5]-delta_sta[5]<<std::endl;
    }
  }

  std::cout<<"Get the deltas for layers"<<std::endl;
  for(int i=0;i<denom_lay.size();i++){
    //      if(i==0||i==11)continue;
    for(  int j=0;j<num_mod[i].size();j++){
      inte=1.;
      if(ntrk_mod[i][j]<500){
	alignparam_output<<(i)/3<<(i)%3<<j<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<std::endl;
      }
      else{
	TMatrixD inv_mod=denom_mod[i][j].Invert(&inte);
	TMatrixD alignmod=inv_mod*num_mod[i][j];
	align_mod[i][j][0]=alignmod[0][0];
	align_mod[i][j][1]=alignmod[1][0];
	align_mod[i][j][2]=alignmod[2][0];
	align_mod[i][j][3]=alignmod[3][0];
	align_mod[i][j][4]=alignmod[4][0];
	align_mod[i][j][5]=alignmod[5][0];
	align_output<<(i)/3<<(i)%3<<j<<" "<<align_mod[i][j][0]<<"+-"<<sqrt(inv_mod[0][0])<<" "<<align_mod[i][j][1]<<"+-"<<sqrt(inv_mod[1][1])<<" "<<align_mod[i][j][2]<<"+-"<<sqrt(inv_mod[2][2])<<" "<<align_mod[i][j][3]<<"+-"<<sqrt(inv_mod[3][3])<<" "<<align_mod[i][j][4]<<"+-"<<sqrt(inv_mod[4][4])<<" "<<align_mod[i][j][5]<<"+-"<<sqrt(inv_mod[5][5])<<" "<<ntrk_mod[i][j]<<std::endl;
	alignparam_output<<(i)/3<<(i)%3<<j<<" "<<align_mod[i][j][0]<<" "<<align_mod[i][j][1]<<" "<<align_mod[i][j][2]<<" "<<align_mod[i][j][3]<<" "<<align_mod[i][j][4]<<" "<<align_mod[i][j][5]<<std::endl;
      }
    }
  }
  std::cout<<"closing the output"<<std::endl;
  alignparam_output.close();
  align_output.close();
}
