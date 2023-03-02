void makenewdb_fix3layers(int iter=0){
  std::ifstream input;
  int Niter=20;
  bool alignmod=true;
  bool alignift=false;
  if(iter%2==0) alignmod=false;
  //iterations
  //0 for the fixed 3 layers
  //1-8 for others
  //9-10 for the fixed 3 layers
  TString inputname="misalign_mc_";
  if(iter==0||iter<0)
    inputname+="ini.txt";
  else{
    inputname+="iter";
    int old = iter -1;
    inputname+=old;
    inputname+=".txt";
  }

//  int level=(iter-1)%3;//0, station, 1, layer, 2, module
  //fix two modules
  //int fix1=101;
  //int fix2=326;
  int fix=1;//middle layer in each station

  input.open(inputname,ios::in);
  TString corrname="all_alignment_input_mc_iter";
  corrname+=iter;
  corrname+=".txt";
  std::ifstream newcorr;
  newcorr.open(corrname,ios::in);
  TString outputname="misalign_mc_iter";
  outputname += iter;
  outputname += ".txt";
  std::cout<<"output file name "<<outputname<<std::endl;
  std::ofstream output(outputname, ios::out);
  std::ofstream outputalign("inputforalign.txt", ios::out);
  int id;
  double x,y,z,rx,ry,rz;
  int idnew;
  double xnew,ynew,znew,rxnew,rynew,rznew;
  int iline=0;
  bool iflay=false;
  bool ifmod=false;
  while(true){
    input>>id>>x>>y>>z>>rx>>ry>>rz;
    std::cout<<"like "<<x<<std::endl;
    if(iter>=0)
      newcorr>>idnew>>xnew>>ynew>>znew>>rxnew>>rynew>>rznew;
    if(input.eof())break;
    if(iter>3&&iter<Niter-4&&!alignift){
	if(id<9&&(!iflay)){
	  if(id==0&&!alignmod){
//	    x+=xnew;
//	    y+=ynew;
//	    z+=znew;
//	    rx+=rxnew;
//	    ry+=rynew;
//	    rz+=rznew;
	  }
	}
	if(iflay&&!ifmod&&!alignmod){
	  if(id%10!=fix&&id>9){
//	    x+=xnew;
	    y+=ynew;
//	    z+=znew;
//	    rx+=rxnew;
//	    ry+=rynew;
	    rz+=rznew;
	  }
      }
	if(ifmod&&alignmod){
	  if((id/10)%10!=fix&&id>99){
	    std::cout<<id<<" add correction "<<x<<" "<<xnew<<std::endl;
	    x+=xnew;
	    //	y+=ynew;
	    rz+=rznew;
	  }
	}
    }
    else if(iter<0&&!alignift){
      //do nothing, use the initial parametesr
      x+=0.;
      y+=0.;
      z+=0.;
      rx+=0.;
      ry+=0.;
      rz+=0.;
    }
    else if(alignift){
	if(id<9&&(!iflay)){
	  if(id==0&&!alignmod){
//	    x+=xnew;
	    y+=ynew;
//	    z+=znew;
//	    rx+=rxnew;
//	    ry+=rynew;
	    rz+=rznew;
	  }
	}
	if(iflay&&!ifmod&&!alignmod){
	  if(id<9){
//	    x+=xnew;
	    y+=ynew;
//	    z+=znew;
//	    rx+=rxnew;
//	    ry+=rynew;
	    rz+=rznew;
	  }
      }
	if(ifmod&&alignmod){
	  if(id<99){
	    std::cout<<id<<" add correction "<<x<<" "<<xnew<<std::endl;
	    x+=xnew;
	    //	y+=ynew;
	    rz+=rznew;
	  }
	}
    }
    else{
      //align the modules at the fixed layers
      if(iflay&&!ifmod&&id%10==fix&&id/10>0&&!alignmod){
//	x+=xnew;
	y+=ynew;
//	z+=znew;
//	rx+=rxnew;
//	ry+=rynew;
	rz+=rznew;
      }
      if(ifmod&&((id/10)%10==fix&&id/100>0)&&alignmod){
	std::cout<<id<<" add correction "<<x<<" "<<xnew<<std::endl;
	x+=xnew;
	rz+=rznew;
      }
    }
    //    output<<id<<" "<<x<<" "<<y<<" "<<z<<" "<<rx<<" "<<ry<<" "<<rz<<std::endl;
    if(iline==0){
      outputalign<<"\""<<id<<"\": ["<<x<<", "<<y<<", "<<z<<", "<<rx<<", "<<ry<<", "<<rz<<"]";
      output<<id<<" "<<x<<" "<<y<<" "<<z<<" "<<rx<<" "<<ry<<" "<<rz<<std::endl;
    }
    else{
      if(iflay&&!ifmod&&id<10){
	outputalign<<", \"0"<<id<<"\": ["<<x<<", "<<y<<", "<<z<<", "<<rx<<", "<<ry<<", "<<rz<<"]";
	output<<"0"<<id<<" "<<x<<" "<<y<<" "<<z<<" "<<rx<<" "<<ry<<" "<<rz<<std::endl;
      }
      else if(ifmod&&id<10){
	outputalign<<", \"00"<<id<<"\": ["<<x<<", "<<y<<", "<<z<<", "<<rx<<", "<<ry<<", "<<rz<<"]";
	output<<"00"<<id<<" "<<x<<" "<<y<<" "<<z<<" "<<rx<<" "<<ry<<" "<<rz<<std::endl;
      }
      else if(ifmod&&id<100){
	outputalign<<", \"0"<<id<<"\": ["<<x<<", "<<y<<", "<<z<<", "<<rx<<", "<<ry<<", "<<rz<<"]";
	output<<"0"<<id<<" "<<x<<" "<<y<<" "<<z<<" "<<rx<<" "<<ry<<" "<<rz<<std::endl;
      }
      else{
	outputalign<<", \""<<id<<"\": ["<<x<<", "<<y<<", "<<z<<", "<<rx<<", "<<ry<<", "<<rz<<"]";
	output<<id<<" "<<x<<" "<<y<<" "<<z<<" "<<rx<<" "<<ry<<" "<<rz<<std::endl;
      }
    }
    iline++;
    if(id==3){
      iflay=true;
    }
    if(id==32)ifmod=true;
  }
  outputalign<<std::endl;
  input.close();
  output.close();
  outputalign.close();
  std::cout<<"in total "<<iline<<" corrections"<<std::endl;
}
