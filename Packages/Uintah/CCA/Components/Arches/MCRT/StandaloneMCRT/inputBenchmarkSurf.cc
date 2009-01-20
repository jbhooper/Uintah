  
   //   *********  the benchmark coupled with surfaces case ***************

int fakeIndex = 0;
double xx, yy, zz;
for ( int k = 0; k < Ncz; k ++ ){
  for ( int j = 0; j < Ncy; j ++) {
    for ( int i = 0; i < Ncx; i ++ ) {
      
	T_Vol[fakeIndex] = 64.80721904; // k
	xx = (X[i] + X[i+1])/2;
	yy = (Y[j] + Y[j+1])/2;
	zz = (Z[k] + Z[k+1])/2;
	
	kl_Vol[fakeIndex] = 0.9 * ( 1 - 2 * abs ( xx ) )
	  * ( 1 - 2 * abs ( yy ) )
	  * ( 1 - 2 * abs ( zz ) ) + 0.1;
	
	a_Vol[fakeIndex] = 1;
	scatter_Vol[fakeIndex] = 0;
	fakeIndex++;
	
    }
  }
 }

// top bottom surface
for ( int i = 0; i < TopBottomNo; i ++ ) { 
  rs_top_surface[i] = 0.02;
  emiss_top_surface[i] = 0.9;
  alpha_top_surface[i] = emiss_top_surface[i]; // for gray diffuse surface??
  rd_top_surface[i] =  1 - rs_top_surface[i] - alpha_top_surface[i];
  T_top_surface[i] = 1200;
  a_top_surface[i] = 1;
  
  rs_bottom_surface[i] = 0.04;
  emiss_bottom_surface[i] = 0.8;
  alpha_bottom_surface[i] = emiss_bottom_surface[i]; // for gray diffuse surface??
  rd_bottom_surface[i] =  1 - rs_bottom_surface[i] - alpha_bottom_surface[i];
  T_bottom_surface[i] = 900;
  a_bottom_surface[i] = 1;
 }
   

// front back surface
for ( int i = 0; i < FrontBackNo; i ++ ) { 
  rs_front_surface[i] = 0.475;
  emiss_front_surface[i] = 0.05;
  alpha_front_surface[i] = emiss_front_surface[i]; // for gray diffuse surface??
  rd_front_surface[i] =  1 - rs_front_surface[i] - alpha_front_surface[i];
  T_front_surface[i] = 1400;
  a_front_surface[i] = 1;
  
  rs_back_surface[i] = 0.19;
  emiss_back_surface[i] = 0.05;
  alpha_back_surface[i] = emiss_back_surface[i]; // for gray diffuse surface??
  rd_back_surface[i] =  1 - rs_back_surface[i] - alpha_back_surface[i];
  T_back_surface[i] = 2000;
  a_back_surface[i] = 1;
 }



// left right surface
for ( int i = 0; i < LeftRightNo; i ++ ) { 
  rs_left_surface[i] = 0.76;
  emiss_left_surface[i] = 0.2;
  alpha_left_surface[i] = emiss_left_surface[i]; // for gray diffuse surface??
  rd_left_surface[i] =  1 - rs_left_surface[i] - alpha_left_surface[i];
  T_left_surface[i] = 600;
  a_left_surface[i] = 1;

  rs_right_surface[i] = 0.76;
  emiss_right_surface[i] = 0.2;
  alpha_right_surface[i] = emiss_right_surface[i]; // for gray diffuse surface??
  rd_right_surface[i] =  1 - rs_right_surface[i] - alpha_right_surface[i];
  T_right_surface[i] = 600;
  a_right_surface[i] = 1;
 }

// ************************* end of this case ********************
   
