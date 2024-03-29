#include "potential.h"
const Real potential::pi = 3.141592653589793238463;
potential::potential(int _N, int _M) :
  //  levmar(_N, _N),
  levmar(_N+1, _N+1),
  N(_N),
  M(_M),
  n1(10),
  theta1(N+1),
  theta2(M),
  //r0(N+1),
  //r1(N+1),
  //r2(N+1),
  //c0(N+1),
  c1(2*N,N+1),
  c2(2*M,M),
  V0(0.972*pi),
  alpha(6.0),
  // F(U), F=U(g=1) is our variable now.
  //  F(0.6),
  //U(-0.6),
  fft1(2*N),
  fft2(2*M),
  //  mu1(N+1),
  //mu2(M),
  //  G(0.0),
  x1(n1),
  x2(n1),
  w1(n1),
  w2(n1),
  //xx(n1),
  //ww(n1),
  xx1_1(n1,N),
  xx1_2(n1,N-1), // upper singularityy l==k
  xx1_3(n1,N-1), // lower singularity l-1==k
  ww1_1(n1,N),
  ww1_2(n1,N-1),
  ww1_3(n1,N-1),
  xx2_1(n1,M+1),
  xx2_2(n1,M),
  xx2_3(n1,M),
  ww2_1(n1,M+1),
  ww2_2(n1,M),
  ww2_3(n1,M),
  m0(10),
  index_B(12),
  index_D(12),
  B(21, 12, 0.0), // last two columns are for computing B*, D* neede for case m>0.9
  D(21, 12, 0.0)
  //A(N+1+M, N+1+M, 0.0),
  //rhs(N+1+M),
  //  _r(0.0, N+1)
{
  initial();
}

void potential::set_theta() {
  // theta1 length N+1
  theta1[0] = 0.0;
  Real temp = pi/N;
  for(int i=1; i<N; i++)
    theta1[i] = theta1[i-1]+temp;
  theta1[N] = pi;
  // theta2 length M
  temp = pi/M;
  theta2[0] = temp/2;
  for(int i=1; i<M; i++)
    theta2[i] = theta2[i-1]+temp;
}

void potential::set_c12() {
  // trignometric coeff of mu1
  mp_mat<Real> temp0(2*N, N+1, 0.0);
#pragma omp parallel for schedule(static)
  for(int i=0; i<=N; i++) {
    temp0(i,i) = 1.0;
    if (2*N-i < 2*N)
      temp0(2*N-i,i) = 1.0;
  }
  fft1.forward(temp0.p, c1.p, N+1, 2*N, 2*N);
  // trignometric coeff of mu2
  temp0.resize(2*M, M, 0.0);
#pragma omp parallel for schedule(static)
  for(int i=0; i<M; i++) 
    temp0(i,i) = temp0(2*M-i-1,i) = 1.0;
  fft2.forward(temp0.p, c2.p, M, 2*M, 2*M);
}


void potential::set_x0() {
  /*    
  int N1 = N/2;
  Real tmp;
  x[0] = -1.0;
  for(int i=1; i<=N1; i++)
    x[i] = 0.9;
  for(int i=N1+1; i<=N; i++) {
    tmp = atan(-20.0/3*tan(theta1[i]));
    x[i] = 36*pow(cos(tmp),2)+0.81*pow(sin(tmp),2);
    x[i] = sqrt(x[i]);
  }
  */
  /*
  int N1 = N/2;
  for(int i=0; i<=15; i++) 
    x[i] = 9.0/5*cos(i*pi/60);
  for(int i=16; i<=N1; i++)
    x[i] = 0.9;
  for(int i=N1+1; i<=N; i++)
    x[i] = x[N-i];
  */

 
  string line1;
  ifstream in("initial");
  for(int i=0; i<=N; i++) {
    getline(in, line1);
    istringstream yan(line1);
    yan >> x[i];
  }
    
    
  // rescale x, y, z of a ball by 3, 5/6, 5/6.
  /*    
  int N1 = N/2;
#pragma omp parallel for schedule(static)
  for(int i=0; i<=N1; i++) {
    Real tmp = atan(18.0/5*tan(theta1[i]));
    x[i] = 9*pow(cos(tmp),2)+25.0/36*pow(sin(tmp),2);
    x[i] = sqrt(x[i]);
  }  
#pragma omp parallel for schedule(static)
  for(int i=N1+1; i<=N; i++)
    x[i] = x[N-i];
  x[0] = -0.6;
  */
}

void potential::compute_c0r12(valarray<Real> &_r0,valarray<cmplx<Real> > &_c0, valarray<Real> &_r1, valarray<Real> &_r2) {
  valarray<Real> temp0(2*N), temp1(2*N);
  //#pragma omp parallel for
  for(int i=0; i<=N; i++)
    temp0[i] = _r0[i];
  //#pragma omp parallel for
  for(int i=N+1; i<2*N; i++) 
    temp0[i] = _r0[2*N-i];
  fft1.forward(&temp0[0], &_c0[0]); // store the fft coeff of R.
  /*
  FILE *fp1 = fopen("new", "w");
  valarray<cmplx<Real> > cc(N+1);
  valarray<Real> xxx(N+1, 0.0);
  for(int i=0; i<=N; i++) {
    cc[i] = _c0[i]*exp(-2*theta1[i]);
  }
  for(int i=0; i<=N; i++)
    xxx[i] = fft1.eval(&cc[0], theta1[i]);
  for(int i=0; i<=N; i++)
    fprintf(fp1, "%23s\n", str(xxx[i], 0));
  */
  /*
  FILE *fp1 = fopen("coeff", "w");
  for(int i=0; i<=N; i++) { // can change                                                                                                                   
    fprintf(fp1, "%23s\n",str(_c0[i].real(), 0));
  }
  */
  fft1.deriv(&temp0[0], &temp1[0]);

  //#pragma omp parallel for
  for(int i=0; i<=N; i++)
    _r1[i] = temp1[i];
  _r1[0] = _r1[N] = temp1[0] = temp1[N] = 0.0;
  
  fft1.deriv(&temp1[0], &temp0[0]);
  //#pagma omp parallel for
  for(int i=0; i<=N; i++)
    _r2[i] = temp0[i];
}

void potential::compute_V0() {
  // using x and trapezoidal rule to set volume
  Real temp0 = 0.0;
#pragma omp parallel for reduction (+: temp0)
  for(int i=1; i<N; i++)
    temp0 += pow(x[i],3)*sin(theta1[i]);

  V0 = temp0 * 2*pi/3 * pi/N;
}

void potential::compute_xw() {
  if (n1==10) {
    x1[0] = 0.013046735741414128;
    x1[1] = 0.067468316655507732;
    x1[2] = 0.16029521585048778;
    x1[3] = 0.28330230293537639;
    x1[4] = 0.42556283050918442;
    x1[5] = 0.57443716949081558;
    x1[6] = 0.71669769706462361;
    x1[7] = 0.83970478414951222;
    x1[8] = 0.93253168334449232;
    x1[9] = 0.98695326425858587;

    w1[0] = 0.033335672154342931;
    w1[1] = 0.074725674575290266;
    w1[2] = 0.10954318125799108;
    w1[3] = 0.13463335965499812;
    w1[4] = 0.14776211235737646;
    w1[5] = 0.14776211235737646;
    w1[6] = 0.13463335965499812;
    w1[7] = 0.10954318125799108;
    w1[8] = 0.074725674575290266;
    w1[9] = 0.033335672154342931;

	  
    x2[0] = 0.482961710689630e-3;
    x2[1] = 0.698862921431577e-2;
    x2[2] = 0.326113965946776e-1;
    x2[3] = 0.928257573891660e-1;
    x2[4] = 0.198327256895404;
    x2[5] = 0.348880142979353;
    x2[6] = 0.530440555787956;
    x2[7] = 0.716764648511655;
    x2[8] = 0.875234557506234;
    x2[9] = 0.975245698684393;

    w2[0] = 0.183340007378985e-2;
    w2[1] = 0.134531223459918e-1;
    w2[2] = 0.404971943169583e-1;
    w2[3] = 0.818223696589036e-1;
    w2[4] = 0.129192342770138;
    w2[5] = 0.169545319547259;
    w2[6] = 0.189100216532996;
    w2[7] = 0.177965753961471;
    w2[8] = 0.133724770615462;
    w2[9] = 0.628655101770325e-1;
#pragma omp parallel
    {
#pragma omp for collapse(2) schedule(static) nowait
      for(int j=0; j<N; j++) 
	for(int i=0; i<n1; i++) {
	  xx1_1(i, j) = pi/N*x1[i]+theta1[j];
	  ww1_1(i, j) = pi/N*w1[i];
	}
#pragma omp for collapse(2) schedule(static) nowait
      for(int j=1; j<N; j++) 
 	for(int i=0; i<n1; i++) {
	  xx1_2(i, j-1) = -pi/N*x2[i]+theta1[j];
	  ww1_2(i, j-1) = pi/N*w2[i];
	}
#pragma omp for collapse(2) schedule(static) nowait
      for(int j=1; j<N; j++) 
	for(int i=0; i<n1; i++) {
	  xx1_3(i, j-1) = pi/N*x2[i]+theta1[j];
	  ww1_3(i, j-1) = pi/N*w2[i];
	}
#pragma omp for schedule(static) nowait
      for(int i=0; i<n1; i++) {
	xx2_1(i, 0) = pi/2/M*x1[i];
	ww2_1(i, 0) = pi/2/M*w1[i];
	xx2_1(i, M) = -pi/2/M*x1[i]+pi;
	ww2_1(i, M) = pi/2/M*w1[i];
      }
#pragma omp for collapse(2) schedule(static) nowait
      for(int j=0; j<M-1; j++) 
	for(int i=0; i<n1; i++) {
	  xx2_1(i, j+1) = pi/M*x1[i]+theta2[j];
	  ww2_1(i, j+1) = pi/M*w1[i];
	}
#pragma omp for schedule(static) nowait
      for(int i=0; i<n1; i++) {
	xx2_2(i, 0) = -pi/M/2*x2[i]+theta2[0];
	ww2_2(i, 0) = pi/M/2*w2[i];
      }
#pragma omp for collapse(2) schedule(static) nowait
      for(int j=1; j<M; j++) 
	for(int i=0; i<n1; i++) {
	  xx2_2(i, j) = -pi/M*x2[i]+theta2[j];
	  ww2_2(i, j) = pi/M*w2[i];
	}
#pragma omp for schedule(static) nowait
      for(int i=0; i<n1; i++) {
	xx2_3(i,M-1) = pi/2/M*x2[i]+theta2[M-1];
	ww2_3(i,M-1) = pi/2/M*w2[i];
      }
#pragma omp for collapse(2) schedule(static) nowait
      for(int j=0; j<M-1; j++) 
	for(int i=0; i<n1; i++) {
	  xx2_3(i, j) = pi/M*x2[i]+theta2[j];
	  ww2_3(i, j) = pi/M*w2[i];
	}
#pragma omp barrier
    }
  }
    // if need higher degree add n==15, 20 case
  
}

void potential::set_BD() {
  m0[0] = 0.05;
  index_B[0] = 12;
  index_D[0] = 12;
  B(0,0) = 0.790401413584395132;
  D(0,0) = 0.800602040206397048;
  B(1,0) = 0.102006266220019155;
  D(1,0) = 0.313994477771767757;
  B(2,0) = 0.039878395558551461;
  D(2,0) = 0.205913118705551955;
  B(3,0) = 0.021737136375982167;
  D(3,0) = 0.157744346538923994;
  B(4,0) = 0.013960979767622058;
  D(4,0) = 0.130595077319933092;
  B(5,0) = 0.009892518822669142;
  D(5,0) = 0.113308474489758567;
  B(6,0) = 0.007484612400663336;
  D(6,0) = 0.101454199173630195;
  B(7,0) = 0.005934625664295474;
  D(7,0) = 0.092918784207297437;
  B(8,0) = 0.004874249053581664;
  D(8,0) = 0.086565380148168087;
  B(9,0) = 0.004114606930310886;
  D(9,0) = 0.081727984665103014;
  B(10,0) = 0.003550452989196177;
  D(10,0) = 0.077990665729107038;
  B(11,0) = 0.003119229959988475;
  D(11,0) = 0.075080426851268007;

  m0[1] = 0.15;
  index_B[1] = 12;
  index_D[1] = 12;
  B(0,1) = 0.801024064452844894;
  D(0,1) = 0.834232667811735098;
  B(1,1) = 0.110695344529634015;
  D(1,1) = 0.360495281619098276;
  B(2,1) = 0.047348746716993718;
  D(2,1) = 0.262379664114505869;
  B(3,1) = 0.028484367255041423;
  D(3,1) = 0.223723944518094276;
  B(4,1) = 0.020277811444003597;
  D(4,1) = 0.206447811775681053;
  B(5,1) = 0.015965005853099119;
  D(5,1) = 0.199809440876486856;
  B(6,1) = 0.013441320273553635;
  D(6,1) = 0.199667451603795275;
  B(7,1) = 0.011871565736951440;
  D(7,1) = 0.204157558868236842;
  B(8,1) = 0.010868363672485521;
  D(8,1) = 0.212387467960572375;
  B(9,1) = 0.010231587232710565;
  D(9,1) = 0.223948914061499360;
  B(10,1) = 0.009849585546666211;
  D(10,1) = 0.238708097425597860;
  B(11,1) = 0.009656606347153765;
  D(11,1) = 0.256707203545463756;
  
  m0[2] = 0.25;
  index_B[2] = 13;
  index_D[2] = 13;
  B(0,2) = 0.812597772919920493;
  D(0,2) = 0.873152581892675550;
  B(1,2) = 0.121109617945510113;
  D(1,2) = 0.420622230667770216;
  B(2,2) = 0.057293376831239877;
  D(2,2) = 0.344231061559450379;
  B(3,2) = 0.038509451602167328;
  D(3,2) = 0.331133021818721762;
  B(4,2) = 0.030783430301775233;
  D(4,2) = 0.345277285052808412;
  B(5,2) = 0.027290564934732527;
  D(5,2) = 0.377945322150393392;
  B(6,2) = 0.025916369289445199;
  D(6,2) = 0.427378012464553881;
  B(7,2) = 0.025847203343361799;
  D(7,2) = 0.494671744307822406;
  B(8,2) = 0.026740923539348855;
  D(8,2) = 0.582685115665646201;
  B(9,2) = 0.028464314554825705;
  D(9,2) = 0.695799207728083165;
  B(10,2) = 0.030995446237278954;
  D(10,2) = 0.840018401472533403;
  B(11,2) = 0.034384369179940976;
  D(11,2) = 1.023268503573606061;
  B(12,2) = 0.038738002072493936;
  D(12,2) = 1.255859085136282496;

  m0[3] = 0.35;
  index_B[3] = 13;
  index_D[3] = 14;
  B(0,3) = 0.825323557983515895;
  B(1,3) = 0.133862116083687790;
  B(2,3) = 0.071011293597988675;
  B(3,3) = 0.054178477417387376;
  B(4,3) = 0.049451744948102993;
  B(5,3) = 0.050222196224107476;
  B(6,3) =  0.054742913171830353;
  B(7,3) = 0.062746257927001699;
  B(8,3) = 0.074669881043476886;
  B(9,3) = 0.091480845177733472;
  B(10,3) = 0.114705092110997824;
  B(11,3) = 0.146571132581439876;
  B(12,3) = 0.190257137333846284;
  D(0,3) = 0.919027039242097348;
  D(1,3) = 0.501002159288247514;
  D(2,3) = 0.468831270566456863;
  D(3,3) = 0.517714227776400015;
  D(4,3) = 0.620843391317303107;
  D(5,3) = 0.782364393786869723;
  D(6,3) = 1.019114535076102913;
  D(7,3) = 1.359345202748496052;
  D(8,3) = 1.845717302358827942;
  D(9,3) = 2.541071703153920729;
  D(10,3) = 3.537404655208041337;
  D(11,3) = 4.969296002977425930;
  D(12,3) = 7.033822870030031126;
  D(13,3) = 10.02004322503447140;

  m0[4] = 0.45;
  index_B[4] = 13;
  index_D[4] = 16;
  B(0,4) = 0.839479570270612971;
  B(1,4) = 0.149916440306396336;
  B(2,4) = 0.090831935819428835;
  B(3,4) = 0.080347033483341786;
  B(4,4) = 0.085638440500470454;
  B(5,4) = 0.101954725932990372;
  B(6,4) = 0.130574811533616015;
  B(7,4) = 0.176105076358849928;
  B(8,4) = 0.246835164402955447;
  B(9,4) = 0.356424476867718855;
  B(10,4) = 0.527002562230102743;
  B(11,4) = 0.794389634259304750;
  B(12,4) = 1.216762532429718021;
  
  D(0,4) = 0.974404366546369673;
  D(1,4) = 0.613246805394160910;
  D(2,4) = 0.671096669502166996;
  D(3,4) = 0.870727620185086140;
  D(4,4) = 1.229542231202690761;
  D(5,4) = 1.826605967544420569;
  D(6,4) = 2.806934530997762740;
  D(7,4) = 4.418789329084028134;
  D(8,4) = 7.083236057478765325;
  D(9,4) = 11.51508812055758294;
  D(10,4) = 18.93151118599927464;
  D(11,4) = 31.41199693820496388;
  D(12,4) = 52.52072945457582854;
  D(13,4) = 88.38485473506529806;
  D(14,4) = 149.5663744939804784;
  D(15,4) = 254.3179084310411743;

  m0[5] = 0.55;
  index_B[5] = 14;
  index_D[5] = 17;
  B(0,5) = 0.855469615156419991;
  B(1,5) = 0.170896072689739584;
  B(2,5) = 0.121335229026948226;
  B(3,5) = 0.128201883574947410;
  B(4,5) = 0.164687281451527560;
  B(5,5) = 0.237418908749381742;
  B(6,5) = 0.369208104716495452;
  B(7,5) = 0.605658733847927717;
  B(8,5) = 1.033705561557812744;
  B(9,5) = 1.818988489363267885;
  B(10,5) = 3.279377651273850938;
  B(11,5) = 6.029888380717536331;
  B(12,5) = 11.26979685557794172;
  B(13,5) = 21.35457785038283450;
  D(0,5) = 1.043455295115133534;
  D(1,5) = 0.779625721928504850;
  D(2,5) = 1.029742360932067582;
  D(3,5) = 1.622037223411353130;
  D(4,5) = 2.787989531185347620;
  D(5,5) = 5.048381487372069147;
  D(6,5) = 9.463277611943484295;
  D(7,5) = 18.18148994942766790;
  D(8,5) = 35.58098059117916870;
  D(9,5) = 70.63393546191445013;
  D(10,5) = 141.8285800834330593;
  D(11,5) = 287.4487512501321663;
  D(12,5) = 587.1153846499230762;
  D(13,5) = 1207.065435225480616;
  D(14,5) = 2495.588727248664223;
  D(15,5) = 5184.692429394806441;
  D(16,5) = 10817.21333690413275;

  m0[6] = 0.65;
  index_B[6] = 16;
  index_D[6] = 18;
  B(0,6) = 0.873920061848643136;
  B(1,6) = 0.199814057482376946;
  B(2,6) = 0.172769615878015213;
  B(3,6) = 0.228106913284202167;
  B(4,6) = 0.370468141118071220;
  B(5,6) = 0.679271252884820555;
  B(6,6) = 1.348008496681757302;
  B(7,6) = 2.827670976853820704;
  B(8,6) = 6.179468250123914084;
  B(9,6) = 13.93568601034281150;
  B(10,6) = 32.21892928105972203;
  B(11,6) = 76.00696295922610103;
  B(12,6) = 182.3214490877540696;
  B(13,6) = 443.5150764411264816;
  B(14,6) = 1091.854722902838829;
  B(15,6) = 2715.765866403819588;

  D(0,6) = 1.133678336575733166;
  D(1,6) = 1.048643173729970391;
  D(2,6) = 1.753465041198464516;
  D(3,6) = 3.523182726803385513;
  D(4,6) = 7.749476413813974582;
  D(5,6) = 17.98645005585073306;
  D(6,6) = 43.25591634623261333;
  D(7,6) = 106.6815344540960170;
  D(8,6) = 268.0984865731174340;
  D(9,6) = 683.6241148502898048;
  D(10,6) = 1763.497085219187407;
  D(11,6) = 4592.374753831163809;
  D(12,6) = 12053.44101904888928;
  D(13,6) = 31846.66302074208170;
  D(14,6) = 84621.22135905680802;
  D(15,6) = 225956.4231829078900;
  D(16,6) = 605941.5172817588600;
  D(17,6) = 1631082.599539268321;

  m0[7] = 0.75;
  index_B[7] = 19;
  index_D[7] = 21;
  B(0,7) = 0.895902820924731621;
  B(1,7) = 0.243140003766786662;
  B(2,7) = 0.273081875594105532;
  B(3,7) = 0.486280007533573324;
  B(4,7) = 1.082747437228230918;
  B(5,7) = 2.743445290986452500;
  B(6,7) = 7.555817828670234627;
  B(7,7) = 22.05194082493752427;
  B(8,7) = 67.15640644740229408;
  B(9,7) = 211.2722537881770962;
  B(10,7) = 681.9037843053270682;
  B(11,7) = 2246.956231592536517;
  B(12,7) = 7531.483865999711792;
  B(13,7) =  25608.51260130241579;
  B(14,7) = 88140.74740089604971;
  B(15,7) = 306564.4242098446591;
  B(16,7) = 1076036.077811072194;
  B(17,7) = 3807218.502573632648;
  B(18,7) = 13566382.24422139551;
  D(0,7) = 1.260612826574911614;
  D(1,7) = 1.548665638082676581;
  D(2,7) = 3.553669411871607615;
  D(3,7) = 9.900444676104398756;
  D(4,7) = 30.32056661745247199;
  D(5,7) = 98.18025865888308915;
  D(6,7) = 329.7710104345570550;
  D(7,7) = 1136.655989742890393;
  D(8,7) = 3993.834335746229798;
  D(9,7) = 14242.72958655527085;
  D(10,7) = 51394.75729168872096;
  D(11,7) = 187246.7029146231521;
  D(12,7) = 687653.0923753899027;
  D(13,7) = 2542385.535653982270;
  D(14,7) = 9453781.219347490272;
  D(15,7) = 35328363.01797091708;
  D(16,7) = 132593262.3833930149;
  D(17,7) = 499544968.1840548215;
  D(18,7) = 1888409347.294438724;
  D(19,7) = 7160267534.478937192;
  D(20,7) = 27223307946.96339622;

  m0[8] = 0.825;
  index_B[8] = 15;
  index_D[8] = 18;
  B(0,8) = 0.915922052601931494;
  B(1,8) = 0.294714252429483394;
  B(2,8) = 0.435776709264636140;
  B(3,8) = 1.067328246493644239;
  B(4,8) = 3.327844118563268085;
  B(5,8) = 11.90406004445092906;
  B(6,8) = 46.47838820224626394;
  B(7,8) = 192.7556002578809477;
  B(8,8) = 835.3356299261900064;
  B(9,8) = 3743.124548343029103;
  B(10,8) = 17219.07731004063094;
  B(11,8) = 80904.60401669850158;
  B(12,8) = 386808.3292751742460;
  B(13,8) = 1876487.670110449342;
  B(14,8) = 9216559.908641567755;
  D(0,8) = 1.402200569110579095;
  D(1,8) = 2.322205897861749447;
  D(2,8) = 7.462158366466719683;
  D(3,8) = 29.43506890797307903;
  D(4,8) = 128.1590924337895775;
  D(5,8) = 591.0807036911982326;
  D(6,8) = 2830.546229607726377;
  D(7,8) = 13917.76431889392230;
  D(8,8) = 69786.10525163921228;
  D(9,8) = 355234.1420341879635;
  D(10,8) = 1830019.186413931054;
  D(11,8) = 9519610.812032515607;
  D(12,8) = 49920868.75574849454;
  D(13,8) = 263567700.9826023474;
  D(14,8) = 1399645765.120061119;
  D(15,8) = 7469935792.837635005;
  D(16,8) = 40041555958.35610574;
  D(17,8) = 215463066814.4966654;

  m0[9] = 0.875;
  index_B[9] = 19;
  index_D[9] = 21;
  B(0,9) = 0.931906061029524828;
  B(1,9) = 0.348448029538453861;
  B(2,9) = 0.666809178846938248;
  B(3,9) = 2.210769135708128663;
  B(4,9) = 9.491765048913406881;
  B(5,9) = 47.09304791027740853;
  B(6,9) = 255.9200460211233087;
  B(7,9) = 1480.029532675805408;
  B(8,9) = 8954.040904734313578;
  B(9,9) = 56052.48220982686950;
  B(10,9) = 360395.7241626000917;
  B(11,9) = 2367539.415273216078;
  B(12,9) = 15829949.57277684102;
  B(13,9) = 107415809.3278511100;
  B(14,9) = 738058546.0239595692;
  B(15,9) = 5126022002.555101497;
  B(16,9) = 35935340655.02416589;
  B(17,9) = 253988125761.2812212;
  B(18,9) = 1808180007145.359570;
  D(0,9) = 1.541690112721819084;
  D(1,9) = 3.379176214579645449;
  D(2,9) = 14.94058385670236672;
  D(3,9) = 81.91773929235074881;
  D(4,9) = 497.4900546551479866;
  D(5,9) = 3205.184010234846235;
  D(6,9) = 21457.32237355321926;
  D(7,9) = 147557.0156564174712;
  D(8,9) = 1035045.290185256525;
  D(9,9) = 7371922.334832212125;
  D(10,9) = 53143443.95142401142;
  D(11,9) = 386882347.5795976313;
  D(12,9) = 2839458401.528033778;
  D(13,9) = 20982661229.43898942;
  D(14,9) = 155961775401.7662418;
  D(15,9) = 1165096220419.884791;
  D(16,9) = 8742012983013.913805;
  D(17,9) = 65847254626723.66919;
  D(18,9) = 497679873706243.4393;
  D(19,9) = 3773018634056605.405;
  D(20,9) = 28682631948378196.60;

  index_B[10] = index_B[11] = 14;
  index_D[10] = index_D[11] = 13;
  B(0,10) = 0.0;
  B(1,10) = -1.0/4;
  B(2,10) = -1.0/32;
  B(3,10) = -3.0/256;
  B(4,10) = -25.0/4096;
  B(5,10) = -245.0/65536;
  B(6,10) = -1323.0/524288;
  B(7,10) = -7623.0/4194304;
  B(8,10) = -184041.0/134217728;
  B(9,10) = -4601025.0/4294967296;
  B(10,10) = -29548805.0/34359738368;
  B(11,10) = -193947611.0/274877906944;
  B(12,10) = -2591845347.0/4398046511104;
  B(13,10) = -35156056117.0/70368744177664;
  
  B(0,11) = 1.0;
  B(1,11) = -1.0/4;
  B(2,11) = 3.0/64;
  B(3,11) = 3.0/128;
  B(4,11) = 665.0/49152;
  B(5,11) = 3437.0/393216;
  B(6,11) = 15981.0/2621440;
  B(7,11) = 188287.0/41943040;
  B(8,11) = 129334777.0/37580963840;
  B(9,11) = 327273375.0/120259084288;
  B(10,11) = 19096474969.0/8658654068736;
  B(11,11) = 631505527133.0/346346162749440;
  B(12,11) = 2224154230753.0/1451355348664320;
  B(13,11) = 181962561086453.0/139330113471774720;

  D(0,10) = 1.0/2;
  D(1,10) = -1.0/8;
  D(2,10) = -3.0/128;
  D(3,10) = -5.0/512;
  D(4,10) = -175.0/32768;
  D(5,10) = -441.0/131072;
  D(6,10) = -4851.0/2097152;
  D(7,10) = -14157.0/8388608;
  D(8,10) = -2760615.0/2147483648;
  D(9,10) = -8690825.0/8589934592;
  D(10,10) = -112285459.0/137438953472;
  D(11,10) = -370263621.0/549755813888;
  D(12,10) = -19870814327.0/35184372088832;
  
  D(0,11) =  -1.0;
  D(1,11) = 0.0;
  D(2,11) = 5.0/128;
  D(3,11) = 31.0/1536;
  D(4,11) = 2365.0/196608;
  D(5,11) = 10409.0/1310720;
  D(6,11) = 117929.0/20971520;
  D(7,11) = 2458621.0/587202560;
  D(8,11) = 194646309.0/60129542144;
  D(9,11) = 5577961675.0/2164663517184;
  D(10,11) = 363577654297.0/173173081374720;
  D(11,11) = 632563423193.0/362838837166080;
  D(12,11) = 102453646108723.0/69665056735887360;
}

void potential::initial() {
  set_theta();
  set_c12();
  set_x0();
  //compute_c0r12(r0, r1, r2);
  //compute_V0();
  compute_xw();
  set_BD();
}

Real potential::compute_B(Real _m, Real _multi) {
  if(_m == 0.0)
    return pi/4;
  if(m == 1.0)
    return 1.0;
  Real temp = 0.0, multi, X, temp1;
  for(int i=0; i<8; i++) {
    if(_m > i*0.1 && _m <= (i+1)*0.1) {
      multi = _m - m0[i];
      for(int j=index_B[i]; j>0; j--) {
	temp *= multi;
	temp += B(j-1,i);
      }
      return temp;
    }
  }
  if(_m<=0.85) {
    multi = _m - m0[8];
    for(int j=index_B[8]; j>0; j--) {
      temp *= multi;
      temp += B(j-1,8);
    }
    return temp;
  }
  else if (_m <= 0.9) {
    multi = _m - m0[9];
    for(int j=index_B[9]; j>0; j--) {
      temp *= multi;
      temp += B(j-1,9);
    }
    return temp;
  }
  else {
    multi = _multi;
    for(int j=index_B[10]; j>0; j--) {
      temp *= multi;
      temp += B(j-1,10);
    }
    temp1 = 0.0;
    for(int j=index_B[11]; j>0; j--) {
      temp1 *= multi;
      temp1 += B(j-1,11);
    }
    X = -log(multi/16);
    temp1 += temp*X;
    temp1 /= _m;
    return temp1;
  } 
}


Real potential::compute_D(Real _m, Real _multi) {
  if(_m == 0.0)
    return pi/4;
  Real temp=0, multi = 0.0, X, temp1;
  for(int i=0; i<8; i++) {
    if(_m > i*0.1 && _m <= (i+1)*0.1) {
      multi = _m - m0[i];
      for(int j=index_D[i]; j>0; j--) {
	temp *= multi;
	temp += D(j-1,i);
      }
      return temp;
    }
  }
  if(_m<=0.85) {
    multi = _m - m0[8];
    for(int j=index_D[8]; j>0; j--) {
      temp *= multi;
      temp += D(j-1,8);
    }
    return temp;
  }
  else if (_m <= 0.9) {
    multi = _m - m0[9];
    for(int j=index_D[9]; j>0; j--) {
      temp *= multi;
      temp += D(j-1,9);
    }
    return temp;
  }
  else {
    multi = _multi;
    for(int j=index_D[10]; j>0; j--) {
      temp *= multi;
      temp += D(j-1,10);
    }
    temp1 = 0.0;
    for(int j=index_D[11]; j>0; j--) {
      temp1 *= multi;
      temp1 += D(j-1,11);
    }
    X = -log(multi/16);
    temp1 += temp*X;
    temp1 /= _m;
    return temp1;
  }
}

Real potential::compute_G(Real _a, Real _b, Real _c, Real _d) {
  // checked with mathematica
  Real G;
  if (_d == 0.0) {
    G = 2*pi*_a/pow(_c, 1.5);
    return G;
  }
  Real b0, d0;
  Real m1 = -2.0*_d/(_c-_d);
  if (m1<0) {
    Real tmp= -m1/(1-m1);
    if (tmp <= 0.9) {
      b0 = compute_B(tmp, 1.0/(1-m1));
      d0 = compute_D(tmp, 1.0/(1-m1));
    }
    else {
      tmp = 2.0*_d/(_c+_d);
      b0 = compute_B(tmp, (_c-_d)/(_c+_d));
      d0 = compute_D(tmp, (_c-_d)/(_c+_d));
    }
      
    tmp = sqrt(1-m1);
    G = 4.0*_b/(_d*sqrt(_c-_d))*(b0+d0)/tmp;
    G -= 4.0*(_b*_c-_a*_d)/((_c+_d)*_d*sqrt(_c-_d))*(b0+d0/(1-m1))*tmp;
  }
  else {
    b0 = compute_B(m1, 1-m1);
    d0 = compute_D(m1, 1-m1);
    G = 4.0*_b/(_d*sqrt(_c-_d))*(b0+d0);
    G -= 4.0*(_b*_c-_a*_d)/((_c+_d)*_d*sqrt(_c-_d))*(b0+(1-m1)*d0);
  }
  return G; 
}

valarray<Real> potential::compute_rr(valarray<Real> &r0, Real U) {
//valarray<Real> potential::compute_rr(valarray<Real> &r0) {
  // compute phi1_theta and phi2_theta
  valarray<cmplx<Real> > c0(N+1);
  valarray<Real> r1(N+1), r2(N+1);
  compute_c0r12(r0,c0,r1,r2);
  
  valarray<Real> mu1(N+1), mu2(M);

  //compute_mu(_r0, mu1, mu2);
  mp_mat<Real> A(N+1+M, N+1+M, 0.0);
  valarray<Real> rhs(N+1+M);

  Real val, y, temp, phi1, phi2;
  Real a, b, c, d, e, f, g, h, G;
  valarray<Real> xx(0.0, n1), ww(0.0, n1);
  valarray<Real> _r(0.0, N+1);
  //valarray<Real> _r(0.0, N); //without volumn
  Real F = -U/0.00025;
  
  rhs[0] = U*r0[0];
  rhs[N] = -U*r0[N];
  A(0,0) = 0.5*r0[0];
  A(N,N) = 0.5*r0[N];
  A(0,0) += pi*(r0[0]-r2[0])*pi/(2*N);
  A(N,N) += pi*(r0[N]-r2[N])*pi/(2*N);

  for(int j=1; j<N; j++) {
    a = -r0[0]*cos(theta1[j])*r0[j]+ pow(r0[0], 2);
    c = r0[j]*r0[j] + pow(r0[0],2) -2*r0[j]*r0[0]*cos(theta1[j]);
    G = compute_G(a, 0, c, 0);
    A(0, j) += pi/N*pow(r0[j],2)*sin(theta1[j])*G;
  }
  for(int j=1; j<N; j++) {
    a = r0[N]*cos(theta1[j])*r0[j] + r0[N]*r0[N];
    c = r0[j]*r0[j]+r0[N]*r0[N]+2*r0[N]*r0[j]*cos(theta1[j]);
    G = compute_G(a,0,c,0);
    A(N, j) += pi/N*r0[j]*r0[j]*sin(theta1[j])*G;
  }
  
  for(int q=0; q<M; q++) {
    y = -tan(theta2[q]-pi/2);
    e = r0[0]*(r0[0]-y);
    g = y*y+r0[0]*r0[0]-2*y*r0[0]+1;
    G = compute_G(e,0,g,0);
    A(0, q+N+1) += (1+y*y)*G*pi/M;

    e = r0[N]*(r0[N]+y);
    g = y*y+r0[N]*r0[N]+2*y*r0[N]+1;
    G = compute_G(e,0,g,0);
    A(N, q+N+1) += (1+y*y)*G*pi/M;
  }
  
  for(int k=1; k<N; k++) {
    temp = sqrt(r0[k]*r0[k]+r1[k]*r1[k]);
    Real common1 = r1[k]*sin(theta1[k])+r0[k]*cos(theta1[k]);
    Real common2 = r0[k]*sin(theta1[k])-r1[k]*cos(theta1[k]);
    rhs[k] = common1*U;
    A(k,k) = 0.5*temp;
    for(int l=1; l<=N; l++) {
      if(k==l) {
	xx = ww1_2.extract_column(k-1);
	ww = ww1_2.extract_column(k-1);
      }
      else if(l== k+1) {
	xx = xx1_3.extract_column(k-1);
	ww = ww1_3.extract_column(k-1);
      }
      else {
	xx = xx1_1.extract_column(l-1);
	ww = ww1_1.extract_column(l-1);
      }
      for(int q=0; q<n1; q++) {
	val = fft1.eval(&c0[0], xx[q]);
	a = -val*cos(xx[q])*common1 + r0[k]*r0[k];
	b = val*sin(xx[q])*common2;
	c = val*val + pow(r0[k],2) -2*val*r0[k]*cos(xx[q])*cos(theta1[k]);
	d = 2*val*r0[k]*sin(xx[q])*sin(theta1[k]);
	G = compute_G(a, b, c, d);
	for(int j=0; j<N+1; j++)
	  A(k, j) += val*val*fft1.eval(&c1.p[2*N*j], xx[q])*sin(xx[q])*G*ww[q];
      }
    }
    // for this second integral, use midpoint rule due to the smooth integrand
     for(int j=0; j<M; j++) {
      y = -tan(theta2[j]-pi/2);
      e = -common1*y+r0[k]*r0[k];
      f = common2;
      g = y*y+pow(r0[k],2)-2*y*r0[k]*cos(theta1[k])+1;
      h = 2*r0[k]*sin(theta1[k]);
      G = compute_G(e, f, g, h);

      val = (1+y*y)*G*pi/M;
      A(k, j+N+1) += val; 
    }
  }
  
  //second set of eqns
  Real x_k;
  for(int k=0; k<M; k++){
    x_k = -tan(theta2[k]-pi/2);
    rhs[N+1+k] = 0;
    A(N+1+k, N+1+k) = 0.5;
    for(int j=1; j<N; j++) {
      // trapezoidal rule
      val = r0[j]*r0[j]*sin(theta1[j]);
      a = -1;
      b = -r0[j]*sin(theta1[j]);
      c = pow(x_k,2)+r0[j]*r0[j]-2*x_k*r0[j]*cos(theta1[j])+1;
      d = 2*r0[j]*sin(theta1[j]);
      G = compute_G(a, b, c, d);
      A(N+1+k, j) += val*G*pi/N;
    }
    for(int l=0; l<=M; l++) {      
      if (l==k) {
	xx = xx2_2.extract_column(k);
	ww = ww2_2.extract_column(k);
      }
      else if (l-1 == k) {
	xx = xx2_3.extract_column(k);
	ww = ww2_3.extract_column(k);
      }
      else {
	xx = xx2_1.extract_column(l);
	ww = ww2_1.extract_column(l);
      }
      for(int q=0; q<n1; q++) {
	y = -tan(xx[q]-pi/2);
	e = -1.0;
	f = -1.0;
	g = pow(y-x_k,2)+2;
	h = 2.0;
	G = compute_G(e, f, g, h);
	for(int j=0; j<M; j++) {
	  val = fft2.eval(&c2.p[2*M*j], xx[q]-pi/(2*M));
	  val *= G*ww[q]*(1+y*y);
	  A(N+1+k, N+1+j) += val;
	}
      }
    }
  }

  //A.dump("AA", 17, '%');
  // solve A\mu = rhs
  int info = 0;
  valarray<int> ipv(N+1+M);
  dgetrf(N+1+M,N+1+M,A.p,N+1+M,&ipv[0],&info);
  dgetrs('N',N+1+M,1,A.p,N+1+M,&ipv[0],&rhs[0],N+1+M,&info);
  //for(int i=0; i<N+1+M; i++)
  //printf("%23s\n", str(rhs[i],0) );
  for(int i=0; i<=N; i++)
    mu1[i] = rhs[i];
  for(int i=0; i<M; i++)
    mu2[i] = rhs[N+1+i];
  
 
  //#pragma omp for 
  for(int k=1; k<N; k++) {
    phi1 = 0.0;
    phi2 = 0.0;
    
    for(int l=1; l<=N; l++) {
      /*
      if (k==l) {
	  xx = -(pi/N)*x2+theta1[k];
	  ww = (pi/N)*w2;
	}
	else if (k==l-1) {
	  xx = (pi/N)*x2+theta1[k];
	  ww = (pi/N)*w2;
	}
	else {
	  xx = pi/N*x1+theta1[l-1];
	  ww = pi/N*w1;
	}
      */
      
      if(k==l) {
	xx = xx1_2.extract_column(k-1);
	ww = ww1_2.extract_column(k-1);
      }
      else if (k==l-1) {
	xx = xx1_3.extract_column(k-1);
	ww = ww1_3.extract_column(k-1);
      }
      else {
	xx = xx1_1.extract_column(l-1);
	ww = ww1_1.extract_column(l-1);
      }
      
      for(int q=0; q<n1; q++) {
	val = fft1.eval(&c0[0], xx[q]);
	a = r0[k]*r1[k]-val*r1[k]*cos(xx[q])*cos(theta1[k])+val*r0[k]*cos(xx[q])*sin(theta1[k]);
	b = val*r1[k]*sin(xx[q])*sin(theta1[k])+val*r0[k]*sin(xx[q])*cos(theta1[k]);
	c = val*val + pow(r0[k], 2) - 2*r0[k]*val*cos(theta1[k])*cos(xx[q]);
	d = 2*r0[k]*val*sin(theta1[k])*sin(xx[q]);
	G = compute_G(a, b, c, d);	  
	temp = 0.0;
	for(int j=0; j<=N; j++)
	  temp += mu1[j]*fft1.eval(&c1.p[2*N*j],xx[q]);
	temp *= val*val*sin(xx[q])*G;
	temp += 2*mu1[k]*r0[k]/(xx[q]-theta1[k]);
	phi1 -= temp*ww[q];
      }
    }
    phi1 += 2*mu1[k]*r0[k]*(log(pi-theta1[k])-log(theta1[k]));
    

      // trapezoidal rule
    for(int j=0; j<M; j++) {
      y = -tan(theta2[j]-pi/2);
      e = r0[k]*r1[k]-y*(r1[k]*cos(theta1[k])-r0[k]*sin(theta1[k]));
      f = r1[k]*sin(theta1[k])+r0[k]*cos(theta1[k]);
      g = y*y+pow(r0[k],2)-2*y*r0[k]*cos(theta1[k])+1;
      h = 2*r0[k]*sin(theta1[k]);
      G = compute_G(e, f, g, h);
      
      phi2 -= mu2[j]*(1+y*y)*G*pi/M;
    }
  
    temp = r0[k]*r0[k]+r1[k]*r1[k];
     _r[k-1] = pow(phi1+phi2+U*r1[k]*cos(theta1[k])-U*r0[k]*sin(theta1[k]),2)/temp;
     _r[k-1] += 2.0/F/F*r0[k]*cos(theta1[k]);
     _r[k-1] += 2.0/alpha*((r0[k]*r0[k]+2*r1[k]*r1[k]-r0[k]*r2[k])/pow(temp,1.5) + (1-r1[k]/r0[k]/tan(theta1[k]))/sqrt(temp));
     _r[k-1] -= 2.0/F/F*r0[0]+4.0/alpha*(r0[0]-r2[0])/pow(r0[0],2);
     //_r[k-1] -= 4.0/alpha*(r0[0]-r2[0])/pow(r0[0],2);
  }
  
  //_r[N-1] = 4.0/alpha*((r0[N]-r2[N])/pow(r0[N],2)-(r0[0]-r2[0])/pow(r0[0],2));
  _r[N-1] = -2.0/F/F*(r0[N]+r0[0])+4.0/alpha*((r0[N]-r2[N])/pow(r0[N],2)-(r0[0]-r2[0])/pow(r0[0],2));
  
  
  temp = 0.0;
  
  //#pragma omp parallel for reduction (+:temp)
  for(int i=1; i<N; i++)
    temp += pow(r0[i],3) * sin(theta1[i]);

   //bubble volume fixed.
   _r[N] = temp*2*pi/3*pi/N - V0;
  
  /*
  for(int i=0; i<=N; i++)
    printf(" %23s\n", str(mu1[i],0));
  for(int i=0; i<M; i++)
    printf("%23s\n", str(mu2[i],0));
  */
  return _r;
}

void potential::compute_r() {
  //update_shape();
  
  valarray<Real> r0(N+1);
  for(int i=1; i<=N; i++)
    r0[i] = x[i];
  r0[0] = 0.9;
  valarray<Real> temp_r(N+1);
  temp_r = compute_rr(r0, x[0]);
  for(int i=0; i<=N; i++)
    r[i] = temp_r[i];  
  /*
  valarray<Real> r0(N+1);
  for(int i=1; i<=N; i++)                                                                                                            
    r0[i] = x[i-1];                                                                                                                  
  r0[0] = 3.0; 
  valarray<Real> temp_r(N);
  temp_r = compute_rr(r0);
  for(int i=0; i<N; i++)
    r[i] = temp_r[i];
  */
}
/*
void potential::update_r(vector<Real> &_r, valarray<Real> _rr) {
  //#pragma omp parallel for
  for(int i=0; i<=N; i++) 
    _r[i] = _rr[i];
}
*/
void potential::compute_J() {
  //update_shape();
  //valarray<Real> g1(N+1), g2(N+1);
  const Real eps = 1.0e-8;//6.0e-6;
  
#pragma omp parallel for schedule(static)
  //for(int i=0; i<N; i++) {
  for(int i=1; i<=N; i++) {
    valarray<Real> g1(N+1), temp_r0(N+1);
    temp_r0[0] = 0.9;
    for(int k=1; k<=N; k++)
      temp_r0[k] = x[k];
    
    temp_r0[i] += eps;
    g1 = compute_rr(temp_r0, x[0]);
    //for(int l=0; l<=N; l++)
    //printf("%d, %23s\nn", l, str(g[l],0));
    
    temp_r0[i] = temp_r0[i] - 2*eps;
    g1 -= compute_rr(temp_r0, x[0]);
    for(int k=0; k<=N; k++) 
      J(k,i) = g1[k]/(2*eps);
  }
  
  
  valarray<Real> g1(N+1), temp_r0(N+1);
  for(int k=1; k<=N; k++)
      temp_r0[k] = x[k];
  temp_r0[0] = 0.9;
  g1 = compute_rr(temp_r0, x[0]+eps);
  g1 -= compute_rr(temp_r0, x[0]-eps);
  for(int k=0; k<=N; k++) {
    J(k,0) = g1[k]/(2*eps);
  }
  //J.dump("JJ_parallel", 17, '%');
}
	  

void potential::interpolate(int l) {
  valarray<cmplx<Real> > c0(N+1);
  valarray<Real> r1(N+1), r2(0.0,N+1);
  valarray<Real> r0(N+1);
  for(int i=0; i<=N; i++)
    r0[i] = x[i];
  compute_c0r12(r0,c0, r1, r2);
  FILE *fp1 = fopen("bigger", "w");
  for(int i=0; i<=l; i++) { // can change 
    fprintf(fp1, "%23s\n",str(fft1.eval(&c0[0], pi*i/l),0) );
  }
}
