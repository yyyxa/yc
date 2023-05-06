//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%     Centrifugal Spinning via DEM Simulation       %%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include <stdio.h>
#include <math.h>     // Must use "gcc *.c -lm" as math.h is included.
#include <stdlib.h>
#include <time.h>
#define Pi 3.1415926   

int main(void)
{  
clock_t t_start, t_end; 
long double t_spent;   // time spent during running
t_start=clock();

int k, i, ui, di, n, nt, Zone, Km, j; 
int material;
int  Ntot, N_in, count; // No. atoms
int Nt;     // Dump data every N timsteps 
int Zone_max=5000;    // Numer of Zone for visualization in Telpot
int Nmax=2000;   // maximum No. of atom 
//int Order[Nmax+1];
int I_end, I_in;
int n_spt[Nmax+1];

double  rpm; 
double  Rdisk;    // radius of disk, [m]
double  a0;       // nozzle radius, [m]
double  lint;     // initial distance between beads, [m]
double  dt;       // dimensionless timestep  
double  Q;        // Flow rate, [ml/min]

double density;  //Density of liquid, [kg/m^3]
double Gamma;    //Surface tention, [N/m]
double mu = 0.0;       // viscosity [Pa-s]
double Lamda = 0.0;    //Stress relaxation time [s]
double Betas = 0.0, Betap;
double Omega;    // rotating rate, [1/s]
double T_air;    // room temperature, [K]
double rho_air[Nmax+1];  // air Density, [kg/m^3]      
double nu_air[Nmax+1];   // air Kinematic Viscosity, [m^2/s]  
double kair[Nmax+1]; // air thermal conductivity, [J/(m s K)]
double Cp_air;       // Heat capacity of air, []
//double Va[3+1];  // air velocity, [m/s]
double Vair[3+1];// air velocity, [--]
double U0;       // Initial velocity at the nozzle exit, [m/s]
double m0;       // Mass of bead, [kg]
double tau0;     // characteristic stress, [N/m^2]
double t0;       // characteristic time, [s]
double T0;       // characteristic temperature, [K]
double t;        // dimensionless time  
double DtDt;     // =(dt)^2
double tinit;
double D, thickness;
 
double  R0;    // dim-less Rdisk
double  epson;  // dim-less lint
double  Re;
double  We;
double  Fr;
double  Rb;
double  De;
double  Re_air;
double  Na;
double  Hf;
double  Pr;
double  Hconv;
  
double  xxf; // Dim-less aerodynamic-drag coefficient (friction)
double  xxp; // Dim-less aerodynamic-drag coefficient (pressure)
double  Xd ; // aerodynamic-drag coefficient (form drag)
double  Fg[3+1] ; //  Gravity force

double  Sp11[Nmax+1]; // Radial Stress
double  dSp11;
double  Sp22[Nmax+1]; // Azimuthal Stress
double  dSp22;
double  Sp12[Nmax+1];
double  bmax;
double  tr_Sp[Nmax+1], tr_QQ[Nmax+1];
double  Stress[Nmax+1], Stress1[Nmax+1], Stress2[Nmax+1]; 
double  X[Nmax+1],     Y[Nmax+1],     Z[Nmax+1];  //     coordinates
double  X_old[Nmax+1], Y_old[Nmax+1], Z_old[Nmax+1];  // old coordinates
double  a[Nmax+1];  // dim-less radius of beads
double  L[Nmax+1],   L_old[Nmax+1];  // dim-less element length
double  dL[Nmax+1];
double  mass[Nmax+1];  // dim-less mass
double  VX[Nmax+1],  VY[Nmax+1],    VZ[Nmax+1], Vr[Nmax+1];  // velocity
double  length;
double  r[Nmax+1];
double  Tin;

double ex[Nmax+1], ey[Nmax+1], ez[Nmax+1];  //the unit vector of element
double rs[3+1], rss[3+1], kc[3+1], vu[3+1];
double term, term1, term2, term3;  // temporary variables used to store some trivial values...
double Xnew, Ynew, Znew; 

double cry[Nmax+1], K, nA, dcry[Nmax+1], Kmax, cry_inf, fa, TD, Tmax, C, delalpha_i, Cop;
double T[Nmax+1], hf, h, Cp, dT[Nmax+1], delH; 
double Tref, fv, alpha; 

double Fvx[Nmax+1], Fvy[Nmax+1], Fvz[Nmax+1];
double Fstx[Nmax+1], Fsty[Nmax+1], Fstz[Nmax+1];
double NormVt, NormVn,Vre[3+1], Vt[3+1], Vn[3+1], Fv1[Nmax+1], Fv2[Nmax+1];
double Fax[Nmax+1], Fay[Nmax+1], Faz[Nmax+1], Fa_ui[3+1], Fa_di[3+1]; 
double AceX[Nmax+1], AceY[Nmax+1], AceZ[Nmax+1], Ace[Nmax+1], Ace_r[Nmax+1];
double Lem, Src[Nmax+1];
double cf, Xf, XP;
double cp(double);

double Ex[Nmax+1],Ey[Nmax+1], Ez[Nmax+1];
double Fex[Nmax+1], Fey[Nmax+1], Fez[Nmax+1];
double q[Nmax+1], qold[Nmax+1];
double I0, ke, Ge, cond, Ea, q0, V0, Cc, h0, H, rij;

double amps, Omega1;
double temp4[Nmax+1];
double attx[1], atty[1];
int collector;

//char Infile[]="/Users/mounicajyothidivvela/Desktop:FS3.in";
FILE *StreamIn;
char Outfile[30];        FILE *StreamOut;
char logfile[30];        FILE *logf;
char tab[30];

// Open files:
//Outfile = CSmelt_test.dat;  
StreamOut=fopen("CSmelt_test.dat","w"); // create a file for Tecplot

//logfile = CSmelt_test.log;  
logf=fopen("CSmelt_test.log","w"); // create a file for reviewing result
  
fprintf (StreamOut,"TITLE = \"Electrospinning Spinning via DEM Simulation\" \n");
fprintf (StreamOut,"VARIABLES = \"id\" \"ContourLength\" \"x\" \"y\" \"z\" \"Radius\" \"Stress\" \"Fez[i]\" \"Fvz[i]\" \"Fstz[i]\" \"Faz[i]\" \"q[i]\" \"T[i]\" \"AceZ[i]\" \n");  


//%%%%%   MATERIAL PARAMETERS   %%%%%
printf("Melt electrospinning of PLA \n"); 

		// PLA 
        Tref =180.0+273.15; 
		mu = 1320; 
        Lamda= 0.05; 
		alpha= 0.015;
        Betas=0.35;
		delH = 9060;
		hf = -93760.0;
		cry_inf = 0.2952;

		a0 = 0.42e-3;
		h0 = 0.05;		
		Tin = 255.0+273.15;
		T_air = 80.0+273.15;
		Q = 0.05e-3/60.0; // Flowrate
		V0 = 10000.0; // Voltage


   lint = 0.5*a0;
   dt = 0.00001;
   
// Characteristic temperature [K]
	T0 = fabs(Tref*Tref/-delH);
	Tref = Tref/T0;
	Tin = Tin/T0;
	
	amps = 0.1;
	Omega1 = 1.0e10;

printf(    "\n***************\n");
printf(     "Input Data: \n");
printf(    "****************\n");
fprintf(logf,"\n**************\n");
fprintf(logf,"Input Data:  \n");
fprintf(logf,"**************\n");
   
// Density [kg/m^3]
        density=1240.0;
// Surface tention [N/m]
        Gamma= 0.0435;
// Ratio of polymer viscosity to zero-shear-rate viscosity     
        Betap=1.0-Betas;
// Applied current [A]
	I0=1.8e-9;
// Coulomb's constant [N.m^2/C^2]
	ke=8.987e9;  
// Electrical Conductivity 
	cond=1e-10;        
// Heat capacity [J/kg-K]
	Cp = 1800;
	Kmax = 0.55;
	Tmax = (65.0+273.15)/T0;
	TD = (60.0+273.15)/T0;

	delalpha_i = 0.0468;
	nA = 1;
	Cop = 1e-8; //[m^2/N]
	C = 1;
     
// OPERATIONAL PARAMETERS   %%%%%  
        printf("\n V0 = %g [V] \n",V0);
        fprintf(logf,"\n RPM = %g \n",V0);
// Disk radius, [m]         
        printf(" h0 = %g [m] \n",h0);
        fprintf(logf," h0 = %g [m] \n",h0);
// Nozzle radius, [m]
        printf(" a0 = %g [m] \n",a0);          
        fprintf(logf," a0 = %g [m] \n",a0);        
// Initial distance between beads, [m]     
        printf(" initial distance between beads = %g [m]\n",lint);
        fprintf(logf," initial distance between beads = %g [m]\n",lint);
// Flow Rate, [ml/min];     
        printf(" Flow Rate = %g [m3/s]\n",Q);
        fprintf(logf," Flow Rate = %g [m3/s]\n",Q);
        U0=Q/(Pi*a0*a0);


  printf(    "\n***************************\n");
  printf(     "Characteristic Values: \n");
  printf(    "***************************\n");
  fprintf(logf,"\n***************************\n");
  fprintf(logf,"Characteristic Values:  \n");
  fprintf(logf,"***************************\n");
   printf(     "  initial velocity at nozzle (characteristic velocity) = %g m/s\n",U0);
  fprintf(logf,"  initial velocity at nozzle (characteristic velocity) = %g m/s\n",U0);
  
// Room temperature [K]
        T_air= T_air/T0;     
// Air velocity [--]
        Vair[1]=0.0/U0;   //x-direction
        Vair[2]=0.0/U0;   //y-direction
        Vair[3]=-10000000000.0;   //z-direction   
// Heat capacity of air
	Cp_air=1012.0; 
        
//%%%%%%%%%% Characterisitc scales %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  m0=density*Pi*a0*a0*lint;    // Mass of bead, [kg]
   printf(    "  Mass of bead =%5.3e kg\n",m0);
  fprintf(logf,"  Mass of bead =%5.3e kg\n",m0);
  
  tau0=mu*U0/a0;               // stress, [N/m^2]
   printf(    "  characteristic Stress = %g N/m^2\n",tau0);
  fprintf(logf,"  characteristic Stress = %g N/m^2\n",tau0);
  
  t0=a0/U0;                    // time, [s]
   printf(    "  characteristic time = %g s\n",t0);
  fprintf(logf,"  characteristic time = %g s\n",t0);

   printf(    "  characteristic temperature = %g K\n\n",T0);
  fprintf(logf,"  characteristic temperature = %g K\n\n",T0);
  
// Characteristic charge [C]
  q0=8.854e-12*2*V0/(a0*log(1+(4*h0/a0)))*2*3.14*a0*lint;
   printf(    "  characteristic charge = %g C\n\n",q0);
  fprintf(logf,"  characteristic charge = %g C\n\n",q0);  

//%%%%%%%%%% Simulation Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  DtDt=dt*dt;            // (dt)^2
  Nt=(int)(0.1/dt);     // Dump data every N timsteps 
//%%%%%%%% Dimensionless Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  printf(    "***************************\n");
  printf(    "Dimensionless Parameters: \n");
  printf(    "***************************\n");
  fprintf(logf,"***************************\n");
  fprintf(logf,"Dimensionless Parameters: \n");
  fprintf(logf,"***************************\n");
  
  epson=lint/a0;
  H = h0/a0;

  length = 2.0*epson;
    printf(    "  length = %g \n",length);
    fprintf(logf,"  length = %g \n",length);

     printf(    "  epson=lint/a0 = %g \n",epson);
    fprintf(logf,"  epson=lint/a0 = %g \n",epson);
  Re=density*U0*a0/mu ;
     printf(    "  Re=density*U0*a0/mu = %g \n",Re);
    fprintf(logf,"  Re=density*U0*a0/mu = %g \n",Re);    
  We=density*U0*U0*a0/Gamma;
     printf(    "  We=density*U0*U0*a0/Gamma = %g \n",We);
    fprintf(logf,"  We=density*U0*U0*a0/Gamma = %g \n",We);
  Fr=U0/sqrt(9.8*a0);
     printf(    "  Fr=U0/sqrt(9.8*a0) = %g \n",Fr);
    fprintf(logf,"  Fr=U0/sqrt(9.8*a0) = %g \n",Fr);
  Rb=1.0/(Omega*t0);
     printf(    "  Rb=1.0/(Omega*t0) = %g \n",Rb);
    fprintf(logf,"  Rb=1.0/(Omega*t0) = %g \n",Rb);
  De=Lamda/t0 ;
     printf(    "  De=Lamda/t0 = %g \n",De);
    fprintf(logf,"  De=Lamda/t0 = %g \n",De);
  Na=tau0/(2.0*density*Cp*T0);
     printf(    "  Na=tau0/(2.0*density*Cp*T0) = %g \n",Na);
    fprintf(logf,"  Na=tau0/(2.0*density*Cp*T0) = %g \n",Na);
  Hf=hf*cry_inf/(Cp*T0);
     printf(    "  Hf=hf*cry_inf/(Cp*T0) = %g \n\n",Hf);
    fprintf(logf,"  Hf=hf*cry_inf/(Cp*T0) = %g \n\n",Hf);
  Ea=q0*V0/(m0*U0*U0);
     printf(    "  Ea=q0*V0/(m0*U0*U0) = %g \n",Ea);
    fprintf(logf,"  Ea=q0*V0/(m0*U0*U0) = %g \n",Ea);
  Cc=cond*a0*a0*m0*U0*Pi/(q0*q0);
     printf(     "  cond*a0*a0*m0*U0*Pi/(q0*q0) = %g \n",Cc);
    fprintf(logf,"  cond*a0*a0*m0*U0*Pi/(q0*q0) = %g \n",Cc);
  Ge=ke*q0*q0/(m0*U0*U0*a0);
     printf(    "  Ge=ke*q0*q0/(m0*U0*U0*a0) = %g \n \n",Ge);
    fprintf(logf,"  Ge=ke*q0*q0/(m0*U0*U0*a0) = %g \n \n",Ge);        

  Fg[1] =  0.0;
  Fg[2] =  0.0;
  Fg[3] = -1.0/Fr/Fr; // Gravity force

//%%%%%%% Initial Conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //N=2;  // number of beads : initial value

int col = 0;
int layer = 0;
collector = 1;
attx[1] = 1.0;
atty[1] = -1.0;

  Ntot = 1;

for (i=1;i<=Ntot;i++)
{
  Sp11[i]=0.0;
  Sp12[i]=0.0;
  Sp22[i]=0.0;
  Stress[i]=0.0;

  X[i] = 0.0; 
  Y[i] = 0.0;    
  Z[i] = (i-Ntot)*epson;
  
  a[i]=((Z[i])*1.07+1.1)/1.1;        
  L[i]=epson;  
  L_old[i]= epson;      
  mass[i]=1.0*a[i]*a[i];
  
  q[i]=100.0;  
  qold[i]=q[i];  

  T[i]=Tin;
  dT[i]=0.0;
  cry[i]=0.0;
  dcry[i]=0.0;

  rho_air[i]=704.0/(T0*(T_air+T_air));
  term = T0*(T_air+T_air);
  kair[i] = 1.03e-4*pow(term,0.866);
  nu_air[i] = 1.517e-9*pow(term,2.5)/(term+240.0);
  
  VX[i] = 0.0;    
  VY[i] = 0.0;           
  VZ[i] = -1.0;	
  
  Vr[i]=1.0;
  n_spt[i]=0.0;

  AceX[i] = 0.0; 
  AceY[i] = 0.0; 
  AceZ[i] = 0.0;
  
  X_old[i]=X[i]-VX[i]*dt;
  Y_old[i]=Y[i]-VY[i]*dt;  
  Z_old[i]=Z[i]-VZ[i]*dt;
  
  r[i] = R0-epson;
  
  Stress[i] = 0.0;
  
  temp4[i] = 0.0;
  
}

  I_end = Ntot;


//%%%%%%%%%%%%%%%=== START LOOP =====%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=0.0;
nt=0;
Zone=0;

while (Zone<=Zone_max && Ntot<Nmax)
{     
    t=t+dt;   // the dimensionless time lapse
    nt=nt+1;

for (i=col+1;i<=Ntot;i++)
{
if(L[i]>100000.0*epson && n_spt[i]<5)
{
n_spt[i]=n_spt[i]+1;
Ntot = Ntot+1;
I_end = Ntot;

for(n=I_end;n>i+1;n--)
{

    X[n]=X[n-1];
    Y[n]=Y[n-1]; 
    Z[n]=Z[n-1];  

    n_spt[n] = n_spt[n-1];   

    VX[n]=VX[n-1];
    VY[n]=VY[n-1];
    VZ[n]=VZ[n-1];

    AceX[n] = AceX[n-1];
    AceY[n] = AceY[n-1];
    AceZ[n] = AceZ[n-1];
    
    T[n]=T[n-1];
  dT[n]=dT[n-1];
  cry[n]=cry[n-1];
  dcry[n]=dcry[n-1];
  
  q[n]=q[n-1];
  qold[n]=q[n-1];

  rho_air[n]=704.0/(T0*(T_air+T[n]));
  term = T0*(T_air+T[n]);
  kair[n] = 1.03e-4*pow(term,0.866);
  nu_air[n] = 1.517e-9*pow(term,2.5)/(term+240.0);

    X_old[n]=X_old[n-1];
    Y_old[n]=Y_old[n-1];
    Z_old[n]=Z_old[n-1];

    L[n]=L[n-1]; 
    L_old[n]=L_old[n-1]; 
    mass[n] = mass[n-1];
    Sp11[n] = Sp11[n-1];
    Sp12[n] = Sp12[n-1];
    Stress[n]=Stress[n-1];

    Sp22[n] = Sp22[n-1];       
    a[n] = a[n-1];
}
L[i]=L[i]*0.5;
L_old[i]=L_old[i]*0.5;
mass[i]=0.5*mass[i];

L[i+1]=L[i];
L_old[i+1]=L_old[i];
a[i+1]=a[i];
mass[i+1] = mass[i];
n_spt[i+1] = n_spt[i];

q[i+1] = q[i];
qold[i+1] = qold[i];

Sp11[i+1]=Sp11[i];
Sp12[i+1]=Sp12[i];
Sp22[i+1]=Sp22[i];
Stress[i+1]=Stress[i];

    X[i+1]=0.5*(X[i]+X[i+2]);
    Y[i+1]=0.5*(Y[i]+Y[i+2]); 
    Z[i+1]=0.5*(Z[i]+Z[i+2]);  

    VX[i+1]=0.5*(VX[i]+VX[i+2]);
    VY[i+1]=0.5*(VY[i]+VY[i+2]);
    VZ[i+1]=0.5*(VZ[i]+VZ[i+2]);
    
      T[i+1]=T[i];
  dT[i+1]=dT[i];
  cry[i+1]=cry[i];
  dcry[i+1]=dcry[i];

  rho_air[i+1]=704.0/(T0*(T_air+T[i+1]));
  term = T0*(T_air+T[i+1]);
  kair[i+1] = 1.03e-4*pow(term,0.866);
  nu_air[i+1] = 1.517e-9*pow(term,2.5)/(term+240.0);

    AceX[i+1] = 0.5*(AceX[i]+AceX[i+2]);
    AceY[i+1] = 0.5*(AceY[i]+AceY[i+2]);
    AceZ[i+1] = 0.5*(AceZ[i]+AceZ[i+2]);

    X_old[i+1]=0.5*(X_old[i]+X_old[i+2]);
    Y_old[i+1]=0.5*(Y_old[i]+Y_old[i+2]);
    Z_old[i+1]=0.5*(Z_old[i]+Z_old[i+2]);

}
}
   
    //% Insert bead inside the nozzle    
    if (Z[I_end]<=-epson)
{
    I_end = I_end+1;

    Ntot = Ntot+1;
    
    i = I_end;
    
    //% position of newly-inserted bead:
    X[I_end]=0.0;
    Y[I_end]=0.0; 
    Z[I_end]=0.0;  

    VX[I_end]=0.0;
    VY[I_end]=0.0;
    VZ[I_end]=-1.0;

    n_spt[I_end]=0;

    AceX[I_end] = 0.0;
    AceY[I_end] = 0.0;
    AceZ[I_end] = 0.0;  

    X_old[I_end]=X[I_end]-VX[I_end]*dt;
    Y_old[I_end]=Y[I_end]-VY[I_end]*dt;
    Z_old[I_end]=Z[I_end]-VZ[I_end]*dt;

    L[I_end]=epson; 
    L_old[I_end]=epson; 
    mass[I_end]=1.0;
    Sp11[I_end] = 0.0;
    Sp12[I_end] = 0.0;
    Stress[I_end]=0.0;

    T[I_end]=Tin;
    dT[I_end]=0.0;
    cry[I_end]=0.0;
    dcry[I_end]=0.0;
    
	q[I_end]=100.0;
    qold[I_end]=q[I_end];

    rho_air[i]=704.0/(T0*(T_air+T_air));
    term = T0*(T_air+T_air);
    kair[i] = 1.03e-4*pow(term,0.866);
    nu_air[i] = 1.517e-9*pow(term,2.5)/(term+240.0);
  
    Sp12[I_end] = -2.0*Betap*0.0;
    Sp11[I_end] = 4.0*Betap*De*0.0;
    Sp22[I_end]=0.0;        
    Stress[I_end] = 0.0; 
    
    a[I_end]=1.0;
    
    temp4[I_end] = 0.0;

    }
    

    // Calculate the location and the velocity for all beads
    for (i=col+1;i<=Ntot;i++)
    {    VX[i]= (X[i]-X_old[i])/dt;             
         VY[i]= (Y[i]-Y_old[i])/dt;
         VZ[i]= (Z[i]-Z_old[i])/dt;       
    }
    
	for (i=col+1; i<=Ntot-1; i++) {
	Vre[1] = Vair[1]-0.5*(VX[i]+VX[ui]);
	Vre[2] = Vair[2]-0.5*(VY[i]+VY[ui]);
	Vre[3] = Vair[3]-0.5*(VZ[i]+VZ[ui]);
	
	term = Vre[1]*ex[i]+Vre[2]*ey[i]+Vre[3]*ez[i];

	Vt[1] = term*ex[i];
	Vt[2] = term*ey[i];
	Vt[3] = term*ez[i];

	Vn[1] = Vre[1]-Vt[1];
	Vn[2] = Vre[2]-Vt[2];
	Vn[3] = Vre[3]-Vt[3];
	NormVt = sqrt(Vt[1]*Vt[1]+Vt[2]*Vt[2]+Vt[3]*Vt[3]);
	NormVn = sqrt(Vn[1]*Vn[1]+Vn[2]*Vn[2]+Vn[3]*Vn[3]);
	
	if (NormVt == 0) {
	term = 1.0;
	}
	else {
	term = 1.0+(64.0*NormVn*NormVn/(NormVt*NormVt));
	}

	Re_air = (2.0*U0*a0/nu_air[i]);
	Pr = rho_air[i]*Cp_air*nu_air[i]/kair[i];    
	Hconv = 0.495*pow(Re_air,0.33)*pow((NormVt/(a[i]*a[i])),0.33)*kair[i]*sqrt(Pr)*pow(term,0.167)/(density*Cp*U0*a0);
	dT[i] = Na*Stress[i]*dL[i]-Hf*dcry[i]-Hconv*(T[i]-T_air)/a[i];	
	T[i] = T[i]+dt*dT[i];
	rho_air[i] = 704/(T0*(T_air+T_air));
	kair[i] = 1.03e-4*pow((T0*(T_air+T_air)),0.866);
	nu_air[i] = 1.517e-9*pow((T0*(T_air+T_air)),2.5)/(T0*(T_air+T_air)+240);
}  

/*	for (i=1; i<=I_end; i++) {
	term = log(1/(1-cry[i]));
	term1 = (nA-1)/nA;
	fa = -Cop*Sp11[i]*tau0/delalpha_i;
	K = Kmax*t0*exp(-4*log(2)*((T[i]-Tmax)/TD)*((T[i]-Tmax)/TD)+C*fa*fa);
	if (cry[i] <=0.99) {
	dcry[i] = 0.0;//nA*K*(1-cry[i])*pow(term,term1);
	cry[i] = 0.0;//cry[i]+dt*dcry[i];
	}
	else {
	dcry[i] = 0.0;
	cry[i] = 0.0;//0.9999;
	}
}*/	

    // Calculate the length of each segment
   
    for (n=col+1;n<=I_end-1;n++)
    {     
        i=n;
        ui=i+1;
        
        //length of segment (use the upstream bead)

        L[i]=sqrt((X[i]-X[ui])*(X[i]-X[ui])+(Y[i]-Y[ui])*(Y[i]-Y[ui])+(Z[i]-Z[ui])*(Z[i]-Z[ui]));

        ex[i]=(X[i]-X[ui])/L[i]; 
        ey[i]=(Y[i]-Y[ui])/L[i]; // unit vector of upstream element   
        ez[i]=(Z[i]-Z[ui])/L[i]; 
             
        a[i]=sqrt(mass[i]*epson/L[i]);  // Calculate radius (no evaporation)
        dL[i]=(L[i]-L_old[i])/L[i]/dt;  // strain rate (upstream)
	    L_old[i]=L[i]; 

	fv = exp(delH*(1.0/T[i]-1.0/Tref)/T0+4.0*(cry[i])*(cry[i]));

	dSp11 = (2.0*Betap*fv*dL[i]-Sp11[i]-(alpha*De*Tref*Sp11[i]*Sp11[i])/(Betap*T[i]))/(De*fv*Tref/T[i])+Sp11[i]*dT[i]/T[i]+2.0*dL[i]*Sp11[i];
	dSp22 = (-Betap*fv*dL[i]-Sp22[i]-(alpha*De*Tref*Sp22[i]*Sp22[i])/(Betap*T[i]))/(De*fv*Tref/T[i])-dL[i]*Sp22[i];
	Sp11[i] = Sp11[i]+dSp11*dt;
	Sp22[i] = Sp22[i]+dSp22*dt; 

	Stress[i] = Sp11[i]-Sp22[i]+3.0*Betas*fv*dL[i];
     }  // END for n..


    for (n=col+1;n<=Ntot;n++)
{
       i=n;
       di = n-1;	
       ui = n+1;

{
       term = sqrt((H-layer*a[1]+Z[i])*(H-layer*a[1]+Z[i])+(attx[1]-X[i])*(attx[1]-X[i])+(atty[1]-Y[i])*(atty[1]-Y[i]));
       Ex[i]=Ea*(attx[1]-X[i])/(term*term)-layer*2.0*Ge*0.00000001*q[1]/(a[1]*L[1])*(attx[1]-X[i])/(term*term);
       Ey[i]=Ea*(atty[1]-Y[i])/(term*term)-layer*2.0*Ge*0.00000001*q[1]/(a[1]*L[1])*(atty[1]-Y[i])/(term*term);
       Ez[i]=-(Ea*(H-layer*a[1]+Z[i]))/(term*term)+layer*2.0*Ge*0.00000001*q[1]/(a[1]*L[1])*(H-layer*a[1]+Z[i])/(term*term);            	
}

       if(temp4[i]==0.0)
{
		Ex[i] = Ex[i];//+amps*2*Ea/((1.0-2.0*Z[i]-Z[i]*Z[i]/H)*log(1+(8.0*H/2.0)))*cos(Omega1*t)*H/(H+Z[i]);
		Ey[i] = Ey[i];//+amps*2*Ea/((1.0-2.0*Z[i]-Z[i]*Z[i]/H)*log(1+(8.0*H/2.0)))*sin(Omega1*t)*H/(H+Z[i]);		
		temp4[i] = 0.0;		
}

	//Ez[i] = Ez[i]+temp4[i];	
	   //if (Z[i] <= -0.6) {
/*       for (j=col+1;j<=Ntot;j++)
{
       if((j!=i) && (i<Ntot))
{
    rij=sqrt((X[i]-X[j])*(X[i]-X[j])+(Y[i]-Y[j])*(Y[i]-Y[j])+(Z[i]-Z[j])*(Z[i]-Z[j]));

	term1=(X[i]-X[j]);
	term2=(Y[i]-Y[j]);
	term3=(Z[i]-Z[j]);
	
    Ex[i]=Ex[i]+Ge*q[j]*term1/(rij*rij*rij);
    Ey[i]=Ey[i]+Ge*q[j]*term2/(rij*rij*rij);
    Ez[i]=Ez[i]+Ge*q[j]*term3/(rij*rij*rij);
} 
}
//}
//}
//	else {
/*	if (i!=1 && i!=Ntot) {
       Ex[i]=Ex[i]+Ge*2.0*0.01*log(H)*(q[i]*ex[i]/L[i]-q[di]*ex[di]/L[di])/(L[i]*a[i]);//+2.1*0.001*log(H)*(Ex[ui]*a[ui]*a[ui]-2.0*Ex[i]*a[i]*a[i]+Ex[di]*a[di]*a[di])/(L[i]*L[i]*2.0);
       Ey[i]=Ey[i]+Ge*2.0*0.01*log(H)*(q[i]*ey[i]/L[i]-q[di]*ey[di]/L[di])/(L[i]*a[i]);//+2.1*0.001*log(H)*(Ey[ui]*a[ui]*a[ui]-2.0*Ey[i]*a[i]*a[i]+Ey[di]*a[di]*a[di])/(L[i]*L[i]*2.0);            
       Ez[i]=Ez[i]+Ge*2.0*0.01*log(H)*(q[i]*ez[i]/L[i]-q[di]*ez[di]/L[di])/(L[i]*a[i]);//+2.1*0.001*log(H)*(Ez[ui]*a[ui]*a[ui]-2.0*Ez[i]*a[i]*a[i]+Ez[di]*a[di]*a[di])/(L[i]*L[i]*2.0);
    }
    {
       Ex[i]=Ex[i]+Ge*2.0*0.01*log(H)*(q[i]*ex[i]/L[i])/(L[i]*a[i]);
       Ey[i]=Ey[i]+Ge*2.0*0.01*log(H)*(q[i]*ey[i]/L[i])/(L[i]*a[i]);
       Ez[i]=Ez[i]+Ge*2.0*0.01*log(H)*(q[i]*ez[i]/L[i])/(L[i]*a[i]);	
    }
 //   }*/
}

      for (n=col+1;n<=Ntot;n++)
{
	i=n;
	di=i-1;

       if((i!=col+1))
{
       term1 = ex[i]*Ex[i]+ey[i]*Ey[i]+ez[i]*Ez[i];
       term2 = ex[di]*Ex[di]+ey[di]*Ey[di]+ez[di]*Ez[di];
       q[i]=qold[i]+dt*((Cc*a[i]*a[i]*term1)-(Cc*a[di]*a[di]*term2));
}

       if((i==col+1) || i==Ntot)
{
       term1 = ex[i]*Ex[i]+ey[i]*Ey[i]+ez[i]*Ez[i];
       q[i] = qold[i]+dt*((Cc*a[i]*a[i]*term1));
}

/*       if((i==N-1 && i!=1))
{
       term2 = ex[di]*Ex[di]+ey[di]*Ey[di]+ez[di]*Ez[di];
       q[i]=qbef[i]+dt*(-(Cc*a[di]*a[di]*term2));
}*/

	qold[i]=q[i];
}	 
 
    // Calculate the applied force on the internal particles
   
    for (n=col+2;n<=I_end;n++)
    {
     i=n;
     ui=i+1;
     di=i-1;       
        
     //Surface tension
       Fstx[i]=2.0*(a[di]*ex[di]-a[i]*ex[i]); 
       Fsty[i]=2.0*(a[di]*ey[di]-a[i]*ey[i]);
       Fstz[i]=2.0*(a[di]*ez[di]-a[i]*ez[i]);
     
     if (i < Ntot) {   
     term=(L[i]+L[di]);
       rs[1] =(X[di]-X[ui])/(term);
       rs[2] =(Y[di]-Y[ui])/(term);
       rs[3] =(Z[di]-Z[ui])/(term);
     
     term=(L[i]*L[di]);
       rss[1]=(X[di]-2.0*X[i]+X[ui])/(term);
       rss[2]=(Y[di]-2.0*Y[i]+Y[ui])/(term);
       rss[3]=(Z[di]-2.0*Z[i]+Z[ui])/(term);
     
       kc[1]= rs[2]*rss[3]-rs[3]*rss[2];
       kc[2]= rs[3]*rss[1]-rs[1]*rss[3];  //cross product
       kc[3]= rs[1]*rss[2]-rs[2]*rss[1];
       kc[0]=sqrt(kc[1]*kc[1] +kc[2]*kc[2] +kc[3]*kc[3]);  //curvature
       
       vu[1]=ex[di]-ex[i];
       vu[2]=ey[di]-ey[i];
       vu[3]=ez[di]-ez[i];
       
       vu[0]=sqrt(vu[1]*vu[1] +vu[2]*vu[2] +vu[3]*vu[3]);
       
   	   if (vu[0] == 0) {
   	   vu[1] = 0.0;
   	   vu[2] = 0.0;
   	   vu[3] = 0.0;
   	   }
   	   else {
       vu[1]=vu[1]/vu[0];
       vu[2]=vu[2]/vu[0];
       vu[3]=vu[3]/vu[0];
       }
   
     term=(a[di]+a[i])*(a[di]+a[i])*kc[0];
       Fstx[i]=Fstx[i]+ 0.25*(term)*vu[1];
       Fsty[i]=Fsty[i]+ 0.25*(term)*vu[2];
       Fstz[i]=Fstz[i]+ 0.25*(term)*vu[3];
   }
     term=(epson*We);
       Fstx[i]=Fstx[i]/(term);
       Fsty[i]=Fsty[i]/(term);
       Fstz[i]=Fstz[i]/(term);
          
     //Viscoelastic force
     term =(a[di]*a[di]* Stress[di]);
     term1=( a[i]*a[i] *  Stress[i]);

     term2=(epson*Re);

       Fvx[i]=((term)*ex[di]-(term1)*ex[i])/(term2);
       Fvy[i]=((term)*ey[di]-(term1)*ey[i])/(term2);
       Fvz[i]=((term)*ez[di]-(term1)*ez[i])/(term2);
      
      if (i < Ntot) {  
     // Aerodynamic drag 
     // --------relative velocity for upstream element
       
       Vre[1]=-0.5*(VX[i]+VX[ui]);
       Vre[2]=-0.5*(VY[i]+VY[ui]); 
       Vre[3]=-0.5*(VZ[i]+VZ[ui]);
                 
     term=Vre[1]*ex[i]+Vre[2]*ey[i]+Vre[3]*ez[i];   // inner product         
       Vt[1]=(term)*ex[i];
       Vt[2]=(term)*ey[i];   // tangential relative velocity
       Vt[3]=(term)*ez[i];
          
       Vn[1]=Vre[1]-Vt[1];
       Vn[2]=Vre[2]-Vt[2];   // Normal relative velocity
       Vn[3]=Vre[3]-Vt[3];
          
       NormVt=sqrt(Vt[1]*Vt[1] +Vt[2]*Vt[2] +Vt[3]*Vt[3]);
       NormVn=sqrt(Vn[1]*Vn[1] +Vn[2]*Vn[2] +Vn[3]*Vn[3]);
        
     term =(2.0*U0*a0*a[i]/nu_air[i]);   term1=(term*NormVt);  term2=(term*NormVn);  
       cf=0.78*pow(term1,-0.61);
       xxf=0.5*(rho_air[i]/density/epson);  // Dim-less aerodynamic-drag coefficient (friction)
       Xf=cf*xxf;
       xxp=xxf/Pi;                // Dim-less aerodynamic-drag coefficient (pressure)
       XP=cp(term2)*xxp;
        
     term=(a[i]*L[i]); 
       Fa_ui[1] = ( Xf*NormVt*Vt[1]+XP*NormVn*Vn[1] )*(term);
       Fa_ui[2] = ( Xf*NormVt*Vt[2]+XP*NormVn*Vn[2] )*(term);
       Fa_ui[3] = ( Xf*NormVt*Vt[3]+XP*NormVn*Vn[3] )*(term);

     //-----------relative velocity for downstream element

       Vre[1]=-0.5*(VX[i]+VX[di]);
       Vre[2]=-0.5*(VY[i]+VY[di]); 
       Vre[3]=-0.5*(VZ[i]+VZ[di]);
          
     term=Vre[1]*ex[di] +Vre[2]*ey[di] +Vre[3]*ez[di];   // inner product         
       Vt[1]=(term)*ex[di];
       Vt[2]=(term)*ey[di];   // tangential relative velocity
       Vt[3]=(term)*ez[di];
          
       Vn[1]=Vre[1]-Vt[1];
       Vn[2]=Vre[2]-Vt[2];   // Normal relative velocity
       Vn[3]=Vre[3]-Vt[3];
          
       NormVt=sqrt(Vt[1]*Vt[1] +Vt[2]*Vt[2] +Vt[3]*Vt[3]);
       NormVn=sqrt(Vn[1]*Vn[1] +Vn[2]*Vn[2] +Vn[3]*Vn[3]);
       
     term =(2.0*U0*a0*a[di]/nu_air[di]);   term1=(term*NormVt);  term2=(term*NormVn);  
       cf=0.78*pow(term1,-0.61);
       xxf=0.5*(rho_air[i]/density/epson);  // Dim-less aerodynamic-drag coefficient (friction)
       Xf=cf*xxf;
       xxp=xxf/Pi;                // Dim-less aerodynamic-drag coefficient (pressure)
       XP=cp(term2)*xxp;
          
     term=(a[di]*L[di]); 
       Fa_di[1] = ( Xf*NormVt*Vt[1]+XP*NormVn*Vn[1] )*(term);
       Fa_di[2] = ( Xf*NormVt*Vt[2]+XP*NormVn*Vn[2] )*(term);
       Fa_di[3] = ( Xf*NormVt*Vt[3]+XP*NormVn*Vn[3] )*(term);
          
       Fax[i]=(Fa_ui[1]+Fa_di[1])/2.0;
       Fay[i]=(Fa_ui[2]+Fa_di[2])/2.0;
       Faz[i]=(Fa_ui[3]+Fa_di[3])/2.0;
       }
       
       else {
       Fax[i]=0.0;
       Fay[i]=0.0;
       Faz[i]=0.0;       
       }
       
       Fex[i] = q[i]*Ex[i];//+Ge*q[i]*a[i]/(L[i]*a[i]*a[i])*(q[i]*ex[i]/(a[i]*L[i])-q[di]*ex[di]/(a[di]*L[di]));
       Fey[i] = q[i]*Ey[i];//+Ge*q[i]*a[i]/(L[i]*a[i]*a[i])*(q[i]*ey[i]/(a[i]*L[i])-q[di]*ey[di]/(a[di]*L[di]));
       Fez[i] = q[i]*Ez[i];//+Ge*q[i]*a[i]/(L[i]*a[i]*a[i])*(q[i]*ez[i]/(a[i]*L[i])-q[di]*ez[di]/(a[di]*L[di]));
          
    // Total applied force
       AceX[i]= Fg[1]+(0.0*Fvx[i]+0.0*Fax[i]+0.0*Fstx[i]+Fex[i])/mass[i]; 
       AceY[i]= Fg[2]+(0.0*Fvy[i]+0.0*Fay[i]+0.0*Fsty[i]+Fey[i])/mass[i];
       AceZ[i]= Fg[3]+(0.0*Fvz[i]+0.0*Faz[i]+0.0*Fstz[i]+Fez[i])/mass[i]; 
              
       if (AceZ[i] > 0.0) {
       		AceZ[i] = 0.0;
       }
       
       if (i==Ntot) {
       		AceX[Ntot] = 0.0;
       		AceY[Ntot] = 0.0;
       		AceZ[Ntot] = 0.001*AceZ[Ntot];
       } 
       
    /*    if (Z[i] < -1.0) {
    	AceZ[i] = Fg[3] + (Fstz[i]+100.0*Fez[i])/mass[i];       	
    }   
       
    if (Z[i] < -1.5) {
    	AceZ[i] = Fg[3] + (Fstz[i]+1000.0*Fez[i])/mass[i];       	
    }*/
       
       Ace[i] = sqrt(AceX[i]*AceX[i]+AceY[i]*AceY[i]+AceZ[i]*AceZ[i]);
     }
   
   {   
    //  Calculate the applied force on the front particle
    // ---- axial capillary force
    if(col==0) {
    i=col+1;
    term=-a[i]*2.0/(We*epson); 
       Fstx[i]=(term)*ex[i]; 
       Fsty[i]=(term)*ey[i];
       Fstz[i]=(term)*ez[i];
       
    term= -a[i]*a[i]*Stress[i];
    term= term/(epson*Re);  
       Fvx[i]=(term)*ex[i];   
       Fvy[i]=(term)*ey[i];
       Fvz[i]=(term)*ez[i];

	if (Ntot > 1) {
    
       Vre[1]=-0.5*(VX[1]+VX[2]);
       Vre[2]=-0.5*(VY[1]+VY[2]);
       Vre[3]=-0.5*(VZ[1]+VZ[2]);
    
    term= Vre[1]*ex[1] + Vre[2]*ey[1] + Vre[3]*ez[1];
       Vt[1]=(term)*ex[i];
       Vt[2]=(term)*ey[i];
       Vt[3]=(term)*ez[i];
       
       Vn[1]=Vre[1]-Vt[1];
       Vn[2]=Vre[2]-Vt[2];   
       Vn[3]=Vre[3]-Vt[3];
       
       NormVt=sqrt(Vt[1]*Vt[1]+Vt[2]*Vt[2]+Vt[3]*Vt[3]);
       NormVn=sqrt(Vn[1]*Vn[1]+Vn[2]*Vn[2]+Vn[3]*Vn[3]);
       
       xxf=0.5*(rho_air[i]/density/epson);  // Dim-less aerodynamic-drag coefficient (friction)
       xxp=xxf/Pi;                // Dim-less aerodynamic-drag coefficient (pressure)
       XP=cp(Re_air*NormVn*a[i])*xxp;
    term=XP*NormVn*a[i]*L[i];
       Fax[i]=(term)*Vn[1];
       Fay[i]=(term)*Vn[2];
       Faz[i]=(term)*Vn[3];
       
    term=Vt[1]*ex[1]+Vt[2]*ey[1]+Vt[3]*ez[1];   
       if (term<0) {
	Xd=0.1*(rho_air[i]/density/epson/Pi)/2.0;  // aerodynamic-drag coefficient (form drag)
           term1=-Xd*NormVt*a[i]*a[i];}
       else
          {cf=0.78*pow(Re_air*NormVt*a[i],-0.61);
	   xxf=0.5*(rho_air[i]/density/epson);  // Dim-less aerodynamic-drag coefficient (friction)
           Xf=cf*xxf;
           term1= Xf*NormVt*a[i]*L[i];}        
       
       Fax[i] = Fax[i] + (term1)*Vt[1];
       Fay[i] = Fay[i] + (term1)*Vt[2];                          
       Faz[i] = Faz[i] + (term1)*Vt[3];      
       }
       
       else {
       Fax[i] = 0.0;
       Fay[i] = 0.0;
       Faz[i] = 0.0;
       } 
       
       Fex[i] = q[i]*Ex[i];
       Fey[i] = q[i]*Ey[i];
       Fez[i] = q[i]*Ez[i];      
              
    AceX[i] = Fg[1] + (0.0*Fvx[i]+0.0*Fax[i]+0.0*Fstx[i]+Fex[i])/mass[i];
    AceY[i] = Fg[2] + (0.0*Fvy[i]+0.0*Fay[i]+0.0*Fsty[i]+Fey[i])/mass[i];
    AceZ[i] = Fg[3] + (0.0*Fvz[i]+0.0*Faz[i]+0.0*Fstz[i]+Fez[i])/mass[i];      

    Ace[i] = sqrt(AceX[i]*AceX[i]+AceY[i]*AceY[i]+AceZ[i]*AceZ[i]);
    
    if (AceZ[i] > 0.0) {
    	AceZ[i] = 0.0;
    }
    }
    
    /*    if (Z[1] < -1.0) {
    	AceZ[1] = Fg[3] + (Fstz[1]+100.0*Fez[1])/mass[1];       	
    }
    
    if (Z[1] < -1.5) {
    	AceZ[1] = Fg[3] + (Fstz[1]+1000.0*Fez[1])/mass[1];       	
    }*/
    }
    
	if (col==1 || (X[col+1] >= 0.95 && collector==1 && col > 1)) {
		attx[1] = 1.0;
		atty[1] = 1.0;
		collector = 2;
		printf("collector = %d \n",collector);
		if (col!=1) {
		layer = layer+1;
}
	}
	
	if ((Y[col+1]>=0.95)  && collector==2) {
		attx[1] = -1.0;
		atty[1] = 1.0;
		collector = 3;
		printf("collector = %d \n",collector);		
	}
	
	if ((X[col+1]<=-0.95)  && collector==3) {
		attx[1] = -1.0;
		atty[1] = -1.0;
		collector = 4;
		printf("collector = %d \n",collector);		
	}	
	
	if (Y[col+1]<=-0.95 && collector==4) {
		attx[1] = 1.0;
		atty[1] = -1.0;
		collector = 1;
		printf("collector = %d \n",collector);		
	}	
	
    
    // Update the position of flying beads
    for (n=col+1;n<=I_end;n++)
    {    i=n;
        Xnew = 2.0*X[i] -X_old[i] + DtDt*AceX[i]; 
        Ynew = 2.0*Y[i] -Y_old[i] + DtDt*AceY[i];
        Znew = 2.0*Z[i] -Z_old[i] + DtDt*AceZ[i];
        
        X_old[i] = X[i];
        Y_old[i] = Y[i];
        Z_old[i] = Z[i];

        X[i] = Xnew;
        Y[i] = Ynew;
        Z[i] = Znew;
    }
    
    n = col;
    for (i=n+1;i<=I_end;i++)
    {
    	if (-Z[i] >=(H-layer*a[1])) {
    		col = i;
    	}
    }
    
    /// Dump data
    if  (nt==Nt )
    {
        Zone=Zone+1;
        nt=0;
 
        //if (I_end>0)
{       Lem=0;
        for (n=Ntot;n>=1;n=n-1)
        {   i=n;      
            Lem=Lem+ L[i];  
            Src[n]=Lem;
        }      

        printf("Zone=%d, N=%d, t=%5.3e[s], a(1)=%5.3e, Src=%5.3e \n",Zone, Ntot,t*t0,a[1],Src[1]);
        fprintf(logf,"Zone=%d, N=%d, t=%5.3e[s], a(1)=%5.3e, Src=%5.3e \n",Zone, Ntot,t*t0,a[1],Src[1]);
        fprintf(StreamOut,"ZONE I=%d, F=POINT \n", Ntot); 
}
                
        for (n=1;n<=Ntot;n++)            
        {
            i=n;            
            fprintf (StreamOut,"%d %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %d %d \n", i,Src[n],X[i],Y[i],Z[i],a[i],Stress[i], Fez[i], Fvz[i], Fstz[i], Faz[i], q[i], T[i], AceZ[i], col, layer);
        }
                  
    } // End of dumping loop
  
    
}//END for while
    


t_end=clock();
t_spent=(t_end-t_start);
t_spent=t_spent/CLOCKS_PER_SEC;
printf("\n Total CPU time elasped= %Lf [s] \n",t_spent);
fprintf(logf,"\n Total CPU time elasped= %Lf [s] \n",t_spent);



fclose(StreamOut);
fclose(StreamIn);
fclose(logf);

return (0);
}//END for main





double cp(double Re)
{

double out, Re_start, Re_end, cp_start, cp_end, slope, power;

if (Re<0.1)
    out=100.0;
else if (Re<1.0 && Re>=0.1)
 {  Re_start=0.1;
    Re_end=1.0;
    cp_start=100.0;
    cp_end=10.0;
    slope=(log10(cp_end)-log10(cp_start))/(log10(Re_end)-log10(Re_start));
    power=log10(cp_start) + slope*(log10(Re)-log10(Re_start));
    out=pow(10.0,power);
 }   
else if (Re>=1.0 && Re<10.0)
 {  Re_start=1.0;
    Re_end=10.0;
    cp_start=10.0;
    cp_end=4.0;
    slope=(log10(cp_end)-log10(cp_start))/(log10(Re_end)-log10(Re_start));
    power=log10(cp_start) + slope*(log10(Re)-log10(Re_start));
    out=pow(10.0,power);
 }   
else if (Re>=10.0 && Re<100.0)
 {  Re_start=10.0;
    Re_end=100.0;
    cp_start=4.0;
    cp_end=2.0;
    slope=(log10(cp_end)-log10(cp_start))/(log10(Re_end)-log10(Re_start));
    power=log10(cp_start) + slope*(log10(Re)-log10(Re_start));
    out=pow(10.0,power);
 }   
else if (Re>=100.0 && Re<1000.0)
 {  Re_start=100.0;
    Re_end=1000.0;
    cp_start=2.0;
    cp_end=1.2;  
    slope=(log10(cp_end)-log10(cp_start))/(log10(Re_end)-log10(Re_start));
    power=log10(cp_start) + slope*(log10(Re)-log10(Re_start));
    out=pow(10.0,power);
 }   
else if (Re>=1.0e3 && Re<1.0e4)
 {  Re_start=1.0e3;
    Re_end=1.0e4;
    cp_start=1.2;
    cp_end=1.5;
    slope=(log10(cp_end)-log10(cp_start))/(log10(Re_end)-log10(Re_start));
    power=log10(cp_start) + slope*(log10(Re)-log10(Re_start));
    out=pow(10.0,power);
 }   
else if (Re>=1.0e4 && Re<1.0e5)
 {  Re_start=1.0e4;
    Re_end=1.0e5;
    cp_start=1.5;
    cp_end=2.0;
    slope=(log10(cp_end)-log10(cp_start))/(log10(Re_end)-log10(Re_start));
    power=log10(cp_start) + slope*(log10(Re)-log10(Re_start));
    out=pow(10.0,power);
 }  
else if (Re>=1.0e5 && Re<2.0e5)
 {  Re_start=1.0e5;
    Re_end=2.0e5;
    cp_start=2.0;
    cp_end=2.2;
    slope=(log10(cp_end)-log10(cp_start))/(log10(Re_end)-log10(Re_start));
    power=log10(cp_start) + slope*(log10(Re)-log10(Re_start));
    out=pow(10.0,power);
 }   
else if (Re>=2.0e5 && Re<3.0e5)
 {  Re_start=2.0e5;
    Re_end=3.0e5;
    cp_start=2.2;
    cp_end=0.38;
    slope=(log10(cp_end)-log10(cp_start))/(log10(Re_end)-log10(Re_start));
    power=log10(cp_start) + slope*(log10(Re)-log10(Re_start));
    out=pow(10.0,power);
 }  
else if (Re>=3.0e5 )
 {  Re_start=3.0e5;
    Re_end=1.0e6;
    cp_start=0.38;
    cp_end=0.8;
    slope=(log10(cp_end)-log10(cp_start))/(log10(Re_end)-log10(Re_start));
    power=log10(cp_start) + slope*(log10(Re)-log10(Re_start));
    out=pow(10.0,power);
 }

else
    out=1.0;



return (out);

}
          
