#include "grid/multigrid3D.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "view.h"
#include "lambda2.h"
#include "maxruntime.h"
#include "navier-stokes/perfs.h"
#include <math.h>

// Constant Parameter
#define density_ratio (1/1000.0)  // Density ratio 
#define viscosity_ratio (1.8e-5/1e-3)  // Viscosity ratio
#define eotvos 27.468           // Eotvos number
#define diameter 1.0          // pipe diameter
#define U1s_init 1.0

// Simulation Parameter 
scalar f0[];                    // volume fraction initially

// Manually input parameter
double max_grid_number;  // Grid number in axial direction
double pipe_length;      // pipe length (Shouldn't be integer - interface intersect the grid)
double forcing;          // forcing term
double h_L_D_init;             // Liquid height to diameter ratio
double U2s_init;          // Inititalize gas superficial velocity
double froude_liquid;
double RUNID;            // For restarting purpose (Look at Shell script)
double reynold_liquid;

// Function declaration
#define PI 3.14159265
#define POPEN(name, mode) fopen (name ".ppm", mode)    // mode corresponding to either read/write/append etc.
double circle_geometry(double h_L_D);   // Function to convert liquid height to area assuming flat interface

// Boundary Condition 
// No slip, no penetration BC on the embedded
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.r[embed] = dirichlet(0.);

int main(int argc, char * argv[]){

  // 7 inputs
  if (argc > 1)
    max_grid_number = atoi(argv[1]);     // grid number x direction
  if (argc > 2)
    pipe_length = atof(argv[2]);         // Pipe axial length
  if (argc > 3)
    forcing = atof(argv[3]);             // Forcing term, unit of g
  if (argc > 4)
    h_L_D_init = atof(argv[4]);          // Initial liquid height unit of D
  if (argc > 5)
    U2s_init = atof(argv[5]);            // Gas superficial velocity unit of sqrt(gD)
  if (argc > 6)
    froude_liquid = atof(argv[6]);            // Gas superficial velocity unit of sqrt(gD)
  if (argc > 7)
    reynold_liquid = atof(argv[7]);            // Gas superficial velocity unit of sqrt(gD)
  
    size(pipe_length);                 // setting the physical size if the x-dir (axial dir)
    dimensions (nx = pipe_length - 0.5 , ny = diameter , nz = diameter );    // domain size 
    init_grid (max_grid_number);               // grid number without refinement
    double center_point = 0.5 * pipe_length/(pipe_length - 0.5);
    origin (0, -center_point, -center_point);      // center point    

    DT = HUGE [0];
    rho1 = 1.0;                                 // Scaled density of phase 1 (Liquid)  
    rho2 = rho1 * density_ratio;                   // Scaled density of phase 2 (Gas)
    mu1 = 1.0/reynold_liquid;                                  // Scaled dynamyic viscosity of phase 1 (Liquid)  
    mu2 = mu1 * viscosity_ratio;                    // Scaled dynamic viscosity of phase 2 (Gas)
 
    f.sigma = 1./eotvos;                           // Surface tension coefficient sigma 
    periodic(right);               // Periodic BC at the axial direction
    run();
}

// Initial condition
event init (t = 0) {
  vertex scalar phi[];
  foreach_vertex()
    {
      phi[] = sq(diameter/2.0)-sq(y)-sq(z);
    }
  fractions(phi,cs,fs);
  fractions_cleanup(cs,fs);

  /*
  fraction(f0, y < (h_L_D_init - diameter/2) ?
	   -sq(y) - sq(z) + pow(diameter/2,2) :-1);    // initialize liquid holdup
  */

  fraction(f0, -sq(y) - sq(z) + pow(diameter/2,2));    // initialize liquid holdup
  
  double aL = 0.;
  foreach(reduction(+:aL)) {
    if (fabs(x) <= 0.5*Delta) {
      // area element in the y–z plane for this cell
      double dA = sq(Delta);

      // cell-open fraction from embedded boundary (1 if not using embed)
      double c = 1.;
      #ifdef EMBED
      c = cs[];
      #endif

      aL += f0[] * c * dA;
    }
  }
  double aT = aL/(PI*sq(0.5));
  
  double init_holdup = circle_geometry(h_L_D_init);
  fprintf(stderr,"At time t=0, liquid holdup: %g , actual %g\n",init_holdup, aT);

  boundary((scalar*){u});
  /*
  foreach() {
    f[] = f0[];                                      // Initialize volume fraction at initial time
    u.x[] =(cs[] > 0.) ?((U1s_init/aT)*f[] + (U2s_init/(1.0-aT))*(1-f[])):0.; // Initialize inlet BC at initial time
  }
  */
  foreach() {
    double r = sqrt(sq(y) + sq(z));  // radius from pipe center
    f[] = f0[];                      // initialize volume fraction
    u.x[] = (cs[] > 0.) ? 2.31658187369 * j0(4.81 * r) : 0.;  // velocity profile
  }

}

event acceleration (i++){
  face vector av = a;                        // acceleration term in vertical direction
  foreach_face(y){
    av.y[] -= sq(U1s_init/froude_liquid)/diameter;         // gravitational acceleration
  }
  foreach_face(x){
    av.x[] += forcing;                        // forcing term 
  }
}

double circle_geometry(double h_L_D){
  double phi_knot = acos(1.0 - 2.0 * h_L_D);
  double initial_holdup = (1.0/PI)*(phi_knot - 0.5*sin(2.0*phi_knot));
  return initial_holdup;
}

/**                                                                                                                                       
  Output video and field*/

event dump_snapshot (t+=1.0)
{
  cs.nodump = true;
  char name[80];
  sprintf (name, "Dump_file/Dump/dump-%05.4lf", t);
  dump (file = name);
}


event movies (t+=0.4) {
  static int frame = 0;  // persistent counter                                                                                             
  scalar omega[];
  vorticity(u, omega);
  clear();
  view(fov=0, tx=0.25, ty=0.25,
       theta = 3.14/2 + 0.5, phi = 0.4, psi = 0.,
       bg = {1,1,1}, samples = 4);
  //  squares ("u.x", linear = true, n = {0,0,1}, alpha = 0);                                                                             
  squares ("u.x", linear = true, n = {-1,0,0}, alpha = - pipe_length);
  squares ("u.x", linear = true, n = {-1,0,0}, alpha = -1);
  squares ("u.x", linear = true, n = {-1,0,0}, alpha = -3);
  
  //  squares ("omega", linear = true, n = {0,0,1}, alpha = 0);                                                                            
  cells   (n = {1,0,0}, alpha = -3.1415);
  draw_vof("f", color = "u.x");
  char s[80];
  sprintf(s, "t = %0.4f", t);
  draw_string(s, size = 30);
  char fname[80];
  sprintf(fname, "Dump_file/Images/3D-%04d.ppm", frame++);
  save(file = fname);
}

static FILE *csvfile = NULL;  // This is for .csv file tracking superficial velocity over time

event superficial (t += 0.2) {
  const int Nslice = 2;                     // planes k = 0,1,2
  const double L = pipe_length - 0.5;
  const double stp = L / Nslice;
  const double R = 0.5*diameter;         // Pipe radius
  const double R2 = R*R;                 // Square of radius
  const double cross_sec = PI*R2;         // for superficial velocities

  // one scalar per slice per quantity (legal reduction targets)
  double A_L0=0, A_L1=0, A_L2=0;
  double A_G0=0, A_G1=0, A_G2=0;
  double U_LS0=0, U_LS1=0, U_LS2=0;
  double U_GS0=0, U_GS1=0, U_GS2=0;

  double V_L=0, V_G=0; 
  
  // accumulate on x-faces so area is Δ^2 and velocity is face-centered
  foreach_face(x,
      reduction(+:A_L0) reduction(+:A_L1) reduction(+:A_L2)
      reduction(+:A_G0) reduction(+:A_G1) reduction(+:A_G2)
      reduction(+:U_LS0) reduction(+:U_LS1) reduction(+:U_LS2)
      reduction(+:U_GS0) reduction(+:U_GS1) reduction(+:U_GS2)) {

    // choose the nearest slice plane
    int k = (int) nearbyint(x / stp);
    if (k < 0 || k > Nslice) continue;
    double xk = k*stp;                           // xk: axial coordinate x of our calculated cross-section

    if (fabs(x - xk) > 0.51*Delta) continue;     // skip x if too far from specified xk 

    // open face area (Δ^2 * mask)
    double area_face = Delta*Delta;     // grid area
    double open = 1.0;                  // open:: cross-section of pipe [0,1]

    #ifdef USE_EMBED_MASK              // if you computed fs via: fractions(phi, cs, fs) for the pipe
      open = fs.x[];
    #else                               // geometric fallback: keep faces whose centers lie inside the circle
      open = (y*y + z*z <= R2) ? 1.0 : 0.0;
    #endif

    area_face *= open;                  // grid area multiply by pipe-cross section 
    if (area_face <= 0.) continue;      // skip outside the pipe (open = 0)

    // phase fraction at the face (average two adjacent cells)
    double f_face = clamp(0.5*(f[] + f[-1]), 0., 1.); // Average to get face value from cell center value (clamp:: ensure in the interval [0,1])
    
    if (k == 0) {
      A_L0  += area_face * f_face;
      A_G0  += area_face * (1. - f_face);
      U_LS0 += u.x[] * area_face * f_face/cross_sec;
      U_GS0 += u.x[] * area_face * (1. - f_face)/cross_sec;
    }
    else if (k == 1) {
      A_L1  += area_face * f_face;
      A_G1  += area_face * (1. - f_face);
      U_LS1 += u.x[] * area_face * f_face/cross_sec;
      U_GS1 += u.x[] * area_face * (1. - f_face)/cross_sec;
    }
    else { // k == 2
      A_L2  += area_face * f_face;
      A_G2  += area_face * (1. - f_face);
      U_LS2 += u.x[] * area_face * f_face/cross_sec;
      U_GS2 += u.x[] * area_face * (1. - f_face)/cross_sec;
    }
  }

  char fname[128];
  sprintf(fname, "Dump_file/Output_field/cell_dump_t%0.4f_pid%d.csv", t, pid());
  FILE *fp = fopen(fname, "w");
  if (!fp) { perror("Cannot open CSV dump file"); return 1; }
  fprintf(fp, "x,y,z,f,u.x,u.y,u.z\n");

  foreach(reduction(+:V_L) reduction(+:V_G)) {
    fprintf(fp, "%g,%g,%g,%g,%g,%g,%g\n",
            x, y, z, f[], u.x[], u.y[], u.z[]);
    
    // restrict to the axial region we care about (same as slices)
    if (x < -0.51*Delta || x > pipe_length + 0.51*Delta) continue;

    // open cell fraction inside the pipe
    double open_cell = 1.0;
    #ifdef USE_EMBED_MASK
      open_cell = cs[];                                // if you have embedded boundary, cs[] is the cell-volume fraction open to flow
    #else
     open_cell = (y*y + z*z <= R2) ? 1.0 : 0.0;        // geometric fallback by cell center (stair-step near wall)
    #endif
    if (open_cell <= 0.) continue;

    double dV = cube(Delta) * open_cell;       // effective cell volume contributing to the pipe

    // phase volumes
    double f_cell = clamp(f[], 0., 1.);
    V_L += f_cell * dV;
    V_G += (1. - f_cell) * dV;
  }

  // analytic pipe volume over [0,pipe_length] for sanity check
  double V_pipe = cross_sec * pipe_length;

  fprintf(stderr, "At time t:%g , i:%d \n",t,i);
  fprintf(stderr, "x=%.9g  A_L=%g  A_G=%g  U_LS=%g  U_GS=%g\n",
          0.0,     A_L0, A_G0, U_LS0, U_GS0);
  fprintf(stderr, "x=%.9g  A_L=%g  A_G=%g  U_LS=%g  U_GS=%g\n",
          stp,     A_L1, A_G1, U_LS1, U_GS1);
  fprintf(stderr, "x=%.9g  A_L=%g  A_G=%g  U_LS=%g  U_GS=%g\n",
          2*stp,   A_L2, A_G2, U_LS2, U_GS2);

  fprintf(stderr, "VOLUMES over [0,pipe_length]:  V_L=%g  V_G=%g  (V_L+V_G=%g,  analytic pipe=%g)\n",
          V_L, V_G, V_L + V_G, V_pipe);

  fclose(fp);

  // Open and write header only once
  if (!csvfile) {
    csvfile = fopen("Dump_file/diagnostics.csv", "w");
    if (!csvfile) {
      perror("Cannot open diagnostics.csv");
      exit(1);
    }
    fprintf(csvfile, "t,x,A_L,A_G,U_LS,U_GS,V_L,V_G,V_total,V_pipe\n");
  }

  // Append rows for each x-location
  fprintf(csvfile, "%.6f,%.9g,%g,%g,%g,%g,,,,\n",
          t, 0.0,   A_L0, A_G0, U_LS0, U_GS0);

  fprintf(csvfile, "%.6f,%.9g,%g,%g,%g,%g,,,,\n",
          t, stp,   A_L1, A_G1, U_LS1, U_GS1);

  fprintf(csvfile, "%.6f,%.9g,%g,%g,%g,%g,,,,\n",
          t, 2*stp, A_L2, A_G2, U_LS2, U_GS2);

  // volumes summary
  fprintf(csvfile, "%.6f,VOLUMES,,,,,%g,%g,%g,%g\n",
          t, V_L, V_G, V_L + V_G, V_pipe);
  fflush(csvfile); // flush to disk each step
}

/**                                                                                                                                      
  Ending simulation*/

event end(t = 6.00){
  if (csvfile) fclose(csvfile);
  fprintf(stderr,"Debugging bitch:t: %g\n",t);
  return 1;
}

