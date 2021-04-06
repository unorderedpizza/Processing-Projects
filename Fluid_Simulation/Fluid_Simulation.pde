
 /* Fluid Simulation
 
 ---

 Introduction
 
  - In this program, a numerical method of fluid simulation is presented, 
    which utilizes the incompressible Navier-Stokes equations:
    
    nabla · u = 0                                                        (1)
    
      |du              |
    ρ |-- + u · nabla u| = - nabla p + µ nabla^2 u + ρ g                 (2)
      |dt              |
    
  - Equation 2 is first satisfied through the direct advection of momentum
    between cells.
  - Equation 1 is then satisfied in the final velocity field through 
    Helmholtz decomposition, where the curl-free velocity field obtained 
    from the intermediate velocity field is subtracted from the intermediate
    velocity field to obtain the final divergence-free velocity field. */

/* --------------------------------------------------------------------------

 Types */

class vector_2d_t {
  
  // variables
  
  double x, y;
  
  // functions
  
  vector_2d_t() {
    x = 0.0;
    y = 0.0;
  }
  vector_2d_t(vector_2d_t a) {
    x = a.x;
    y = a.y;
  }
  vector_2d_t(double a, double b) {
    x = a;
    y = b;
  }
  
  vector_2d_t add(vector_2d_t a) {
    vector_2d_t temp = new vector_2d_t(x + a.x, y + a.y);
    return temp;
  }
  vector_2d_t scale(double a) {
    vector_2d_t temp = new vector_2d_t(a * x, a * y);
    return temp;
  }
}

/* --------------------------------------------------------------------------

 Variables */
 
double RHO = 1.0;
double MU = 0.06;
double G = 0.0;

double CELL_SIZE = 0.1;
double FRAME_DT = 1.0 / 60.0;
boolean ADAPTIVE_STEP = true;
double TIME_STEP_MULT = 0.7;
int GAUSS_S_N_STEPS = 3;

double SQ_CELL_SIZE = sq(CELL_SIZE);
double INV_CELL_SIZE = 1.0 / CELL_SIZE;
double INV_SQ_CELL_S = 1.0 / sq(CELL_SIZE);
double VISC_ACC_MULT = MU / RHO;
vector_2d_t G_ACC = new vector_2d_t(0, G);

vector_2d_t u[]; // discretized velocities
vector_2d_t temp_u[]; // temporary storage
double a[]; // an arbitrary quantity
double temp_a[]; // temporary storage
double div[]; // discretized divergence values
double phi[]; // scalar potential values for helmholtz decomposition

/* --------------------------------------------------------------------------

 Initialization */
 
void setup() {
  size(256, 128);
  background(0);
  
  u = new vector_2d_t[width * height];
  temp_u = new vector_2d_t[width * height];
  for (int i = 0; i < width * height; i ++) {
    u[i] = new vector_2d_t();
    temp_u[i] = new vector_2d_t();
  }
  
  a = new double[width * height];
  temp_a = new double[width * height];
  div = new double[width * height];
  phi = new double[width * height];
}

/* --------------------------------------------------------------------------

 Main Loop */
 
void draw() {
  
  // ensure covergence
  
  double time_step;
  if (ADAPTIVE_STEP) {
    double max_c = 0.000001;
    for (int i = 0; i < width; i ++) {
      for (int j = 0; j < height; j ++) {
        if (abs(u[i * height + j].x) > max_c) {
          max_c = abs(u[i * height + j].x);
        }
        if (abs(u[i * height + j].y) > max_c) {
          max_c = abs(u[i * height + j].y);
        }
      }
    }
    time_step = TIME_STEP_MULT * CELL_SIZE / max_c;
  } else {
    time_step = FRAME_DT;
  }
  
  // loop until everything is complete
  
  int n_steps = 0;
  if (FRAME_DT < time_step) {
    time_step = FRAME_DT;
    n_steps = 1;
  } else {
    n_steps = (int) (FRAME_DT / time_step);
  }
  for (int i = 0; i < n_steps; i ++) {
    
    // reset temporary storage
    
    for (int j = 0; j < width; j ++) {
      for (int k = 0; k < height; k ++) {
        temp_u[j * height + k] = new vector_2d_t();
        temp_a[j * height + k] = 0.0;
      }
    }
    
    // apply boundary conditions
    
    if (mousePressed) {
      int mouse_x = constrain(mouseX, 0, width - 1);
      int mouse_y = constrain(mouseY, 0, height - 1);
      
      temp_u[mouse_x * height + mouse_y].x = -32.0;
      temp_u[mouse_x * height + mouse_y].y = -32.0;
      
      if (mouseButton == LEFT) {
        temp_a[mouse_x * height + mouse_y] = time_step * 1000.0;
      }
    }
    for (int j = 0; j < width; j ++) {
      temp_u[j * height] = new vector_2d_t();
      temp_u[j * height + height - 1] = new vector_2d_t();
    }
    for (int j = 1; j < height - 1; j ++) {
      temp_u[j] = new vector_2d_t();
      temp_u[(width - 1) * height + j] = new vector_2d_t();
    }
    
    for (int j = 1; j < width - 1; j ++) {
      for (int k = 1; k < height - 1; k ++) {
        
        // local attributes
        
        vector_2d_t cell_loc = new vector_2d_t(j * CELL_SIZE, k * CELL_SIZE);
        vector_2d_t visc_acc = (new vector_2d_t(
          (u[(j + 1) * height + k].x + u[(j - 1) * height + k].x 
          + u[j * height + k + 1].x + u[j * height + k - 1].x
          - 4 * u[j * height + k].x) * INV_SQ_CELL_S, 
          (u[(j + 1) * height + k].y + u[(j - 1) * height + k].y 
          + u[j * height + k + 1].y + u[j * height + k - 1].y 
          - 4 * u[j * height + k].y) * INV_SQ_CELL_S)).scale(VISC_ACC_MULT);
        vector_2d_t u_prime = u[j * height + k].add(visc_acc.add(G_ACC)
          .scale(time_step));
        
        // bilinear interpolation
        
        vector_2d_t dest_loc = cell_loc.add(u[j * height + k]
          .scale(time_step));
        vector_2d_t dest_ind = dest_loc.scale(INV_CELL_SIZE);
        
        int dest_ind_x = (int) dest_ind.x;
        int dest_ind_y = (int) dest_ind.y;
        
        double fr_v = dest_ind.y - dest_ind_y;
        double fr_h = dest_ind.x - dest_ind_x;
        
        double c0 = fr_v * fr_h;
        double c1 = fr_v * (1 - fr_h);
        double c2 = (1 - fr_v) * (1 - fr_h);
        double c3 = (1 - fr_v) * fr_h;
        
        // advect properties between cells
        
        dest_ind_x = constrain(dest_ind_x, 0, width - 2);
        dest_ind_y = constrain(dest_ind_y, 0, height - 2);
        temp_u[(dest_ind_x + 1) * height + dest_ind_y + 1] = 
          temp_u[(dest_ind_x + 1) * height + dest_ind_y + 1]
          .add(u_prime.scale(c0));
        temp_u[dest_ind_x * height + dest_ind_y + 1] = 
          temp_u[dest_ind_x * height + dest_ind_y + 1]
          .add(u_prime.scale(c1));
        temp_u[dest_ind_x * height + dest_ind_y] = 
          temp_u[dest_ind_x * height + dest_ind_y].add(u_prime.scale(c2));
        temp_u[(dest_ind_x + 1) * height + dest_ind_y] = 
          temp_u[(dest_ind_x + 1) * height + dest_ind_y]
          .add(u_prime.scale(c3));
        
        temp_a[(dest_ind_x + 1) * height + dest_ind_y + 1] = 
          temp_a[(dest_ind_x + 1) * height + dest_ind_y + 1] + c0 
          * a[j * height + k];
        temp_a[dest_ind_x * height + dest_ind_y + 1] = 
          temp_a[dest_ind_x * height + dest_ind_y + 1] + c1 
          * a[j * height + k];
        temp_a[dest_ind_x * height + dest_ind_y] = 
          temp_a[dest_ind_x * height + dest_ind_y] + c2 
          * a[j * height + k];
        temp_a[(dest_ind_x + 1) * height + dest_ind_y] = 
          temp_a[(dest_ind_x + 1) * height + dest_ind_y] + c3 
          * a[j * height + k];
      }
    }
    
    // remove divergence
    
    remove_div();
    
    // copy to actual arrays
    
    for (int j = 1; j < width - 1; j ++) {
      for (int k = 1; k < height - 1; k ++) {
        u[j * height + k] = new vector_2d_t(temp_u[j * height + k]);
        a[j * height + k] = temp_a[j * height + k];
      }
    }
  }
  
  // graphics
  
  loadPixels();
  color c;
  for (int i = 1; i < width - 1; i ++) {
    for (int j = 1; j < height - 1; j ++) {
      c = color(constrain((int) (255 * a[i * height + j]), 0, 255));
      pixels[j * width + i] = c;
    } 
  }
  updatePixels();
}

/* --------------------------------------------------------------------------

 Functions */
 
void remove_div() { // helmholtz decomposition
  
  // compute divergence values
  
  for (int i = 1; i < width - 1; i ++) {
    for (int j = 1; j < height - 1; j ++) {
      div[i * height + j] = (temp_u[(i + 1) * height + j].x
        - temp_u[(i - 1) * height + j].x + temp_u[i * height + (j + 1)].y
        - temp_u[i * height + (j - 1)].y) * 0.5 * INV_CELL_SIZE;
    }
  }
  
  // use the gauss-seidel method to solve for phi values
  // values are slightly offsetted for better convergence
  
  for (int i = 0; i < GAUSS_S_N_STEPS; i ++) {
    
    // match divergence
    
    for (int j = 1; j < width - 1; j ++) {
      for (int k = 1; k < height - 1; k ++) {
        phi[j * height + k] = - (SQ_CELL_SIZE * div[j * height + k]
          - phi[(j + 1) * height + k] - phi[(j - 1) * height + k]
          - phi[j * height + k + 1] - phi[(j + 1) * height + k - 1]) * 0.249;
      }
    }
    
    // match vertical velocity at boundaries
    
    for (int j = 1; j < width - 1; j ++) {
      phi[j * height] = - 0.99 * (CELL_SIZE * temp_u[j * height].y 
        - phi[j * height + 1]);
      phi[j * height + height - 1] = 0.99 * CELL_SIZE 
        * temp_u[j * height + height - 1].y 
        + phi[j * height + height - 2];
    }
    
    // match horizontal velocity at boundaries
    
    for (int j = 1; j < height - 1; j ++) {
      phi[j] = - 0.99 * (CELL_SIZE * temp_u[j].x - phi[height + j]);
      phi[(width - 1) * height + j] = 0.99 * CELL_SIZE 
        * temp_u[(width - 1) * height + j].x 
        + phi[(width - 2) * height + j];
    }
  }
  
  // correct intermediate velocity field
  
  for (int j = 1; j < width - 1; j ++) {
    for (int k = 1; k < height - 1; k ++) {
      temp_u[j * height + k].x -= (phi[(j + 1) * height + k] 
        - phi[(j - 1) * height + k]) * 0.5 * INV_CELL_SIZE;
      temp_u[j * height + k].y -= (phi[j * height + k + 1] 
        - phi[j * height + k - 1]) * 0.5 * INV_CELL_SIZE;
    }
  }
}

double abs(double a) {
  return (a < 0)? - a : a;
}
double sq(double a) {
  return a * a;
}
