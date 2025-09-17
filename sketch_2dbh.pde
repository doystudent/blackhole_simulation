//Made by Matthew

double c = 299792458.0;
double G = 6.6743e-11;

float scaleFactor; // meters per pixel

class BlackHole {
  PVector position;
  double mass;
  double r_s;  // Schwarzschild radius in meters

  BlackHole(double mass, PVector position) {
    this.mass = mass;
    this.position = position;
    this.r_s = (2 * G * mass) / (c * c);
  }

  void draw() {
    fill(255);
    noStroke();
    float diameter = (float)(2 * r_s / scaleFactor); // meters to pixels
    ellipse(position.x, position.y, diameter, diameter);
  }
}

class Ray {
  // Position in meters (cartesian)
  double x, y;

  // Polar coordinates in meters, radians
  double r, phi;

  // Velocities in polar coords (meters per affine parameter unit)
  double dr, dphi;

  //conserved quantities
  double E, L;
  
  ArrayList<PVector> trail;

  Ray(double x, double y, PVector dir) {
    this.x = x;
    this.y = y;

    this.r = Math.sqrt(x*x + y*y);
    this.phi = Math.atan2(y, x);

    // Velocity vector (meters per affine parameter)
    double vx = dir.x * c;
    double vy = dir.y * c;

    // Project velocity onto polar coordinates
    this.dr = vx * Math.cos(phi) + vy * Math.sin(phi);
    this.dphi = (-vx * Math.sin(phi) + vy * Math.cos(phi)) / r;
    
    this.L = r*r * dphi;
    double f = 1.0 - bh.r_s/r;
    double dt_dλ = Math.sqrt( (dr*dr)/(f*f) + (r*r*dphi*dphi)/f );
    this.E = f * dt_dλ;
    
    trail = new ArrayList<PVector>();
  }

  void step(double r_s, double dλ) {
    // If inside event horizon, stop
    if (r < r_s) {
      return;
    }
    
    rk4Step(this, dλ, bh.r_s);

    // Convert back to cartesian (meters)
    x = r * Math.cos(phi);
    y = r * Math.sin(phi);

    // Add pixel coords to trail
    float px = (float)(x / scaleFactor);
    float py = (float)(y / scaleFactor);
    trail.add(new PVector(px, py));

    // Limit trail length
    if (trail.size() > 1000) {
      trail.remove(0);
    }
  }

  void draw() {
    stroke(255, 255, 0);
    strokeWeight(2);
    // Draw current point
    point((float)(x / scaleFactor), (float)(y / scaleFactor));

    noFill();
    stroke(255, 255, 100, 100);
    beginShape();
    for (PVector p : trail) {
      vertex(p.x, p.y);
    }
    endShape();
  }
}

Ray copyWithState(Ray original, double[] state) {
  Ray copy = new Ray(0, 0, new PVector(1, 0)); // dummy values for direction
  copy.r = state[0];
  copy.phi = state[1];
  copy.dr = state[2];
  copy.dphi = state[3];
  copy.E = original.E;
  copy.L = original.L;
  return copy;
}

void geodesicRHS(Ray ray, double[] rhs, double rs) {
  double r    = ray.r;
  double dr   = ray.dr;
  double dphi = ray.dphi;
  double E    = ray.E;
  
  double f = 1.0 - rs/r;
  
  rhs[0] = dr;
  rhs[1] = dphi;
  
  double dt_dλ = E / f;
  
  rhs[2] = 
  - (rs/(2*r*r)) * f * (dt_dλ*dt_dλ)
  + (rs/(2*r*r*f)) * (dr*dr)
  + (r - rs) * (dphi*dphi);
  
  rhs[3] = -2.0 * dr * dphi / r;
}

void addState(double[] a, double[] b, double factor, double[] out) {
  for (int i = 0; i < 4; i++) {
    out[i] = a[i] + b[i] * factor;
  }   
}

void rk4Step(Ray ray, double dλ, double rs) {
  double[] y0 = { ray.r, ray.phi, ray.dr, ray.dphi };
  double[] k1 = new double[4], 
           k2 = new double[4], 
           k3 = new double[4], 
           k4 = new double[4], 
           temp = new double[4];

  // k1
  geodesicRHS(ray, k1, rs);

  // k2
  addState(y0, k1, dλ / 2.0, temp);
  Ray r2 = copyWithState(ray, temp);
  geodesicRHS(r2, k2, rs);

  // k3
  addState(y0, k2, dλ / 2.0, temp);
  Ray r3 = copyWithState(ray, temp);
  geodesicRHS(r3, k3, rs);

  // k4
  addState(y0, k3, dλ, temp);
  Ray r4 = copyWithState(ray, temp);
  geodesicRHS(r4, k4, rs);

  // Final update
  ray.r    += (dλ / 6.0) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
  ray.phi  += (dλ / 6.0) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
  ray.dr   += (dλ / 6.0) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);
  ray.dphi += (dλ / 6.0) * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]);
}


// Globals
BlackHole bh;
ArrayList<Ray> rays;

void setup() {
  size(800, 800);

  // Sagittarius A* mass (kg)
  bh = new BlackHole(8.55e36, new PVector(0, 0));

  // Scale so black hole diameter is 100 pixels
  scaleFactor = (float)((2 * bh.r_s) / 100.0);

  rays = new ArrayList<Ray>();

  // Start rays at x = -500 r_s meters left of BH, with slight vertical spread, moving right (+x)
  for (int i = 0; i < 100; i++) {
    double startX = -5 * bh.r_s;
    double startY = (-500 + i * 10) * bh.r_s * 0.01; // small vertical spread in meters
    rays.add(new Ray(startX, startY, new PVector(1, 0))); // moving right
  }
}

void draw() {
  background(0);
  translate(width/2, height/2); // center origin

  bh.draw();

  for (Ray ray : rays) {
    ray.step(bh.r_s, 1); // pass Schwarzschild radius in meters
    ray.draw();
  }
}
