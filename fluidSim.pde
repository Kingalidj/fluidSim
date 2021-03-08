simulation fluid;
float kernelRadius = 40;
float timeStep = 2;
float stepsPerFrame = 20;

void setup() {
	size(1200, 1200, P2D);
	colorMode(HSB);

	fluid = new simulation(timeStep, kernelRadius);
	for (float x = width * 0.4; x<(width * 0.6); x+=kernelRadius*0.25)
		for (float y = height * 0.2; y<(height * 0.8); y+=kernelRadius*0.25)
			fluid.particles.add(new particle(x, y, false));
	createBorders();
}

void draw() {
	background(255);
	if (mousePressed && mouseButton == LEFT)addDroplet(mouseX, mouseY, 10, 20);
	if (mousePressed && mouseButton == RIGHT)fluid.particles.add(new particle(mouseX, mouseY, true));
	for (int i = 0; i < stepsPerFrame; i++)fluid.step();
	fluid.show();
}

void addDroplet(float x, float y, float n, float r) {
	for (int i = 0; i < n; i++)
		fluid.particles.add(new particle(x + random(-r, r), y + random(-r, r), false));
}

void createBorders() {
	for (float x=0; x < width; x+= kernelRadius * 0.1) {
		fluid.particles.add(new particle(x, height, true));
		fluid.particles.add(new particle(x, 0, true));
	}
	for (float y=0; y < height; y+= kernelRadius * 0.1) {
		fluid.particles.add(new particle(0, y, true));
		fluid.particles.add(new particle(width, y, true));
	}
}

class simulation {
	ArrayList<particle> particles;
	grid Grid;
	PVector g;
	float dt;
	float h;
	float k;
	float nk;
	float rho0;

	simulation(float dt, float h) {
		particles = new ArrayList<particle>();
		Grid = new grid(width, height, h);
		g = new PVector(0, 0.01);
		g.mult(dt);
		this.dt = dt;
		this.h = h;
		k = 0.008;
		nk = 0.01;
		rho0 = 10;
	}

	void step() {
		for (particle p : particles) {
			if (!p.rigid) {
				p.vel.add(g);
			}
		}

		for (particle p : particles) {
			p.ppos = p.pos.copy();
			p.pos.add(PVector.mult(p.vel, dt));
		}

		Grid.clearGrid();
		for (particle p : particles) {
			Grid.addParticle(p);
		}

		for (particle p : particles) {
			p.getNeighbors(Grid);
			p.rho = 0;
			p.nrho = 0;

			for (particle pNeighbor : p.neighbors) {
				if (pNeighbor != p) {
					PVector rIJVec = PVector.sub(pNeighbor.pos, p.pos);
					float rIJ = rIJVec.mag();
					float q = rIJ / h;
					if (q < 1) {
						float temp = 1 - q;
						float powOfTemp = temp * temp;
						p.rho += powOfTemp;
						powOfTemp *= temp;
						p.nrho += powOfTemp;
					}
				}
			}

			p.pressure = k * (p.rho - rho0);
			p.npressure = nk * p.nrho;

			p.dx.mult(0);
			for (particle pNeighbor : p.neighbors) {
				PVector rIJVec = PVector.sub(pNeighbor.pos, p.pos);
				float rIJ = rIJVec.mag();
				float q = rIJ / h;
				if (q < 1) {
					rIJVec.normalize();
					PVector d = rIJVec;
					d.mult(dt * dt * (p.pressure * (1 - q) + p.npressure * (1 - q) * (1 - q)));
					d.mult(0.5);
					if (!pNeighbor.rigid)pNeighbor.pos.add(d);
					p.dx.sub(d);
				}
			}
			if (!p.rigid)p.pos.add(p.dx);
		}

		for (particle p : particles) {
			if (p.pos.x < 0)p.pos.x = 0;
			if (p.pos.x > width)p.pos.x = width;
			if (p.pos.y < 0)p.pos.y = 0;
			if (p.pos.y > height)p.pos.y = height;
		}

		for (particle p : particles) {
			p.vel = PVector.sub(p.pos, p.ppos);
			p.vel.div(dt);
		}
	}

	void show() {
		for (particle p : particles)p.show();

		//int res = 80;
		//float [][] avgPressure = new float[res][res];
		//for (int i = 0; i < res; i++)
		//  for (int j = 0; j < res; j++) {
		//    float x = (float) i / res * width;
		//    float y = (float) j / res * height;
		//    float d = (float) width / res;
		//    int counter = 0;
		//    for (particle p : particles) {
		//      if (dist(p.pos.x, p.pos.y, x + d/2, y + d/2) <= d) {
		//        avgPressure[i][j] += p.pressure; 
		//        counter++;
		//      }
		//    }
		//    fill(255);
		//    if (counter != 0) {
		//      avgPressure[i][j] /= counter;
		//      fill(150 - avgPressure[i][j] * 2000, 255, 200);
		//    }
		//    rect(x, y, d, d);
		//  }
	}
}

class particle {
	PVector pos;
	PVector vel;
	PVector ppos;
	PVector dx;
	ArrayList <particle> neighbors; 

	float rho;
	float nrho;
	float pressure;
	float npressure;
	boolean rigid;

	particle(float x, float y, boolean r) {
		pos = new PVector(x, y);
		vel = new PVector();
		neighbors = new ArrayList<particle>();
		rho = 0;
		nrho = 0;
		dx = new PVector(0, 0);
		rigid = r;
	}

	void show() {
		noStroke();
		if (rigid)fill(0);
		else fill(150 - pressure * 2000, 255, 200);
		ellipse(pos.x, pos.y, kernelRadius / 2, kernelRadius / 2);
	}

	void getNeighbors(grid g) {
		subGrid sG = g.getSubGrid(pos.x, pos.y);
		if (sG != null)neighbors = sG.particles;
	}
}

class grid {
	subGrid [] field;
	int nW, nH;
	float size;

	grid(float w, float h, float s) {
		size = s;
		nW = ceil(w / s);
		nH = ceil(h / s);
		field = new subGrid[nW * nH];
		for (int i = 0; i < field.length; i++) {
			field[i] = new subGrid();
		}
	}

	subGrid getSubGrid(float x, float y) {
		int i = floor(x / size);
		int j = floor(y / size);
		if ((i >= 0) && (i < nW) && (j >= 0) && (j < nH)) {
			int index = i * nH + j;
			return field[index];
		} else return null;
	}

	void addParticle(particle p) {
		for (int i = -1; i < 2; i++) {
			for (int j = -1; j < 2; j++) {
				float x = p.pos.x + size * float(i);
				float y = p.pos.y + size * float(j);
				subGrid sG = getSubGrid(x, y);
				if (sG != null)sG.particles.add(p);
			}
		}
	}

	void clearGrid() {
		for (subGrid sG : field)sG.particles.clear();
	}
}

class subGrid {
	ArrayList<particle> particles;

	subGrid() {
		particles = new ArrayList<particle>();
	}

	float getAvgPressure() {
		if (particles.size() == 0)return 0;
		float avgPressure = 0;
		for (particle p : particles)avgPressure += p.pressure;
		avgPressure /= particles.size();
		return avgPressure;
	}
}
