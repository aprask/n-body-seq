#include <iostream>
#include <cmath>
#include <string>
#include <random>
#include <chrono>

using namespace std;
using namespace std::chrono;

#define ARGS 5 // num of particles, time step size, num of iterations, how often to dump state
const double G = 6.674*pow(10,-11); // G for grav force
const double SOFTENING_FACTOR = 0.000000001; // this is just an arbitary float (I wasn't sure what "softening factor" actuall was...)

double calculateForce(double mass1, double coordinate1, double mass2, double coordinate2);
double calculateAcceleration(double force, double mass);
double calculatePosition(double coordinate, double velocity, double delta_t);
double calculateVelocity(double velocity, double acceleration, double delta_t);
double calculateDistance(double coordinate1, double coordinate2);

struct acceleration {
    double ax;
    double ay;
    double az;
};

struct position {
    double x;
    double y;
    double z;
};

struct velocity {
    double vx;
    double vy;
    double vz;
};

struct force {
    double fx;
    double fy;
    double fz;
};

struct particleNode {
    const double mass = rand() % 10000000;
    struct position position;
    struct velocity velocity;
    struct force force;
    struct acceleration acceleration;
};

int main (int argc, char* argv[]) {
    size_t size = 2;
    struct particleNode* particleField = (struct particleNode*)malloc(sizeof(struct particleNode) * size);
    if (!particleField) {
        cerr << "Cannot dynamically allocate mem for particle field" << endl;
        return -1;
    }
    for (int i = 0; i < size; ++i) {
        (particleField+i)->velocity.vx = rand() % 1000;
        (particleField+i)->velocity.vy = rand() % 1000;
        (particleField+i)->velocity.vz = rand() % 1000;
        (particleField+i)->force.fx = rand() % 1000;
        (particleField+i)->force.fy = rand() % 1000;
        (particleField+i)->force.fz = rand() % 1000;
        (particleField+i)->position.x = rand() % 1000;
        (particleField+i)->position.y = rand() % 1000;
        (particleField+i)->position.z = rand() % 1000;
    }
    for (int i = 0; i < size; ++i) {
        auto t_initial = steady_clock::now();
        for (int j = 0; j < size; ++j) {
            if (i == j) continue;
            auto delta_t = duration_cast<seconds>(steady_clock::now() - t_initial).count();

            (particleField+i)->force.fx = calculateForce((particleField+i)->mass, (particleField+i)->force.fx, (particleField+i)->mass, (particleField+j)->force.fx);
            (particleField+i)->force.fy = calculateForce((particleField+i)->mass, (particleField+i)->force.fy, (particleField+i)->mass, (particleField+j)->force.fy);
            (particleField+i)->force.fz = calculateForce((particleField+i)->mass, (particleField+i)->force.fz, (particleField+i)->mass, (particleField+j)->force.fz);

            (particleField+i)->acceleration.ax = calculateAcceleration((particleField+i)->force.fx, (particleField+i)->mass);
            (particleField+i)->acceleration.ay = calculateAcceleration((particleField+i)->force.fy, (particleField+i)->mass);
            (particleField+i)->acceleration.az = calculateAcceleration((particleField+i)->force.fz, (particleField+i)->mass);
            
            (particleField+i)->velocity.vx = calculateVelocity((particleField+i)->velocity.vx, (particleField+i)->acceleration.ax, delta_t);
            (particleField+i)->velocity.vy = calculateVelocity((particleField+i)->velocity.vy, (particleField+i)->acceleration.ay, delta_t);
            (particleField+i)->velocity.vz = calculateVelocity((particleField+i)->velocity.vz, (particleField+i)->acceleration.az, delta_t);

            (particleField+i)->position.x = calculatePosition((particleField+i)->position.x, (particleField+i)->velocity.vx, delta_t);
            (particleField+i)->position.y = calculatePosition((particleField+i)->position.y, (particleField+i)->velocity.vy, delta_t);
            (particleField+i)->position.z = calculatePosition((particleField+i)->position.z, (particleField+i)->velocity.vz, delta_t);
        }
    }
    return 0;
}

double calculateDistance(double coordinate1, double coordinate2) {
    return coordinate2-coordinate1;
}

double calculateForce(double mass1, double coordinate1, double mass2, double coordinate2) {
    double totalMass = mass1*mass2;
    double distanceDifference = calculateDistance(coordinate1, coordinate2);
    double r = distanceDifference/abs(distanceDifference);
    double rSquared = pow(r,2);
    return r*G*(totalMass/(rSquared+SOFTENING_FACTOR));
}

double calculatePosition(double coordinate, double velocity, double delta_t) {
    return coordinate + velocity*delta_t;
}

double calculateVelocity(double velocity, double acceleration, double delta_t) {
    return velocity + acceleration*delta_t;
}

double calculateAcceleration(double force, double mass) {
    return force/mass;
}