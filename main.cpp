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
double calculateAcceleration(struct particleNode* p, double force);
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
    double mass;
    struct position position;
    struct velocity velocity;
    struct force force;
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
        (particleField+i)->mass = rand() % 1000;
    }
    for (int i = 0; i < size; ++i) {
        auto t_initial = steady_clock::now();
        for (int j = 0; j < size; ++j) {
            if (i == j) continue;
            auto delta_t = duration_cast<seconds>(steady_clock::now() - t_initial).count(); // change in time (t2-t1)
            
            
            // double changeInX = calculateChangeInPosition((particleField+i)->position.x, (particleField+i), delta_t); // new x
            // double changeInY = calculateChangeInPosition((particleField+i)->position.y, (particleField+i), delta_t); // new y
            // double changeInZ = calculateChangeInPosition((particleField+i)->position.z, (particleField+i), delta_t); // new z
            // double forceBetweenTwoParticles = calculateForce((particleField+i), (particleField+j));
            // double acceleration = calculateAcceleration((particleField+i), forceBetweenTwoParticles);
            // double newVelocity = calculateInstVelocity((particleField+i), delta_t, acceleration);
            // (particleField+i)->force += forceBetweenTwoParticles; // accumulating force (total force acting on particle i)
            // (particleField+i)->velocity = newVelocity;
            // (particleField+i)->position.x = changeInX;
            // (particleField+i)->position.y = changeInY;
            // (particleField+i)->position.z = changeInZ;
        }
    }
    return 0;
}

double calculateDistance(double coordinate1, double coordinate2) {\
    return coordinate2-coordinate1;
}

double calculateForce(double mass1, double coordinate1, double mass2, double coordinate2) {
    double totalMass = mass1 * mass2;
    cout << "Total Mass: " << totalMass << endl;
    double distance = calculateDistance(coordinate1, coordinate2);
    cout << "Distance: " << distance << endl;
    double distancedSquared = pow(distance, 2);
    cout << "Distance Squared: " << distancedSquared << endl;
    double directionOfForce = abs(distance)/distance;
    cout << "Direction of Force: " << directionOfForce << endl;
    double gravitationalForce = G * (totalMass/((distancedSquared + SOFTENING_FACTOR)));
    return directionOfForce*gravitationalForce;
}

// double calculateInstVelocity(struct particleNode* p, double delta_t, double acceleration) {
//     return p->velocity + (acceleration*(delta_t));
// }

// double calculateChangeInPosition(double coordinate, struct particleNode* p, double delta_t) {
//     double velocityTimesTime = p->velocity*delta_t;
//     return (coordinate + velocityTimesTime);
// }

double calculateAcceleration(struct particleNode* p, double force) {
    return force/p->mass;
}