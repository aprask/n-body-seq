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

double calculateForce(struct particleNode* nodeA, struct particleNode* nodeB);
double calculateDistance(double x1, double y1, double z1, double x2, double y2, double z2);
double calculateDistanceMagnitude(double distance);
double calculateAverageVelocity(struct particleNode* p, double delta_t);
double calculateInstVelocity(struct particleNode* p, double delta_t);
double calculateChangeInPosition(double coordinate, struct particleNode* p, double delta_t);
double calculateAcceleration(struct particleNode* p, double force);
struct position {
    double x;
    double y;
    double z;
};

struct particleNode {
    double mass;
    struct position position;
    double velocity;
    double force;
    double acceleration;
};

int main (int argc, char* argv[]) {
    size_t size = 2;
    struct particleNode* particleField = (struct particleNode*)malloc(sizeof(struct particleNode) * size);
    if (!particleField) {
        cerr << "Cannot dynamically allocate mem for particle field" << endl;
        return -1;
    }
    for (int i = 0; i < size; ++i) {
        (particleField+i)->velocity = rand() % 1000;
        (particleField+i)->acceleration = rand() % 1000;
        (particleField+i)->force = rand() % 1000;
        (particleField+i)->position.x = rand() % 1000;
        (particleField+i)->position.y = rand() % 1000;
        (particleField+i)->position.z = rand() % 1000;
        (particleField+i)->mass = rand() % 1000;
    }
    for (int i = 0; i < size; ++i) {
        for (int j = 0; i < size; ++j) (particleField+j)->force = 0; // resetting the force
        auto t_initial = steady_clock::now();
        for (int j = 0; j < size; ++j) {
            if (i == j) continue;
            auto delta_t = duration_cast<seconds>(steady_clock::now() - t_initial).count(); // change in time (t2-t1)
            double changeInX = calculateChangeInPosition((particleField+i)->position.x, (particleField+i), delta_t); // new x
            double changeInY = calculateChangeInPosition((particleField+i)->position.y, (particleField+i), delta_t); // new y
            double changeInZ = calculateChangeInPosition((particleField+i)->position.z, (particleField+i), delta_t); // new z

        }
    }
    return 0;
}

double calculateDistance(double x1, double y1, double z1, double x2, double y2, double z2) {
    return sqrt(pow((x2-x1), 2) + pow((y2-y1), 2) + pow((z2-z1), 2));
}

double calculateForce(struct particleNode* nodeA, struct particleNode* nodeB) {
    double totalMass = nodeA->mass * nodeB->mass; // distance between particles in space
    cout << "Total Mass: " << totalMass << endl;
    double distance = calculateDistance(
        nodeA->position.x,
        nodeA->position.y,
        nodeA->position.z,
        nodeB->position.x,
        nodeB->position.y,
        nodeB->position.z
    );
    cout << "Distance: " << distance << endl;
    double distancedSquared = pow(distance, 2);
    cout << "Distance Squared: " << distancedSquared << endl;
    double directionOfForce = abs(distance)/distance;
    cout << "Direction of Force: " << directionOfForce << endl;
    double gravitationalForce = G * (totalMass/((distancedSquared + SOFTENING_FACTOR)));
    return directionOfForce*gravitationalForce;
}

double calculateInstVelocity(struct particleNode* p, double delta_t) {
    return p->velocity + (p->acceleration*(delta_t));
}

double calculateChangeInPosition(double coordinate, struct particleNode* p, double delta_t) {
    double velocityTimesTime = p->velocity*delta_t;
    return (coordinate + velocityTimesTime);
}

double calculateAcceleration(struct particleNode* p, double force) {
    return force/p->mass;
}