#include <iostream>
#include <cmath>
#include <ctime>
#include <string>
#include <random>

using namespace std;

#define ARGS 5 // num of particles, time step size, num of iterations, how often to dump state
const double G = 6.674*pow(10,-11); // G for grav force
const double SOFTENING_FACTOR = 0.000000001; // this is just an arbitary float (I wasn't sure what "softening factor" actuall was...)

double calculateForce(struct particleNode* nodeA, struct particleNode* nodeB);
double calculateDistance(double x1, double y1, double z1, double x2, double y2, double z2);
double calculateDistanceMagnitude(double distance);
double calculateAverageVelocity(struct particleNode* p);
double calculateInstVelocity(struct particleNode* p);
double calculateChangeInPosition(struct particleNode* p);
double calculateAcceleration(struct particleNode* p, double force);
struct position {
    double x;
    double y;
    double z;
};

struct particleNode {
    double mass;
    struct position oldPosition;
    struct position position;
    double velocity;
    double force;
    double acceleration;
    time_t begTime;
    time_t endTime;
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

double calculateAverageVelocity(struct particleNode* p) {
    double seconds = difftime(p->endTime, p->begTime); // reference: https://en.cppreference.com/w/c/chrono/difftime
    double displacement = calculateDistance(
        p->oldPosition.x,
        p->oldPosition.y,
        p->oldPosition.z,
        p->position.x, 
        p->position.y,
        p->position.z
    );
    return displacement/seconds;
}

double calculateInstVelocity(struct particleNode* p) {
    return p->velocity + (p->acceleration*(p->endTime-p->begTime));
}

double calculateChangeInPosition(struct particleNode* p) {
    double seconds = difftime(p->endTime, p->begTime);
    double velocityTimesTime = p->velocity*seconds;
    double changeInOldPosition = calculateDistance(
        p->oldPosition.x,
        p->oldPosition.y,
        p->oldPosition.z,
        p->position.x,
        p->position.y,
        p->position.z
    );
    return (changeInOldPosition + velocityTimesTime);
}

double calculateAcceleration(struct particleNode* p, double force) {
    return force/p->mass;
}