#include <iostream>
#include <cmath>

using namespace std;

#define G  6.674*pow(10,-11);
#define ARGS 5 // num of particles, time step size, num of iterations, how often to dump state

double calculateForce(struct particleNode* nodeA, struct particleNode* nodeB);
double calculateDistance(double x1, double y1, double z1, double x2, double y2, double z2);
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
    struct particleNode* nextParticleInSeq;
    struct particleNode* previousParticleInSeq;
};

int main (int argc, char* argv[]) {

    return 0;
}

double calculateDistance(double x1, double y1, double z1, double x2, double y2, double z2) {
    return sqrt(pow((x2-x1), 2) + pow((y2-y1), 2) + pow((z2-z1), 2));
}

double calculateForce(struct particleNode* nodeA, struct particleNode* nodeB) {
    double totalMass = nodeA->mass * nodeB->mass; // distance between particles in space
    double distance = calculateDistance(
        nodeA->position.x,
        nodeA->position.y,
        nodeA->position.z,
        nodeB->position.x,
        nodeB->position.y,
        nodeB->position.z
    );
    double distancedSquared = pow(distance, 2);
    double directionOfForce = distance/abs(distance);
    return (directionOfForce)*G(totalMass/distancedSquared);
}