#include <iostream>
#include <cmath>

using namespace std;

#define G  6.674*pow(10,-11);
#define ARGS 5 // num of particles, time step size, num of iterations, how often to dump state

double calculateForce(struct particleNode* nodeA, struct particleNode* nodeB);
struct particleNode {
    double mass;
    double position;
    double velocity;
    double force;
    struct particleNode* nextParticleInSeq;
    struct particleNode* previousParticleInSeq;
};

int main (int argc, char* argv[]) {

    return 0;
}

double calculateForce(struct particleNode* nodeA, struct particleNode* nodeB) {
    double totalMass = nodeA->mass * nodeB->mass; // distance between particles in space
    double distanceBetweenParticles = nodeB->position - nodeA->position; // r
    distanceBetweenParticles = pow(distanceBetweenParticles, 2); // r squared
    return G(totalMass/distanceBetweenParticles);
}