#include <iostream>
#include <cmath>

using namespace std;

#define ARGS 5 // num of particles, time step size, num of iterations, how often to dump state
const double G = 6.674*pow(10,-11); // G for grav force

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
    struct particleNode* p1 = (struct particleNode*)malloc(sizeof(struct particleNode));
    struct particleNode* p2 = (struct particleNode*)malloc(sizeof(struct particleNode));
    p1->mass = 5453535;
    p1->position.x = 33;
    p1->position.y = 33;
    p1->position.z = 54;
    
    p2->mass = 454524;
    p2->position.x = 34;
    p2->position.y = 2;
    p2->position.z = 3;

    printf("Force is %f\n", calculateForce(p1, p2));
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
    double gravitationalForce = G * (totalMass/(distancedSquared));
    return directionOfForce*gravitationalForce;
}