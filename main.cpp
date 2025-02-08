#include <iostream>
#include <cmath>
#include <ctime>

using namespace std;

#define ARGS 5 // num of particles, time step size, num of iterations, how often to dump state
const double G = 6.674*pow(10,-11); // G for grav force
const double SOFTENING_FACTOR = 0.000000001; // this is just an arbitary float (I wasn't sure what "softening factor" actuall was...)

double calculateForce(struct particleNode* nodeA, struct particleNode* nodeB);
double calculateDistance(double x1, double y1, double z1, double x2, double y2, double z2);
double caclulateDisplacement(double x1, double y1, double z1, double x2, double y2, double z2);
double calculateAverageVelocity(struct particleNode* p);
double calculateInstVelocity(struct particleNode* p);
double calculateChangeInPosition(struct particleNode* p);
struct position {
    double x;
    double y;
    double z;
};

struct particleNode {
    double mass;
    struct position oldPosition;
    struct position position;
    double oldVelocity;
    double velocity;
    double force;
    double acceleration;
    time_t begTime;
    time_t endTime;
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

    p1->oldVelocity = -5;
    p1->acceleration = 10;
    p1->begTime = 1;
    p1->endTime = 20;

    double instVelocity = calculateInstVelocity(p1);
    cout << "Inst Velocity: " << instVelocity << endl;
    return 0;
}

double calculateDistance(double x1, double y1, double z1, double x2, double y2, double z2) {
    return sqrt(pow((x2-x1), 2) + pow((y2-y1), 2) + pow((z2-z1), 2));
}

double caclulateDisplacement(double x1, double y1, double z1, double x2, double y2, double z2) {
    return (x2-x1)+(y2-y1)+(z2-z1);
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
    double seconds = difftime(p->begTime, p->endTime); // reference: https://en.cppreference.com/w/c/chrono/difftime
    double displacement = caclulateDisplacement(
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
    return p->oldVelocity + (p->acceleration*(p->endTime-p->begTime));
}

double calculateChangeInPosition(struct particleNode* p) {
    double seconds = difftime(p->begTime, p->endTime);
    double velocityTimesTime = p->velocity*seconds;
    double changeInOldPosition = caclulateDisplacement(
        p->oldPosition.x,
        p->oldPosition.y,
        p->oldPosition.z,
        p->position.x,
        p->position.y,
        p->position.z
    );
    return (changeInOldPosition + velocityTimesTime);
}