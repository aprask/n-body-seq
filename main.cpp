#include <iostream>
#include <cmath>
#include <ctime>

using namespace std;

#define ARGS 5 // num of particles, time step size, num of iterations, how often to dump state
const double G = 6.674*pow(10,-11); // G for grav force
const double SOFTENING_FACTOR = 0.000000001;

double calculateForce(struct particleNode* nodeA, struct particleNode* nodeB);
double calculateDistance(double x1, double y1, double z1, double x2, double y2, double z2);
double caclulateDisplacement(double x1, double y1, double z1, double x2, double y2, double z2);
double calculateAverageVelocity(struct particleNode* p);
struct position {
    double x;
    double y;
    double z;
};

struct particleNode {
    double mass;
    struct position started;
    struct position position;
    double velocity;
    double force;
    time_t begTime;
    struct tm elaspedTime;
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

    p1->started.x = 1;
    p1->started.y = 2;
    p1->started.z = 9;

    time_t testTime = time(0);

    struct tm elapsedTime = *localtime(&testTime);

    elapsedTime.tm_hour = 1;
    elapsedTime.tm_sec = 49;

    p1->begTime = testTime;
    p1->elaspedTime = elapsedTime;

    printf("Avg Velocity is %f\n", calculateAverageVelocity(p1));
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
    double seconds = difftime(p->begTime, mktime(&p->elaspedTime)); // reference: https://en.cppreference.com/w/c/chrono/difftime
    double displacement = caclulateDisplacement(
        p->started.x,
        p->started.y,
        p->started.z,
        p->position.x, 
        p->position.y, 
        p->position.z
    );
    return displacement/seconds;
}