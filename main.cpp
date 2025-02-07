#include <iostream>
#include <cmath>

using namespace std;

#define G  6.674(10e-11)
#define ARGS 5 // num of particles, time step size, num of iterations, how often to dump state

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