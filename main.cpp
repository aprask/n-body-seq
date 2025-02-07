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
};

int main (int argc, char* argv[]) {
    struct particleNode* testParticle1 = (struct particleNode*)malloc(sizeof(struct particleNode));
    testParticle1->velocity = 500.25;
    struct particleNode* testParticle2 = (struct particleNode*)malloc(sizeof(struct particleNode));
    testParticle2->velocity = 8900.25;
    struct particleNode* testParticle3 = (struct particleNode*)malloc(sizeof(struct particleNode));
    testParticle3->velocity = 23323.25;

    testParticle1->nextParticleInSeq = testParticle2;
    testParticle2->nextParticleInSeq = testParticle3;

    cout << "particle1's next particle's velocity: " << testParticle1->nextParticleInSeq->velocity << endl;;
    cout << "particle1's next next particle's velocity: " << testParticle1->nextParticleInSeq->nextParticleInSeq->velocity << endl;

    free(testParticle1);
    free(testParticle2);
    free(testParticle3);
    return 0;
}