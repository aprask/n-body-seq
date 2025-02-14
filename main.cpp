#include <iostream>
#include <cmath>
#include <string>
#include <random>
#include <fstream>
#include <ctype.h>
#include <cstring>
#include <filesystem>
#include <sstream>
#include <chrono>
#include <ctime>

#define ARGS 5

class Particle {
    private:
        const long double G = 6.674*pow(10,-11);
        const long double SOFTENING_FACTOR = 0.000001;
        struct acceleration {
            double ax;
            double ay;
            double az;
        };
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

    public:
        struct force force;
        struct velocity velocity;
        struct acceleration acceleration;
        struct position position;
        long double mass;

        Particle() {
            this->force.fx = 0;
            this->force.fy = 0;
            this->force.fz = 0;
            this->velocity.vx = rand() % 10000000;
            this->velocity.vy = rand() % 10000000;
            this->velocity.vz = rand() % 10000000;
            this->acceleration.ax = rand() % 10000000;
            this->acceleration.ay = rand() % 10000000;
            this->acceleration.az = rand() % 10000000;
            this->position.x = rand() % 10000000;
            this->position.y = rand() % 10000000;
            this->position.z = rand() % 10000000;
            this->mass = rand() % 1000000000000;
        }

        Particle(
            long double forceFx, 
            long double forceFy, 
            long double forceFz,
            long double velocityVx,
            long double velocityVy,
            long double velocityVz,
            long double accelerationAx,
            long double accelerationAy,
            long double accelerationAz,
            long double positionX,
            long double positionY,
            long double positionZ,
            long double mass
        ) {
            this->force.fx = forceFx;
            this->force.fy = forceFy;
            this->force.fz = forceFz;
            this->velocity.vx = velocityVx;
            this->velocity.vy = velocityVy;
            this->velocity.vz = velocityVz;
            this->acceleration.ax = accelerationAx;
            this->acceleration.ay = accelerationAy;
            this->acceleration.az = accelerationAz;
            this->position.x = positionX;
            this->position.y = positionY;
            this->position.z = positionZ;
            this->mass = mass;
        }

        ~Particle() {}

        long calculateForceDirectionFx(const Particle& otherParticle) {
            long double forceMagnitude = this->calculateForceMagnitude(otherParticle);
            long double distanceVector = this->calculateEuclideanDistance(otherParticle);
            long double linearDifference = this->position.x - otherParticle.position.x;
            return forceMagnitude*(linearDifference/distanceVector);
        }

        long calculateForceDirectionFy(const Particle& otherParticle) {
            long double forceMagnitude = this->calculateForceMagnitude(otherParticle);
            long double distanceVector = this->calculateEuclideanDistance(otherParticle);
            long double linearDifference = this->position.y - otherParticle.position.y;
            return forceMagnitude*(linearDifference/distanceVector);
        }

        long calculateForceDirectionFz(const Particle& otherParticle) {
            long double forceMagnitude = this->calculateForceMagnitude(otherParticle);
            long double distanceVector = this->calculateEuclideanDistance(otherParticle);
            long double linearDifference = this->position.z - otherParticle.position.z;
            return forceMagnitude*(linearDifference/distanceVector);
        }

        long double calculateForceMagnitude(const Particle& otherParticle) {
            long double totalMassScalar = this->mass * otherParticle.mass;
            long double distanceVector = this->calculateEuclideanDistance(otherParticle);
            return G*(totalMassScalar/(pow(distanceVector,2)+SOFTENING_FACTOR));
        }

        long double calculateEuclideanDistance(const Particle& otherParticle) {
            return sqrt((pow(this->position.x-otherParticle.position.x,2)
            + pow(this->position.y-otherParticle.position.y,2),
            + pow(this->position.z-otherParticle.position.z,2)));
        }

        long double calculatePositionX(const long double delta_t) {
            return this->position.x * (this->velocity.vx*delta_t);
        }

        long double calculatePositionY(const long double delta_t) {
            return this->position.y * (this->velocity.vy*delta_t);
        }

        long double calculatePositionZ(const long double delta_t) {
            return this->position.z * (this->velocity.vz*delta_t);
        }

        long double calculateAccelerationX() {
            return this->force.fx/this->mass;
        }

        long double calculateAccelerationY() {
            return this->force.fy/this->mass;
        }

        long double calculateAccelerationZ() {
            return this->force.fz/this->mass;
        }

        long double calculateVelocityX(const long double delta_t) {
            return this->velocity.vx + this->acceleration.ax*delta_t;
        }

        long double calculateVelocityY(const long double delta_t) {
            return this->velocity.vy + this->acceleration.ay*delta_t;
        }

        long double calculateVelocityZ(const long double delta_t) {
            return this->velocity.vz + this->acceleration.az*delta_t;
        }

        long double getForceFX() const {
            return this->force.fx;
        }

        void setForceFx(long double forceFx) {
            this->force.fx = forceFx;
        }

        long double getForceFy() const {
            return this->force.fy;
        }

        void setForceFy(long double forceFy) {
            this->force.fy = forceFy;
        }

        long double getForceFz() const {
            return this->force.fz;
        }

        void setForceFz(long double forceFz) {
            this->force.fz = forceFz;
        }

        long double getVelocityVx() const {
            return this->velocity.vx;
        }

        void setVelocityVx(long double velocityVx) {
            this->velocity.vx = velocityVx;
        }

        long double getVelocityVy() const {
            return this->velocity.vy;
        }

        void setVelocityVy(long double velocityVy) {
            this->velocity.vy = velocityVy;
        }

        long double getVelocityVz() const {
            return this->velocity.vz;
        }

        void setVelocityVz(long double velocityVz) {
            this->velocity.vz = velocityVz;
        }

        long double getAccelerationAx() const {
            return this->acceleration.ax;
        }

        void setAccelerationAx(long double accelerationAx) {
            this->acceleration.ax = accelerationAx;
        }

        long double getAccelerationAy() const {
            return this->acceleration.ay;
        }

        void setAccelerationAy(long double accelerationAy) {
            this->acceleration.ay = accelerationAy;
        }
        
        long double getAccelerationAz() const {
            return this->acceleration.az;
        }

        void setAccelerationAz(long double accelerationAz) {
            this->acceleration.az = accelerationAz;
        }

        long double getPositionX() const {
            return this->position.x;
        }

        void setPositionX(long double positionX) {
            this->position.x = positionX;
        }

        long double getPositionY() const {
            return this->position.y;
        }

        void setPositionY(long double positionY) {
            this->position.y = positionY;
        }

        long double getPositionZ() const {
            return this->position.z;
        }

        void setPositionZ(long double positionZ) {
            this->position.z = positionZ;
        }

        long double getMass() const {
            return this->mass;
        }

        void setMass(long double mass) {
            this->mass = mass;
        }
};


void parseFile(const std::string& fileName, size_t& N, std::vector<Particle*>* particleField);
void particleInitFromFile(const int& tokenIdx, std::vector<Particle*>* particleField, const int& classIdx, std::string token);


int main(int argc, char* argv[]) {
    bool readFromFile = false;
    size_t N = 0;
    std::vector<Particle*> particleField;
    if (argc != ARGS) {
        std::cerr << "Invalid number of arguments <" << argc << ">" << ". Expected <" << ARGS << ">" << std::endl;
        return 1;
    } else {
        for (int i = 1; i < ARGS; ++i) {
            for (int j = 0; j < strlen(argv[i]); ++j) {
                if (!isdigit(argv[i][j])) {
                    if (i == 1) {
                        std::string fileName = argv[i];
                        fileName = fileName + ".tsv";
                        std::ifstream sourceFile(fileName);
                        if (!sourceFile.is_open()) {
                            std::cerr << "Failed to open file: " << argv[i] << std::endl;
                            return 1;
                        } else {
                            readFromFile = true;
                            std::string line;
                            while(getline(sourceFile, line)) {
                                std::stringstream s(line);
                                std::string token;
                                int col = 0;
                                while(getline(s, token, '\t')) {
                                    col++;
                                    if(!(col % 10)) N++;
                                }
                            }
                            sourceFile.close();
                            std::cout << "I am here" << std::endl;                    
                            parseFile(fileName, N, &particleField);
                            break;
                        }
                    } else {
                        std::cerr << "Invalid numerical argument type passed: " << argv[i] << std::endl;
                        return 1;
                    }
                }
            }
        }
    }
    return 0;
}

void particleInitFromFile(const int& tokenIdx, std::vector<Particle*>* particleField, const int& classIdx, std::string token) {
    switch (tokenIdx) {
        case 0: // we skip the first col (which is particles)
            break;
        case 1:
            (*particleField)[classIdx]->setMass(stold(token));
            std::cout << "Mass: " << (*particleField)[classIdx]->getMass() << std::endl;
            break; // mass
        case 2:
            (*particleField)[classIdx]->setPositionX(stold(token));
            std::cout << "X: " << (*particleField)[classIdx]->getPositionX() << std::endl;
            break; // x
        case 3:
            (*particleField)[classIdx]->setPositionY(stold(token));
            std::cout << "Y: " << (*particleField)[classIdx]->getPositionY() << std::endl;
            break; // y
        case 4:
            (*particleField)[classIdx]->setPositionZ(stold(token));
            std::cout << "Z: " << (*particleField)[classIdx]->getPositionZ() << std::endl;
            break; // z
        case 5:
            (*particleField)[classIdx]->setVelocityVx(stold(token));
            std::cout << "Vx: " << (*particleField)[classIdx]->getVelocityVx() << std::endl;
            break; // vx
        case 6:
            (*particleField)[classIdx]->setVelocityVy(stold(token));
            std::cout << "Vy: " << (*particleField)[classIdx]->getVelocityVy() << std::endl;
            break; // vy
        case 7:
            (*particleField)[classIdx]->setVelocityVz(stold(token));
            std::cout << "Vy: " << (*particleField)[classIdx]->getVelocityVz() << std::endl;
            break; // vz
        case 8:
            (*particleField)[classIdx]->setForceFx(stold(token));
            std::cout << "Fx: " << (*particleField)[classIdx]->getForceFX() << std::endl;
            break; // fx
        case 9:
            (*particleField)[classIdx]->setForceFy(stold(token));
            std::cout << "Fy: " << (*particleField)[classIdx]->getForceFy() << std::endl;
            break; // fy
        case 10:
            (*particleField)[classIdx]->setForceFz(stold(token));
            std::cout << "Fz: " << (*particleField)[classIdx]->getForceFz() << std::endl;
            break; // fz
        default:
            break;
    }                                    
}

void parseFile(const std::string& fileName, size_t& N, std::vector<Particle*>* particleField) {
    std::cout << "Size of N: " << N << std::endl;
    std::cout << "I am here" << std::endl;
    for (int i = 0; i < N; ++i) (*particleField).push_back(new Particle());
    std::string line;
    std::ifstream sourceFile(fileName);
    if (!sourceFile.is_open()) {
        std::cerr << "Failed to open file" << std::endl;
        exit(1);
    }
    std::cout << "I am here" << std::endl;
    while (getline(sourceFile, line)) {
        std::stringstream s(line);
        std::string token;
        int tokenIdx = 0;
        int i = 0;
        while (getline(s, token, '\t')) {
            if (tokenIdx > 10) {
                tokenIdx = 1;
                i++;
            }
            particleInitFromFile(tokenIdx, particleField, i, token);
            tokenIdx++;
        }
        sourceFile.close();
    }
}