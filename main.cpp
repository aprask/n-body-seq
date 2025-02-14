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
        const double G = 6.674*pow(10,-11);
        const double SOFTENING_FACTOR = 0.000001;
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
        double mass;

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
            double forceFx, 
            double forceFy, 
            double forceFz,
            double velocityVx,
            double velocityVy,
            double velocityVz,
            double accelerationAx,
            double accelerationAy,
            double accelerationAz,
            double positionX,
            double positionY,
            double positionZ,
            double mass
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

        void printProperties() {
            std::cout << "Mass: " << this->mass << std::endl;

            std::cout << "X: " << this->position.x << std::endl;
            std::cout << "Y: " << this->position.y << std::endl;
            std::cout << "Z: " << this->position.z << std::endl;

            std::cout << "Vx: " << this->velocity.vx << std::endl;
            std::cout << "Vy: " << this->velocity.vy << std::endl;
            std::cout << "Vz: " << this->velocity.vz << std::endl;

            std::cout << "Ax: " << this->acceleration.ax << std::endl;
            std::cout << "Ay: " << this->acceleration.ay << std::endl;
            std::cout << "Az: " << this->acceleration.az << std::endl;

            std::cout << "Fx: " << this->force.fx << std::endl;
            std::cout << "Fy: " << this->force.fy << std::endl;
            std::cout << "Fz: " << this->force.fz << std::endl;
        }

        void clearForce() {
            this->force.fx = 0;
            this->force.fy = 0;
            this->force.fz = 0;
        }

        void updatePosition(const size_t& DELTA_T) {
            this->acceleration.ax = this->calculateAccelerationX();
            this->acceleration.ay = this->calculateAccelerationY();
            this->acceleration.az = this->calculateAccelerationZ();

            this->velocity.vx = this->calculateVelocityX(DELTA_T);
            this->velocity.vy = this->calculateVelocityY(DELTA_T);
            this->velocity.vz = this->calculateVelocityZ(DELTA_T);

            this->position.x = this->calculatePositionX(DELTA_T);
            this->position.y = this->calculatePositionY(DELTA_T);
            this->position.z = this->calculatePositionZ(DELTA_T);
        }

        double calculateForceDirectionFx(const Particle* otherParticle) {
            double forceMagnitude = this->calculateForceMagnitude(otherParticle);
            double distanceVector = this->calculateEuclideanDistance(otherParticle);
            double linearDifference = this->position.x - otherParticle->position.x;
            return forceMagnitude*(linearDifference/distanceVector);
        }

        double calculateForceDirectionFy(const Particle* otherParticle) {
            double forceMagnitude = this->calculateForceMagnitude(otherParticle);
            double distanceVector = this->calculateEuclideanDistance(otherParticle);
            double linearDifference = this->position.y - otherParticle->position.y;
            return forceMagnitude*(linearDifference/distanceVector);
        }

        double calculateForceDirectionFz(const Particle* otherParticle) {
            double forceMagnitude = this->calculateForceMagnitude(otherParticle);
            double distanceVector = this->calculateEuclideanDistance(otherParticle);
            double linearDifference = this->position.z - otherParticle->position.z;
            return forceMagnitude*(linearDifference/distanceVector);
        }

        double calculateForceMagnitude(const Particle* otherParticle) {
            double totalMassScalar = this->mass * otherParticle->mass;
            double distanceVector = this->calculateEuclideanDistance(otherParticle);
            return G*(totalMassScalar/(pow(distanceVector,2)+SOFTENING_FACTOR));
        }

        double calculateEuclideanDistance(const Particle* otherParticle) {
            return sqrt((pow(this->position.x-otherParticle->position.x,2)
            + pow(this->position.y-otherParticle->position.y,2) + pow(this->position.z-otherParticle->position.z,2)));
        }

        double calculatePositionX(const double delta_t) {
            return this->position.x + (this->velocity.vx*delta_t);
        }

        double calculatePositionY(const double delta_t) {
            return this->position.y + (this->velocity.vy*delta_t);
        }

        double calculatePositionZ(const double delta_t) {
            return this->position.z + (this->velocity.vz*delta_t);
        }

        double calculateAccelerationX() {
            return this->force.fx/this->mass;
        }

        double calculateAccelerationY() {
            return this->force.fy/this->mass;
        }

        double calculateAccelerationZ() {
            return this->force.fz/this->mass;
        }

        double calculateVelocityX(const double delta_t) {
            return this->velocity.vx + this->acceleration.ax*delta_t;
        }

        double calculateVelocityY(const double delta_t) {
            return this->velocity.vy + this->acceleration.ay*delta_t;
        }

        double calculateVelocityZ(const double delta_t) {
            return this->velocity.vz + this->acceleration.az*delta_t;
        }

        double getForceFX() const {
            return this->force.fx;
        }

        void setForceFx(const double& forceFx) {
            this->force.fx += forceFx;
        }

        void updateOtherFx(Particle* otherParticle) {
            double forceMagnitude = otherParticle->calculateForceMagnitude(this);
            double distanceVector = otherParticle->calculateEuclideanDistance(this);
            double linearDifference = otherParticle->position.x - this->position.x;
            double forceComp = forceMagnitude*(linearDifference/distanceVector);
            otherParticle->force.fx -= forceComp;
        }

        double getForceFy() const {
            return this->force.fy;
        }

        void setForceFy(const double& forceFy) {
            this->force.fy += forceFy;
        }

        void updateOtherFy(Particle* otherParticle) {
            double forceMagnitude = otherParticle->calculateForceMagnitude(this);
            double distanceVector = otherParticle->calculateEuclideanDistance(this);
            double linearDifference = otherParticle->position.y - this->position.y;
            double forceComp = forceMagnitude*(linearDifference/distanceVector);
            otherParticle->force.fy -= forceComp;
        }

        double getForceFz() const {
            return this->force.fz;
        }

        void setForceFz(const double& forceFz) {
            this->force.fz += forceFz;
        }

        void updateOtherFz(Particle* otherParticle) {
            double forceMagnitude = otherParticle->calculateForceMagnitude(this);
            double distanceVector = otherParticle->calculateEuclideanDistance(this);
            double linearDifference = otherParticle->position.z - this->position.z;
            double forceComp = forceMagnitude*(linearDifference/distanceVector);
            otherParticle->force.fz -= forceComp;
        }

        double getVelocityVx() const {
            return this->velocity.vx;
        }

        void setVelocityVx(const double& velocityVx) {
            this->velocity.vx = velocityVx;
        }

        double getVelocityVy() const {
            return this->velocity.vy;
        }

        void setVelocityVy(const double& velocityVy) {
            this->velocity.vy = velocityVy;
        }

        double getVelocityVz() const {
            return this->velocity.vz;
        }

        void setVelocityVz(const double& velocityVz) {
            this->velocity.vz = velocityVz;
        }

        double getAccelerationAx() const {
            return this->acceleration.ax;
        }

        void setAccelerationAx(const double& accelerationAx) {
            this->acceleration.ax = accelerationAx;
        }

        double getAccelerationAy() const {
            return this->acceleration.ay;
        }

        void setAccelerationAy(const double& accelerationAy) {
            this->acceleration.ay = accelerationAy;
        }
        
        double getAccelerationAz() const {
            return this->acceleration.az;
        }

        void setAccelerationAz(const double& accelerationAz) {
            this->acceleration.az = accelerationAz;
        }

        double getPositionX() const {
            return this->position.x;
        }

        void setPositionX(const double& positionX) {
            this->position.x = positionX;
        }

        double getPositionY() const {
            return this->position.y;
        }

        void setPositionY(const double& positionY) {
            this->position.y = positionY;
        }

        double getPositionZ() const {
            return this->position.z;
        }

        void setPositionZ(const double& positionZ) {
            this->position.z = positionZ;
        }

        double getMass() const {
            return this->mass;
        }

        void setMass(double mass) {
            this->mass = mass;
        }
};


void parseFile(const std::string& fileName, const size_t& N, std::vector<Particle*>* particleField);
void particleInitFromFile(const size_t& tokenIdx, std::vector<Particle*>* particleField, const size_t& classIdx, const std::string& token);
void nBodySim(const size_t& N, const size_t& TIME_STEPS, const size_t& DELTA_T, const size_t& DUMP_RATE, std::vector<Particle*>* particleField);
void randInit(const size_t& N, std::vector<Particle*>* particleField);
void dumpData(const size_t& N, std::vector<Particle*>* particleField, std::ofstream* file);

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
                            parseFile(fileName, N, &particleField);
                            break;
                        }
                    } else {
                        std::cerr << "Invalid numerical argument type passed: " << argv[i] << std::endl;
                        return 1;
                    }
                } else {
                    if (readFromFile) break;
                    N = std::stoi(argv[1]);
                    randInit(N, &particleField);
                }
            }
        }
    }
    const size_t DELTA_T = std::stoi(argv[2]);
    const size_t TIME_STEPS = std::stoi(argv[3]);
    const size_t DUMP_RATE = std::stoi(argv[4]);
    nBodySim(N, DELTA_T, TIME_STEPS, DUMP_RATE, &particleField);
    return 0;
}

void dumpData(const size_t& N, std::vector<Particle*>* particleField, std::ofstream* file) {
    if (particleField == nullptr) {
        std::cerr << "Particle field pointer is NULL" << std::endl;
        exit(1);
    }
    *(file) << N << "\t";
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            *(file) << (*particleField)[j]->getMass() << "\t";
            *(file) << (*particleField)[j]->getPositionX() << "\t";
            *(file) << (*particleField)[j]->getPositionY() << "\t";
            *(file) << (*particleField)[j]->getPositionZ() << "\t";
            *(file) << (*particleField)[j]->getVelocityVx() << "\t";
            *(file) << (*particleField)[j]->getVelocityVy() << "\t";
            *(file) << (*particleField)[j]->getVelocityVz() << "\t";
            *(file) << (*particleField)[j]->getForceFX() << "\t";
            *(file) << (*particleField)[j]->getForceFy() << "\t";
            *(file) << (*particleField)[j]->getForceFz() << "\t";
        }
    }
    *(file) << std::endl;
}

void randInit(const size_t& N, std::vector<Particle*>* particleField) {
    if (particleField == nullptr) {
        std::cerr << "Particle field pointer is NULL" << std::endl;
        exit(1);
    }
    for (int i = 0; i < N; ++i) particleField->push_back(new Particle());
}

void nBodySim(const size_t& N, const size_t& TIME_STEPS, const size_t& DELTA_T, const size_t& DUMP_RATE, std::vector<Particle*>* particleField) {
    if (particleField == nullptr) {
        std::cerr << "Particle field pointer is NULL" << std::endl;
        exit(1);
    }
    std::ofstream dataFile;
    char* buffer = new char[40]; // something arbitrary but big enough to hold a file
    if (!buffer) {
        std::cerr << "Failed to allocate memory for buffer" << std::endl;
        exit(1);
    }
    // https://cplusplus.com/reference/ctime/localtime/
    std::time_t now; // time val
    struct tm* date; // contains components as it relates to time
    std::time(&now);
    date = std::localtime(&now);
    // https://cplusplus.com/reference/ctime/strftime/
    std::strftime(buffer, 80, "results_%B_%A_%H_%M_%S.tsv", date); // I needed a way to make unique files (date seemed appropriate)
    dataFile.open(buffer, std::ios::trunc);
    auto global_time = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < TIME_STEPS; ++i) {
        for (size_t j = 0; j < N; ++j) {
            (*particleField)[j]->clearForce();
        }
        for (size_t j = 0; j < N; ++j) {
            for (size_t k = j + 1; k < N; ++k) {
                // computing forces
                (*particleField)[j]->setForceFx((*particleField)[j]->calculateForceDirectionFx((*particleField)[k]));
                (*particleField)[j]->setForceFy((*particleField)[j]->calculateForceDirectionFy((*particleField)[k]));
                (*particleField)[j]->setForceFz((*particleField)[j]->calculateForceDirectionFz((*particleField)[k]));
                (*particleField)[j]->updateOtherFx((*particleField)[k]);
                (*particleField)[j]->updateOtherFy((*particleField)[k]);
                (*particleField)[j]->updateOtherFz((*particleField)[k]);
            }
            // calculate positions
            for (int k = 0; k < N; ++k) (*particleField)[k]->updatePosition(DELTA_T);
        }
        if(!(i % DUMP_RATE)) dumpData(N, particleField, &dataFile);
    }

    dataFile.close();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - global_time).count();
    std::ofstream logFile;
    logFile.open("data.csv", std::ios::app);
    logFile << N << "," << DELTA_T << "," << TIME_STEPS << "," << DUMP_RATE << "," << duration << std::endl;
    delete[] buffer;
}

void particleInitFromFile(const size_t& tokenIdx, std::vector<Particle*>* particleField, const size_t& classIdx, const std::string& token) {
    if (particleField == nullptr) {
        std::cerr << "Particle field pointer is NULL" << std::endl;
        exit(1);
    }
    switch (tokenIdx) {
        case 0: // we skip the first col (which is particles)
            break;
        case 1:
            (*particleField)[classIdx]->setMass(stod(token));
            break; // mass
        case 2:
            (*particleField)[classIdx]->setPositionX(stod(token));
            break; // x
        case 3:
            (*particleField)[classIdx]->setPositionY(stod(token));
            break; // y
        case 4:
            (*particleField)[classIdx]->setPositionZ(stod(token));
            break; // z
        case 5:
            (*particleField)[classIdx]->setVelocityVx(stod(token));
            break; // vx
        case 6:
            (*particleField)[classIdx]->setVelocityVy(stod(token));
            break; // vy
        case 7:
            (*particleField)[classIdx]->setVelocityVz(stod(token));
            break; // vz
        case 8:
            (*particleField)[classIdx]->setForceFx(stod(token));
            break; // fx
        case 9:
            (*particleField)[classIdx]->setForceFy(stod(token));
            break; // fy
        case 10:
            (*particleField)[classIdx]->setForceFz(stod(token));
            break; // fz
        default:
            break;
    }                                    
}

void parseFile(const std::string& fileName, const size_t& N, std::vector<Particle*>* particleField) {
    if (particleField == nullptr) {
        std::cerr << "Particle field pointer is NULL" << std::endl;
        exit(1);
    }
    for (size_t i = 0; i < N; ++i) (*particleField).push_back(new Particle());
    std::string line;
    std::ifstream sourceFile(fileName);
    if (!sourceFile.is_open()) {
        std::cerr << "Failed to open file" << std::endl;
        exit(1);
    }
    while (getline(sourceFile, line)) {
        std::stringstream s(line);
        std::string token;
        size_t tokenIdx = 0;
        size_t classIdx = 0;
        while (getline(s, token, '\t')) {
            if (tokenIdx > 10) {
                tokenIdx = 1;
                classIdx++;
            }
            particleInitFromFile(tokenIdx, particleField, classIdx, token);
            tokenIdx++;
        }
        sourceFile.close();
    }
}