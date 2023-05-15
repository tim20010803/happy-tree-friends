#include <iostream>
#include <vector>
#include <string>

#define X_MIN 100
#define Y_MIN 100
#define X_MAX 200
#define Y_MAX 200

enum EQuadrant{
    NE = 0,
    NW,
    SW,
    SE,
    NONE
};



void BuildTree() {
    //ResetTree();

    for (auto& particle : particles) {
        rootNode->InsertToNode(particle);
    }
}

class Particle {
public:
    std::vector<float> posi;
    std::vector<float> velocity;
    float mass;
};

class TreeNode {
public:
    TreeNode() : numParticles(0), particle(nullptr) {
        for (int i = 0; i < 4; ++i) {
            subnodes[i] = nullptr;
        }
        numParticles = 0;
        parent = nullptr;
    }

    void InsertToNode(Particle* newParticle) {
        if (numParticles > 1) {
            int quadrant = GetQuadrant(newParticle);
            if (subnodes[quadrant] == nullptr) {
                subnodes[quadrant] = new TreeNode;
                subnodes[quadrant]->Chooseboundary(quadrant);
            }
            subnodes[quadrant]->InsertToNode(newParticle);
        }
        else if (numParticles == 1) {
            int quadrant = GetQuadrant(particle);
            if (subnodes[quadrant] == nullptr) {
                subnodes[quadrant] = new TreeNode;
                subnodes[quadrant]->Chooseboundary(quadrant);
                parent = this;
            }
            subnodes[quadrant]->InsertToNode(particle);

            quadrant = GetQuadrant(newParticle);
            if (subnodes[quadrant] == nullptr) {
                subnodes[quadrant] = new TreeNode;
                subnodes[quadrant]->Chooseboundary(quadrant);
            }
            subnodes[quadrant]->InsertToNode(newParticle);
            numParticles = 2;
            particle = nullptr;
        }
        else if (numParticles == 0) {
            particle = newParticle;
            numParticles = 1;
        }
    }

private:
    int numParticles;
    Particle* particle;
    TreeNode* subnodes[4];
    TreeNode* parent;
    float mx, my, Mx, My;
    float midx = (mx+Mx)/2;
    float midy = (my+My)/2;

    void Chooseboundary(int quadrant){

        switch(quadrant)
        {
            case 0:
                subnodes[quadrant]->Mx = parent->midx;
                subnodes[quadrant]->My = parent->midy;
                break;
            case 1:
                subnodes[quadrant]->mx = parent->midx;
                subnodes[quadrant]->My = parent->midy;
                break;
            case 2:
                subnodes[quadrant]->Mx = parent->midx;
                subnodes[quadrant]->my = parent->midy;
                break;
            case 3:
                subnodes[quadrant]->mx = parent->midx;
                subnodes[quadrant]->my = parent->midy;
                break;
            }

    }

    int GetQuadrant(Particle* particle) {
        
        if (particle->posi[0] <= midx && particle->posi[1] <= midy)
        {
            return SW;
        }
        else if (particle->posi[0] <= midx && particle->posi[1] >= midy)
        {
            return NW;
        }
        else if (particle->posi[0] >= midx && particle->posi[1] >= midy)
        {
            return NE;
        }
        else if (particle->posi[0] >= midx && particle->posi[1] <= midy)
        {
            return SE;
        }
        return 0;
    }
};



int main() {

    TreeNode* rootNode;
    std::vector<Particle> particles;
    rootNode = new TreeNode;

    rootNode->mx = X_MIN;
    rootNode->my = Y_MIN;
    rootNode->Mx = X_MAX;
    rootNode->My = Y_MAX;

    BuildTree();
-
    
/*
    void ResetTree() {

    }
    */
    return 0;
};
