// C++ code
#include <iostream>
#include <sstream>
#include <queue>
#include <vector>
#include <string>
struct Particle{
    std::vector<float> posi;
    std::vector<float> velocity;
    std::vector<float> acceleration;
    float mass;
};

class QuadrupleTree;
class TreeNode{
// private:
public:                                // set all of parameter be public for testing. you should only access all those parameters by function in QuadrupleTree.
    TreeNode *parent;                  // point to the parent node
    TreeNode *NW;                      // northwest child node
    TreeNode *NE;                      // northeast child node
    TreeNode *SW;                      // southwest child node
    TreeNode *SE;                      // southeast child node
    Particle *ptclPtr;                 // the pointer of single particle if this node is a leaf(without any child node), otherwise it's a null pointer.  
    int level;                         // the depth level of this node (which is zero for root node)
    bool leaf;                         // if leaf is true then this node is a leaf the point toward the particle is ptclPtr

    TreeNode():NW(NULL),NE(NULL),SW(NULL),SE(NULL),parent(NULL),level(0),leaf(false),ptclPtr{NULL}{};                   // constuctor of TreeNode (which creats a TreeNode object and initialize parameters) 
    TreeNode(Particle *newPtcl):NW(NULL),NE(NULL),SW(NULL),SE(NULL),parent(NULL),level(0),leaf(true),ptclPtr(newPtcl){};// constuctor of TreeNode which stores the particle's location into pointer and set leaf==true (since this node is a particle)
    void PrintNode(){                                                                                    // print the information of node parent, self and children are the memeory locations of each node
        std::cout << "parent: " <<  parent <<  "\n";
        std::cout << "self: " <<  this <<  "\n";
        std::cout << "children: " <<"NW:"<< NW <<" NE:"<< NE <<" SW:"<< SW <<" SE:"<< SE << "\n";
        std::cout << "level: " << level << "\n";
        if (leaf == true){                                                                               // If the node is a particle, print the parameter of it
            std::cout << "this node is a Particle\n";
            std::cout << "position: (" << ptclPtr->posi[0] << "," << ptclPtr->posi[1] << "," << ptclPtr->posi[2] << ")\n";
            std::cout << "velocity: (" << ptclPtr->velocity[0] << "," << ptclPtr->velocity[1] << ","<<ptclPtr->velocity[2] << ")\n";
            std::cout << "mass: " << ptclPtr->mass << "\n\n";
        }
        else{std::cout << "\n\n";}
        };
        friend class QuadrupleTree;   // give the QuadrupleTree class access to those private members(such as private functions and parameters)
};

class QuadrupleTree{
private:

    float minX{0};float minY{0};float minZ{0}; // the boundary of space (minimum and maximum of x,y,z of all particles)
    float maxX{1};float maxY{1};float maxZ{1};
    TreeNode *TwoParticleSubtree(TreeNode *ptcTree1, Particle &ptc2,float tempminX,float tempminY,float tempminZ,float tempmaxX,float tempmaxY, float tempmaxZ); 
    // TwoParticleSubtree function creats a tree which connects two node particle,
    // if a particle in the same region of existing particle in a node 
    // in the original level of first particle
public:

    TreeNode *root;
    QuadrupleTree():root(NULL){};
    QuadrupleTree(Particle &firstPtc,float mX,float mY,float mZ,float MX,float MY, float MZ);                // take one particle and boundary of all particle to ininitaize the tree, and mX is minX(minimum x), MX is maxX(maxmum x);
    QuadrupleTree(std::vector<Particle> &Particles,float mX,float mY,float mZ,float MX,float MY, float MZ);  // take several particles(with type of std::vector) and boundary of all particle to ininitaize the tree, and mX is minX(minimum x), MX is maxX(maxmum x);
    ~QuadrupleTree(){delete []root;};                     // desturctor (to destroy the whole tree and release the memory space it takes)
    void Insert(Particle& newPtc);                        // insert one particle in the tree
    float *Monople(TreeNode *subtreeNode);              // total mass and center of mass of this subtree: (mass, x, y, z) (not done yet)
};
QuadrupleTree::QuadrupleTree(Particle &firstPtc,float mX,float mY,float mZ,float MX,float MY, float MZ){                 
    root = new TreeNode;                                  // allocate memory for root
    maxX = MX; maxY = MY; maxZ = MZ;                      // inintialize boundary of this tree
    minX = mX; minY = mY; minZ = mZ; 
    Insert(firstPtc);
}
QuadrupleTree::QuadrupleTree(std::vector<Particle> &Particles,float mX,float mY,float mZ,float MX,float MY, float MZ){  
    root = new TreeNode;                                  // allocate memory for root
    maxX = MX; maxY = MY; maxZ = MZ;                      // inintialize boundary of this tree
    minX = mX; minY = mY; minZ = mZ; 
    for (int i = 0; i < Particles.size(); i++){           // insert all particles into tree
        Insert(Particles[i]);
    }
}

TreeNode *QuadrupleTree::TwoParticleSubtree(TreeNode *ptcTree1, Particle &ptc2,float tempminX,float tempminY,float tempminZ,float tempmaxX,float tempmaxY, float tempmaxZ){
    float mX = tempminX;float mY = tempminY;float mZ = tempminZ;        // rename new boundary and cut it into four subregion
    float MX = tempmaxX;float MY = tempmaxY;float MZ = tempmaxZ;
    float midX = (mX + MX)/2.;float midY = (mY + MY)/2.;float midZ = (mZ + MZ)/2.;

    bool firstPtlN = ptcTree1->ptclPtr->posi[1] > midY;                 // determine each Particle is in which part of subregion
    bool firstPtlE = ptcTree1->ptclPtr->posi[0] > midX;
    bool secPtlN = ptc2.posi[1] > midY;
    bool secPtlE = ptc2.posi[0] > midX;

    TreeNode *output = new TreeNode;                                    // allocate a new node for two particle node
    output->level = ptcTree1->level - 1;                                // set the levels of output node and particle one node  
    ptcTree1->level += 1;

    TreeNode *tempOut = output;                                         // new TreeNode pointer for while loop
    while (true){                                                       // For each run, determine if two particle are in the same subregion. 
        if ( (firstPtlE == secPtlE) and (firstPtlN == secPtlN)){        // If two particles are in the same subregion, divide the region and go to the next run. 
            tempOut->leaf = false;
            tempOut->ptclPtr = NULL;
            if (firstPtlE and firstPtlN){                               // if two particles are in the northeast region 
                mX = midX;                                              // redraw the boundary of region into its subregion
                mY = midY;
                tempOut->NE = new TreeNode;                             // creats a new node and connect it to correct pointer of given subregion
                tempOut->NE->parent = tempOut;                          // update parameters of new node
                tempOut->NE->level = tempOut->level + 1;
                tempOut = tempOut->NE;                                  // move tempnode to new node for next run
            }
            else if (firstPtlE == false and firstPtlN){
                MX = midX;
                mY = midY;
                tempOut->NW = new TreeNode;
                tempOut->NW->parent = tempOut;
                tempOut->NW->level = tempOut->level + 1;
                tempOut = tempOut->NW;
            }
            else if (firstPtlE and firstPtlN  == false){
                mX = midX;
                MY = midY;
                tempOut->SE = new TreeNode;
                tempOut->SE->parent = tempOut;
                tempOut->SE->level = tempOut->level + 1;
                tempOut = tempOut->SE;
            }
            else if (firstPtlE == false and firstPtlN == false){
                MX = midX;
                MY = midY;
                tempOut->SW = new TreeNode;
                tempOut->SW->parent = tempOut;
                tempOut->SW->level = tempOut->level + 1;
                tempOut = tempOut->SW;
            }
            midX = (mX + MX)/2.; midY = (mY + MY)/2.; midZ = (mZ + MZ)/2.;    // update new mid-line of new region 
            firstPtlN = ptcTree1->ptclPtr->posi[1] > midY;                    // determin each Particle in which part of subregion
            firstPtlE = ptcTree1->ptclPtr->posi[0] > midX;
            secPtlN = ptc2.posi[1] > midY;
            secPtlE = ptc2.posi[0] > midX;
        }
        else{                                           // if two particles are not in the same region                                  
            tempOut->leaf = false;
            tempOut->ptclPtr = NULL;
            if (firstPtlE and firstPtlN){               // connect ptcTree1 to tempOut
                tempOut->NE = ptcTree1;
                ptcTree1->parent = tempOut;
                ptcTree1->level = tempOut->level +1;
            }
            else if (firstPtlE == false and firstPtlN){
                tempOut->NW = ptcTree1;
                ptcTree1->parent = tempOut;
                ptcTree1->level = tempOut->level +1;
            }
            else if (firstPtlE and firstPtlN  == false){
                tempOut->SE = ptcTree1;
                ptcTree1->parent = tempOut;
                ptcTree1->level = tempOut->level +1;

            }
            else if (firstPtlE == false and firstPtlN == false){
                tempOut->SW = ptcTree1;
                ptcTree1->parent = tempOut;
                ptcTree1->level = tempOut->level +1;
            }                                            // creat node for second paticle and connect it to tempOut
            if (secPtlE and secPtlN){
                tempOut->NE = new TreeNode(&ptc2);
                tempOut->NE->parent = tempOut;
                tempOut->NE->level = tempOut->level + 1;
                tempOut->NE->leaf = true;
            }
            else if (secPtlE == false and secPtlN){
                tempOut->NW = new TreeNode(&ptc2);
                tempOut->NW->parent = tempOut;
                tempOut->NW->level = tempOut->level + 1;
                tempOut->NW->leaf = true;
            }
            else if (secPtlE and secPtlN  == false){
                tempOut->SE = new TreeNode(&ptc2);
                tempOut->SE->parent = tempOut;
                tempOut->SE->level = tempOut->level + 1;
                tempOut->SE->leaf = true;
            }
            else if (secPtlE == false and secPtlN == false){
                tempOut->SW = new TreeNode(&ptc2);
                tempOut->SW->parent = tempOut;
                tempOut->SW->level = tempOut->level + 1;
                tempOut->SW->leaf = true;
            }
            break;
        }
    }
    return output;
}

void QuadrupleTree::Insert(Particle& newPtc){    
    std::queue<TreeNode*> q;                                        // store node in first-in-first-out order to search 
    TreeNode *current = root;
    float tempMaxX{maxX};float tempMaxY{maxY};float tempMaxZ{maxZ}; // copy thw boundary of tree (as search gets deeper,
    float tempMinX{minX};float tempMinY{minY};float tempMinZ{minZ}; // the boundary shrinks for determining which node should be search)
    float tempMidX = (tempMaxX + tempMinX) / 2.;                    // mid-line of boundary (changes as boundary being changed)
    float tempMidY = (tempMaxY + tempMinY) / 2.;
    while (current) {
        if (newPtc.posi[0] < (tempMidX) and newPtc.posi[1] > (tempMidY)){ // the new particle is in northwest subregion
            if (current->NW != NULL ){                  // if current node's child node NW already exists
                if (current->NW->leaf == false){        // if this child node is not leaf (which means it's not a particle node and current node has grandchild) 
                    q.push(current->NW);                // push it to the tail of queue (search will continue, starting from this node)
                    current->leaf = false;              // mark that current node is not a particle
                }
                else{                                   // current node's child(NW) is a particle node(two particle are in the same region in this level)
                    current->NW = TwoParticleSubtree(current->NW, newPtc, tempMinX,tempMinY,tempMinZ,tempMaxX,tempMaxY,tempMaxZ)->NW; // creat a subtree contains both particles and connect it to the original node
                    delete [] current->NW->parent;
                    current->NW->parent = current;
                    break;
                }
            }
            else{                                       // current node does not have this(NW) child node(empty)
                TreeNode *new_node = new TreeNode(&newPtc);   // creat a new node and put new particle's pointer in it then connect it to the current node (insertion complete)
                new_node->parent = current;
                new_node->level = current->level + 1;
                new_node->leaf = true;
                current->leaf = false;
                current->NW = new_node;
                break;                         
            }
            tempMaxX = (tempMidX);                      // update boundary of new region for next run
            tempMinY = (tempMidY);
        }
        else if (newPtc.posi[0] >= (tempMidX) and newPtc.posi[1] > (tempMidY)){ // the new particle is in northeast subregion
            if (current->NE != NULL ){                  // if current node's child node NE already exists
                if (current->NE->leaf == false){        // if this child node is not leaf (which means it's not a particle node and current node has grandchild) 
                    q.push(current->NE);                // push it to the tail of queue (search will continue, starting from this node)
                    current->leaf = false;              // mark that current node is not a particle
                }
                else{                                   // current node's child(NE) is a particle node(two particle are in the same region in this level)
                    current->NE = TwoParticleSubtree(current->NE, newPtc, tempMinX,tempMinY,tempMinZ,tempMaxX,tempMaxY,tempMaxZ)->NE; // creat a subtree contains both particles and connect it to the original node
                    delete [] current->NE->parent;
                    current->NE->parent = current;
                    break;
                }
            }
            else{                                       // current node does not have this(NE) child node(empty)
                TreeNode *new_node = new TreeNode(&newPtc);   // creat a new node and put new particle's pointer in it then connect it to the current node (insertion complete)
                new_node->parent = current;
                new_node->level = current->level + 1;
                new_node->leaf = true;
                current->leaf = false;
                current->NE = new_node;
                break;
            }
            tempMinX = (tempMidX);                      // update boundary of new region for next run
            tempMinY = (tempMidY);
        }
        else if (newPtc.posi[0] < (tempMidX) and newPtc.posi[1] <= (tempMidY)){ // the new particle is in southwest subregion
            if (current->SW != NULL ){                  // if current node's child node SW already exists
                if (current->SW->leaf == false){        // if this child node is not leaf (which means it's not a particle node and current node has grandchild) 
                    q.push(current->SW);                // push it to the tail of queue (search will continue, starting from this node)
                    current->leaf = false;              // mark that current node is not a particle
                }
                else{                                   // current node's child(SW) is a particle node(two particle are in the same region in this level)
                    current->SW = TwoParticleSubtree(current->SW, newPtc, tempMinX,tempMinY,tempMinZ,tempMaxX,tempMaxY,tempMaxZ)->SW; // creat a subtree contains both particles and connect it to the original node
                    delete [] current->SW->parent;
                    current->SW->parent = current;
                    break;
                }
            }
            else{                                       // current node does not have this(SW) child node(empty)
                TreeNode *new_node = new TreeNode(&newPtc);   // creat a new node and put new particle's pointer in it then connect it to the current node (insertion complete)
                new_node->parent = current;
                new_node->level = current->level + 1;
                new_node->leaf = true;
                current->leaf = false;
                current->SW = new_node;
                break;                
            }
            tempMaxX = (tempMidX);                      // update boundary of new region for next run
            tempMaxY = (tempMidY);
        }
        else if (newPtc.posi[0] >= (tempMidX) and newPtc.posi[1] <= (tempMidY)){ // the new particle is in southeast subregion
            if (current->SE != NULL ){                  // if current node's child node SE already exists且該分支不是particle
                if (current->SE->leaf == false){        // if this child node is not leaf (which means it's not a particle node and current node has grandchild) 
                    q.push(current->SE);                // push it to the tail of queue (search will continue, starting from this node)
                    current->leaf = false;              // mark that current node is not a particle
                }
                else{                                   // current node's child(SE) is a particle node(two particle are in the same region in this level)
                    current->SE = TwoParticleSubtree(current->SE, newPtc, tempMinX,tempMinY,tempMinZ,tempMaxX,tempMaxY,tempMaxZ)->SE; // creat a subtree contains both particles and connect it to the original node
                    delete [] current->SE->parent;
                    current->SE->parent = current;
                    break;
                }
            }
            else{                                       // current node does not have this(SE) child node(empty)
                TreeNode *new_node = new TreeNode(&newPtc);   // creat a new node and put new particle's pointer in it then connect it to the current node (insertion complete)
                new_node->parent = current;
                new_node->level = current->level + 1;
                new_node->leaf = true;
                current->leaf = false;
                current->SE = new_node;
                break;                         
            }
            tempMinX = (tempMidX);                      // update boundary of new region for next run
            tempMaxY = (tempMidY);
        }
        tempMidX = (tempMaxX + tempMinX) / 2.;          // uodate new mid-line for next run
        tempMidY = (tempMaxY + tempMinY) / 2.;
        current = q.front();                            // change current node to the next (deeper node)
        q.pop();
    }
}









// main function is for testing
int main() {
    Particle a;
    a.posi = {1.,6.,4.};
    a.velocity = {4.,3.,2.};
    a.mass = {12.};
    Particle b;
    b.posi = {2.,7.,8.};
    b.velocity = {1.,6.,7.};
    b.mass = {23.};
    Particle c;
    c.posi = {3.,5.,8.};
    c.velocity = {1.,6.,7.};
    c.mass = {212.};
    Particle d;
    d.posi = {4.,6.,8.};
    d.velocity = {8.,6.,7.};
    d.mass = {62.};
    // TreeNode* AnodePtr = new TreeNode (&a);
    // TreeNode* BnodePtr = new TreeNode (&b);
    // TreeNode* CnodePtr = new TreeNode (&c);
    // AnodePtr->PrintNode();



    // TreeNode* BinA = TwoParticleSubtree(AnodePtr,b,0.,0.,0.,10.,10.,10.);
    // std::cout<<"\n\n" ;
    // if (BinA != NULL)
    // {
    // BinA->PrintNode();
    //     if (BinA->SE != NULL)
    //     {
    //         std::cout <<"SE:\n";
    //         BinA->SE->PrintNode();
    //     }
    //     if (BinA->NE != NULL)
    //     {
    //         std::cout <<"NE:\n";
    //         BinA->NE->PrintNode();
    //     }
    //     if (BinA->SW != NULL)
    //     {
    //         std::cout <<"SW:\n";
    //         BinA->SW->PrintNode();
    //     }
    //     if (BinA->NW != NULL)
    //     {
    //         std::cout <<"NW:\n";
    //         BinA->NW->PrintNode();
    //         BinA->NW->SW->PrintNode();
    //         BinA->NW->SW->NE->PrintNode();
    //         BinA->NW->SW->SW->PrintNode();
    //     }

    // }
    
    
    // QuadrupleTree T(a,0.,0.,0.,10.,10.,10.); 


    // T.Insert(b);
    // T.Insert(c);
    // T.Insert(d);
    // T.root->PrintNode();
    // T.root->SW->PrintNode();
    // T.root->NW->PrintNode();
    // T.root->NW->SE->PrintNode();
    // T.root->NW->SW->PrintNode();
    // T.root->NW->SW->SW->PrintNode();
    // T.root->NW->SW->NE->PrintNode();
    std::vector<Particle> Pvec = {a,b,c,d};
    QuadrupleTree T(Pvec,0.,0.,0.,10.,10.,10.); 



    T.root->PrintNode();
    T.root->SW->PrintNode();
    T.root->NW->PrintNode();
    T.root->NW->SE->PrintNode();
    T.root->NW->SW->PrintNode();
    T.root->NW->SW->SW->PrintNode();
    T.root->NW->SW->NE->PrintNode();
    return 0;
}
