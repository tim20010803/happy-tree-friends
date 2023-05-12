// C++ code
#include <iostream>
#include <queue>
#include <vector>
#include <cmath>

//constant setting
#define THETA 1.0
#define G_CONST 6.67428e-11

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
    float *monople;                    // storege the monople of all subtree
    float MaxXBoundary{0};             // the boundary of each Node (I haven't test if all the nodes own their right boundary.)
    float minXBoundary{0};
    float MaxYBoundary{0};
    float minYBoundary{0};
    std::vector<float> CalculateForce(TreeNode *Node, Particle &tarPtc); //Compute the force(acceleration) acting from this node to a particle p

    TreeNode():NW(NULL),NE(NULL),SW(NULL),SE(NULL),parent(NULL),level(0),leaf(false),ptclPtr{NULL},monople{NULL}{};                   // constuctor of TreeNode (which creats a TreeNode object and initialize parameters) 
    TreeNode(Particle *newPtcl):NW(NULL),NE(NULL),SW(NULL),SE(NULL),parent(NULL),level(0),leaf(true),ptclPtr(newPtcl),monople{NULL}{};// constuctor of TreeNode which stores the particle's location into pointer and set leaf==true (since this node is a particle)
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
        else{
            std::cout << "monople: (" << *monople<<","  << *(monople+1)<<","  << *(monople+2)<<","  << *(monople+3) <<")\n";
            std::cout << "\n\n";
            }
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
    ~QuadrupleTree();                  // desturctor (to destroy the whole tree and release the memory space it takes)
    void DeleteNode(TreeNode *Node);                      // delete the Node and its all descendent
    void Insert(Particle& newPtc);                        // insert one particle in the tree
    void Trim(TreeNode *Node);                            // delete all empty subnodes and itself if the input node turn out to be a empty (those node without any child and not a particle node)
    float *Monople(TreeNode *Node);                       // return total mass and center of mass of this subtree: (mass, x, y, z) and initialize the monople at each tree (it don't modify any monople in the subtree, it just read it.)
    void TotalForce(Particle &Ptc);                       // calculate the total force(acceleration) of a given particle
    
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
    root->MaxXBoundary = MX; root->minXBoundary = mX;     // save the boundary of root
    root->MaxYBoundary = MY; root->minYBoundary = mY;
    for (int i = 0; i < Particles.size(); i++){           // insert all particles into tree
        Insert(Particles[i]);
    }
    Monople(root);                                        //initialize all monople of nodes in whole tree
    for (int i = 0; i < Particles.size(); i++){           // insert all particles into tree
        TotalForce(Particles[i]);
    }
}
QuadrupleTree::~QuadrupleTree(){
    if (root != NULL){
        if(root->NE != NULL){DeleteNode(root->NE); root->NE = NULL;};
        if(root->NW != NULL){DeleteNode(root->NW); root->NW = NULL;};
        if(root->SE != NULL){DeleteNode(root->SE); root->SE = NULL;};
        if(root->SW != NULL){DeleteNode(root->SW); root->SW = NULL;};
        // delete root;
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
                tempOut->MaxXBoundary = tempmaxX; tempOut->minXBoundary = tempminX; //save the boundary of the node
                tempOut->MaxYBoundary = tempmaxY; tempOut->minYBoundary = tempminY;
            }
            else if (firstPtlE == false and firstPtlN){
                MX = midX;
                mY = midY;
                tempOut->NW = new TreeNode;
                tempOut->NW->parent = tempOut;
                tempOut->NW->level = tempOut->level + 1;
                tempOut = tempOut->NW;
                tempOut->MaxXBoundary = tempmaxX; tempOut->minXBoundary = tempminX;
                tempOut->MaxYBoundary = tempmaxY; tempOut->minYBoundary = tempminY;
            }
            else if (firstPtlE and firstPtlN  == false){
                mX = midX;
                MY = midY;
                tempOut->SE = new TreeNode;
                tempOut->SE->parent = tempOut;
                tempOut->SE->level = tempOut->level + 1;
                tempOut = tempOut->SE;
                tempOut->MaxXBoundary = tempmaxX; tempOut->minXBoundary = tempminX;
                tempOut->MaxYBoundary = tempmaxY; tempOut->minYBoundary = tempminY;
            }
            else if (firstPtlE == false and firstPtlN == false){
                MX = midX;
                MY = midY;
                tempOut->SW = new TreeNode;
                tempOut->SW->parent = tempOut;
                tempOut->SW->level = tempOut->level + 1;
                tempOut = tempOut->SW;
                tempOut->MaxXBoundary = tempmaxX; tempOut->minXBoundary = tempminX;
                tempOut->MaxYBoundary = tempmaxY; tempOut->minYBoundary = tempminY;
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
                ptcTree1->MaxXBoundary = tempmaxX; ptcTree1->minXBoundary = tempminX; // save the boundary of the node
                ptcTree1->MaxYBoundary = tempmaxY; ptcTree1->minYBoundary = tempminY;
            }
            else if (firstPtlE == false and firstPtlN){
                tempOut->NW = ptcTree1;
                ptcTree1->parent = tempOut;
                ptcTree1->level = tempOut->level +1;
                ptcTree1->MaxXBoundary = tempmaxX; ptcTree1->minXBoundary = tempminX;
                ptcTree1->MaxYBoundary = tempmaxY; ptcTree1->minYBoundary = tempminY;
            }
            else if (firstPtlE and firstPtlN  == false){
                tempOut->SE = ptcTree1;
                ptcTree1->parent = tempOut;
                ptcTree1->level = tempOut->level +1;
                ptcTree1->MaxXBoundary = tempmaxX; ptcTree1->minXBoundary = tempminX;
                ptcTree1->MaxYBoundary = tempmaxY; ptcTree1->minYBoundary = tempminY;
            }
            else if (firstPtlE == false and firstPtlN == false){
                tempOut->SW = ptcTree1;
                ptcTree1->parent = tempOut;
                ptcTree1->level = tempOut->level +1;
                ptcTree1->MaxXBoundary = tempmaxX; ptcTree1->minXBoundary = tempminX;
                ptcTree1->MaxYBoundary = tempmaxY; ptcTree1->minYBoundary = tempminY;
            }                                            // creat node for second paticle and connect it to tempOut
            if (secPtlE and secPtlN){
                tempOut->NE = new TreeNode(&ptc2);
                tempOut->NE->parent = tempOut;
                tempOut->NE->level = tempOut->level + 1;
                tempOut->NE->leaf = true;
                tempOut->NE->MaxXBoundary = tempmaxX; tempOut->NE->minXBoundary = tempminX;
                tempOut->NE->MaxYBoundary = tempmaxY; tempOut->NE->minYBoundary = tempminY;
            }
            else if (secPtlE == false and secPtlN){
                tempOut->NW = new TreeNode(&ptc2);
                tempOut->NW->parent = tempOut;
                tempOut->NW->level = tempOut->level + 1;
                tempOut->NW->leaf = true;
                tempOut->NW->MaxXBoundary = tempmaxX; tempOut->NW->minXBoundary = tempminX;
                tempOut->NW->MaxYBoundary = tempmaxY; tempOut->NW->minYBoundary = tempminY;
            }
            else if (secPtlE and secPtlN  == false){
                tempOut->SE = new TreeNode(&ptc2);
                tempOut->SE->parent = tempOut;
                tempOut->SE->level = tempOut->level + 1;
                tempOut->SE->leaf = true;
                tempOut->SE->MaxXBoundary = tempmaxX; tempOut->SE->minXBoundary = tempminX;
                tempOut->SE->MaxYBoundary = tempmaxY; tempOut->SE->minYBoundary = tempminY;
            }
            else if (secPtlE == false and secPtlN == false){
                tempOut->SW = new TreeNode(&ptc2);
                tempOut->SW->parent = tempOut;
                tempOut->SW->level = tempOut->level + 1;
                tempOut->SW->leaf = true;
                tempOut->SW->MaxXBoundary = tempmaxX; tempOut->SW->minXBoundary = tempminX;
                tempOut->SW->MaxYBoundary = tempmaxY; tempOut->SW->minYBoundary = tempminY;
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
            current->NW->MaxXBoundary = tempMaxX; current->NW->minXBoundary = tempMinX;  //save the boundary of the node
            current->NW->MaxYBoundary = tempMaxY; current->NW->minYBoundary = tempMinY;
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
            current->NE->MaxXBoundary = tempMaxX; current->NE->minXBoundary = tempMinX;
            current->NE->MaxYBoundary = tempMaxY; current->NE->minYBoundary = tempMinY;
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
            current->SW->MaxXBoundary = tempMaxX; current->SW->minXBoundary = tempMinX;
            current->SW->MaxYBoundary = tempMaxY; current->SW->minYBoundary = tempMinY;
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
            current->SE->MaxXBoundary = tempMaxX; current->SE->minXBoundary = tempMinX;
            current->SE->MaxYBoundary = tempMaxY; current->SE->minYBoundary = tempMinY;
        }
        tempMidX = (tempMaxX + tempMinX) / 2.;          // uodate new mid-line for next run
        tempMidY = (tempMaxY + tempMinY) / 2.;
        current = q.front();                            // change current node to the next (deeper node)
        q.pop();
    }
}
void QuadrupleTree::DeleteNode(TreeNode *Node){
    if(Node != root){           // cut the pointer (which is point to self from its parent)
        if (Node->parent->monople != NULL){//unfinished
            for (int i = 1; i < 3; i++){
                *(Node->parent->monople + i) *= *(Node->parent->monople);
                *(Node->parent->monople + i) -= (*(Node->monople+i) * (*(Node->monople))  ); 
            }
            *(Node->parent->monople) -= *(Node->monople);
            for (int i = 1; i < 3; i++){
                *(Node->parent->monople + i) /= *(Node->parent->monople);
            }
        }
        if(Node->parent->NE == Node){Node->parent->NE = NULL;};
        if(Node->parent->NW == Node){Node->parent->NW = NULL;};
        if(Node->parent->SE == Node){Node->parent->SE = NULL;};
        if(Node->parent->SW == Node){Node->parent->SW = NULL;};
    }
    if (Node->leaf){            // if the node is a particle node, disconnect to particle and delete self
        Node->ptclPtr = NULL;
        delete Node->monople;
        delete Node;
        return;
    }
    else{                       // delete all the children
        if(Node->NE != NULL){DeleteNode(Node->NE);};
        if(Node->NW != NULL){DeleteNode(Node->NW);};
        if(Node->SE != NULL){DeleteNode(Node->SE);};
        if(Node->SW != NULL){DeleteNode(Node->SW);};
        delete Node->monople;
        Node->monople = NULL;
        delete Node;
        return;
    }
}
void QuadrupleTree::Trim(TreeNode *Node){
    if(Node->leaf){return;}            // if this node is a particle return
    if((Node->leaf == false and Node->NE == NULL) and (Node->NW == NULL) and (Node->SE == NULL) and (Node->SW == NULL)){
        DeleteNode(Node);
        return;
    }
    else if(Node->leaf == false){      // trim all the children
        if(Node->NE != NULL){Trim(Node->NE);};
        if(Node->NW != NULL){Trim(Node->NW);};
        if(Node->SE != NULL){Trim(Node->SE);};
        if(Node->SW != NULL){Trim(Node->SW);};
        Trim(Node);
        return;
    }
}
float *QuadrupleTree::Monople(TreeNode *Node){
    if (Node->monople != NULL)
    {
        return Node->monople;
    }
    float *monoParaPtr = new float;     // {mass, x, y, z}
    for (int i = 0; i < 3; i++){        // initialize the value which pointer point to
        *(monoParaPtr + i)= 0.;        
    }
    if (Node->leaf == true){            // if this is a leaf(particle node) output the particle
        *monoParaPtr = Node->ptclPtr->mass;
        *(monoParaPtr + 1) = Node->ptclPtr->posi[0];
        *(monoParaPtr + 2) = Node->ptclPtr->posi[1];
        // *(monoParaPtr + 3) = Node->ptclPtr->posi[2]; // for z direction 
    }
    else{                               // if this node has children then compute the monople of children
        if(Node->NE != NULL){
            float *monoNEptr = Monople(Node->NE);
            *monoParaPtr += *monoNEptr;
            *(monoParaPtr + 1) += *(monoNEptr + 1) * (*monoNEptr); // add center of mass times mass  
            *(monoParaPtr + 2) += *(monoNEptr + 2) * (*monoNEptr); // add center of mass times mass  
            // *(monoParaPtr + 3) += *(monoNEptr + 3) * (*monoNEptr); // add center of mass times mass  
            // delete[] monoNEptr;
        }
        if(Node->NW != NULL){
            float *monoNWptr = Monople(Node->NW);
            *monoParaPtr += *monoNWptr;
            *(monoParaPtr + 1) += *(monoNWptr + 1) * (*monoNWptr); // add center of mass times mass  
            *(monoParaPtr + 2) += *(monoNWptr + 2) * (*monoNWptr); // add center of mass times mass  
            // *(monoParaPtr + 3) += *(monoNWptr + 3) * (*monoNWptr); // add center of mass times mass  
            // delete[] monoNWptr;
        }
        if(Node->SE != NULL){
            float *monoSEptr = Monople(Node->SE);
            *monoParaPtr += *monoSEptr;
            *(monoParaPtr + 1) += *(monoSEptr + 1) * (*monoSEptr); // add center of mass times mass  
            *(monoParaPtr + 2) += *(monoSEptr + 2) * (*monoSEptr); // add center of mass times mass  
            // *(monoParaPtr + 3) += *(monoSEptr + 3) * (*monoSEptr); // add center of mass times mass  
            // delete[] monoSEptr;
        }
        if(Node->SW != NULL){
            float *monoSWptr = Monople(Node->SW);
            *monoParaPtr += *monoSWptr;
            *(monoParaPtr + 1) += *(monoSWptr + 1) * (*monoSWptr); // add center of mass times mass  
            *(monoParaPtr + 2) += *(monoSWptr + 2) * (*monoSWptr); // add center of mass times mass  
            // *(monoParaPtr + 3) += *(monoSWptr + 3) * (*monoSWptr); // add center of mass times mass  
            // delete[] monoSWptr;
        }
        for (int i = 1; i < 3; i++){        // normalize center of mass with total mass
            *(monoParaPtr + i) /= (*monoParaPtr);        
        }
    }
    Node->monople = monoParaPtr;
    return monoParaPtr;
}


//calculate the force(acceleration) acting from Node to the target particle
std::vector<float> TreeNode::CalculateForce(TreeNode *Node, Particle &tarPtc){

    float r{0};float a{0}; 
    std::vector<float> acc;
    // r is the distance between the Node and the particle
    r = sqrt((tarPtc.posi[0] - *(Node->monople +1))*(tarPtc.posi[0] - *(Node->monople +1))+(tarPtc.posi[1] - *(Node->monople +2))*(tarPtc.posi[1] - *(Node->monople +2)));
    // a is the acceleration of the particle using the Newton's law
    // a is divided by one more r because the next step is divided by one less r
    a = G_CONST * *(Node->monople) / (r * r * r);
    acc.push_back( a * (tarPtc.posi[0] - *(Node->monople +1))); // the acceleration of x component ( a * x/r )
    acc.push_back( a * (tarPtc.posi[1] - *(Node->monople +2))); // the acceleration of y component ( a * y/r )

    return acc; // include x and y axis
}

//calculate the total force acting on the given particle
void QuadrupleTree::TotalForce(Particle &Ptc){
    std::queue<TreeNode*> q;
    TreeNode *current = root;
    std::vector<float> accSum {0.0, 0.0}; // the sum of acceleration, including x and y components
    std::vector<TreeNode*>section; // change NE, NW, SE, SW into a vector to use for-loop
    std::vector<float> acc; // to store the computing acceleration
    float r{0};float d{0};
    while (current){
        r = sqrt((Ptc.posi[0] - *(current->monople +1))*(Ptc.posi[0] - *(current->monople +1))+(Ptc.posi[1] - *(current->monople +2))*(Ptc.posi[1] - *(current->monople +2)));
        d = current->MaxXBoundary  - current->minXBoundary;
        // this is the node where the particle itself exists
        if (r <= 0){
            break; 
        }
        // test if Multipole-Acceptance-Criterion(MAC) can be used in this node
        if (d / r <= THETA){
            acc = current->CalculateForce(current, Ptc); // calculate the force(acceleration)
            accSum[0] += acc[0]; // x component
            accSum[1] += acc[1]; // y component
        }
        else{
            section.push_back(current->NE);// to store the computing acceleration
            section.push_back(current->NW);
            section.push_back(current->SE);
            section.push_back(current->SW);
            // look down and check if MAC can be used in the children nodes
            for (int i = 0; i < 4; i++){
                if(section[i] != NULL){
                    q.push(section[i]);
                }
            }
            section.clear();
        current = q.front();
        q.pop();
        }
    }
    Ptc.acceleration[0] = accSum[0];
    Ptc.acceleration[1] = accSum[1];
    return;
}


// main function is for testing
int main() {
    Particle a;
    a.posi = {1.,6.,4.};
    a.velocity = {4.,3.,2.};
    a.mass = {12.};
    a.acceleration = {0., 0., 0.};
    Particle b;
    b.posi = {2.,7.,8.};
    b.velocity = {1.,6.,7.};
    b.mass = {23.};
    b.acceleration = {0., 0., 0.};
    Particle c;
    c.posi = {3.,5.,8.};
    c.velocity = {1.,6.,7.};
    c.mass = {212.};
    c.acceleration = {0., 0., 0.};
    Particle d;
    d.posi = {4.,6.,8.};
    d.velocity = {8.,6.,7.};
    d.mass = {62.};
    d.acceleration = {0., 0., 0.};
    
    std::vector<Particle> Pvec = {a,b,c,d};
    QuadrupleTree T(Pvec,0.,0.,0.,10.,10.,10.); 
    
    //std::cout << "Particle a: " <<  a.acceleration[0] <<  "\n";

    T.root->PrintNode();
    T.root->SW->PrintNode();
    T.root->NW->PrintNode();
    T.root->NW->SE->PrintNode();
    T.root->NW->SW->PrintNode();
    T.root->NW->SW->SW->PrintNode();
    T.root->NW->SW->NE->PrintNode();
    T.root->NW->PrintNode();
    return 0;
}