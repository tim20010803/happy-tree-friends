// C++ code
#include <iostream>
#include <queue>
#include <vector>
#include <cmath>
#include <iomanip>
#include "orbit_integration.h"
#include "quadrupleTree.h"
//constant setting
// #define THETA 1.0

TreeNode::~TreeNode(){
    if(parent!=0 and level>0){
        if(parent->NE == this){ parent->NE = NULL;};
        if(parent->NW == this){ parent->NW = NULL;};
        if(parent->SE == this){ parent->SE = NULL;};
        if(parent->SW == this){ parent->SW = NULL;};
        parent = NULL;};
    if(NE != NULL){(delete NE);NE = NULL;};
    if(NW != NULL){(delete NW);NW = NULL;};
    if(SE != NULL){(delete SE);SE = NULL;};
    if(SW != NULL){(delete SW);SW = NULL;};
    if(ptclPtr != NULL){ptclPtr = NULL;};
    if(monople != NULL){(delete monople);monople = NULL;};
}
QuadrupleTree::QuadrupleTree(Particle &firstPtc,double mX,double mY,double mZ,double MX,double MY, double MZ){                 
    root = new TreeNode;                                  // allocate memory for root
    root->parent = nullptr;
    root->level =0;
    maxX = MX; maxY = MY; maxZ = MZ;                      // inintialize boundary of this tree
    minX = mX; minY = mY; minZ = mZ; 
    Insert(firstPtc);
}
QuadrupleTree::QuadrupleTree(std::vector<Particle> &Particles,double mX,double mY,double mZ,double MX,double MY, double MZ){  
    PtcVectorPtr = &Particles;
    root = new TreeNode;                                  // allocate memory for root
    root->parent = nullptr;
    root->level =0;
    maxX = MX; maxY = MY; maxZ = MZ;                      // inintialize boundary of this tree
    minX = mX; minY = mY; minZ = mZ;
    for (int i = 0; i < Particles.size(); i++){           // insert all particles into tree
        if(Particles[i].posi[0]>MX or Particles[i].posi[0]<mX or Particles[i].posi[1]>MY or Particles[i].posi[1]<mY){
            std::cout<<"particle out of boundry";
            std::exit(0);
        }
    }
    for (int i = 0; i < Particles.size(); i++){           // insert all particles into tree
        Insert(Particles[i]);
    }
    Monople(root);                                        //initialize all monople of nodes in whole tree

}
QuadrupleTree::QuadrupleTree(double theta,std::vector<Particle> &Particles,double mX,double mY,double mZ,double MX,double MY, double MZ){  
    PtcVectorPtr = &Particles;
    THETA = theta;
    root = new TreeNode;                                  // allocate memory for root
    root->parent = nullptr;
    root->level =0;
    maxX = MX; maxY = MY; maxZ = MZ;                      // inintialize boundary of this tree
    minX = mX; minY = mY; minZ = mZ;
    for (int i = 0; i < Particles.size(); i++){           // insert all particles into tree
        if(Particles[i].posi[0]>MX or Particles[i].posi[0]<mX or Particles[i].posi[1]>MY or Particles[i].posi[1]<mY){
            std::cout<<"particle out of boundry";
            std::exit(0);
        }
    }
    for (int i = 0; i < Particles.size(); i++){           // insert all particles into tree
        Insert(Particles[i]);
    }
    Monople(root);                                        //initialize all monople of nodes in whole tree

}
QuadrupleTree::~QuadrupleTree(){
    // if (root != NULL){
    //     if(root->NE != NULL){DeleteNode(root->NE); root->NE = NULL;};
    //     if(root->NW != NULL){DeleteNode(root->NW); root->NW = NULL;};
    //     if(root->SE != NULL){DeleteNode(root->SE); root->SE = NULL;};
    //     if(root->SW != NULL){DeleteNode(root->SW); root->SW = NULL;};
    //     // delete root;
    // }
    // if(root->NE != NULL){delete root->NE ;root->NE = NULL;};
    // if(root->NW != NULL){delete root->NW ;root->NE = NULL;};
    // if(root->SE != NULL){delete root->SE ;root->NE = NULL;};
    // if(root->SW != NULL){delete root->SW ;root->NE = NULL;};
    if (root!=nullptr)
    {
        delete root;
    }
    
    
    root=nullptr;
    PtcVectorPtr = NULL;
}
void TreeNode::PrintNode(){                                                                                    // print the information of node parent, self and children are the memeory locations of each node
        std::cout << "parent: " <<  parent <<  "\n";
        std::cout << "self: " <<  this <<  "\n";
        std::cout << "children: " <<"NW:"<< NW <<" NE:"<< NE <<" SW:"<< SW <<" SE:"<< SE << "\n";
        std::cout << "level: " << level << "\n";
        if (leaf == true){                                                                               // If the node is a particle, print the parameter of it
            std::cout << "this node is a Particle\n";
            std::cout << "position: (" << ptclPtr->posi[0] << "," << ptclPtr->posi[1] << "," << ptclPtr->posi[2] << ")\n";
            std::cout << "velocity: (" << ptclPtr->velocity[0] << "," << ptclPtr->velocity[1] << ","<<ptclPtr->velocity[2] << ")\n";
            std::cout << "acceleration: (" << ptclPtr->acceleration[0] << "," << ptclPtr->acceleration[1] << "," << ptclPtr->acceleration[2] << ")\n";
            std::cout << "mass: " << ptclPtr->mass << "\n\n";
        }
        else{
            std::cout << "monople: (" << *monople<<","  << *(monople+1)<<","  << *(monople+2)<<","  << *(monople+3) <<")\n";
            std::cout << "\n\n";
            }
        };
TreeNode *QuadrupleTree::TwoParticleSubtree(TreeNode *ptcTree1, Particle &ptc2,double tempminX,double tempminY,double tempminZ,double tempmaxX,double tempmaxY, double tempmaxZ){
    double mX = tempminX;double mY = tempminY;double mZ = tempminZ;        // rename new boundary and cut it into four subregion
    double MX = tempmaxX;double MY = tempmaxY;double MZ = tempmaxZ;
    double midX = (mX + MX)/2.;double midY = (mY + MY)/2.;double midZ = (mZ + MZ)/2.;

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
    double tempMaxX{maxX};double tempMaxY{maxY};double tempMaxZ{maxZ}; // copy thw boundary of tree (as search gets deeper,
    double tempMinX{minX};double tempMinY{minY};double tempMinZ{minZ}; // the boundary shrinks for determining which node should be search)
    double tempMidX = (tempMaxX + tempMinX) / 2.;                    // mid-line of boundary (changes as boundary being changed)
    double tempMidY = (tempMaxY + tempMinY) / 2.;
    while (current) {
        if (newPtc.posi[0] < (tempMidX) and newPtc.posi[1] > (tempMidY)){ // the new particle is in northwest subregion
            if (current->NW != NULL ){                  // if current node's child node NW already exists
                if (current->NW->leaf == false){        // if this child node is not leaf (which means it's not a particle node and current node has grandchild) 
                    q.push(current->NW);                // push it to the tail of queue (search will continue, starting from this node)
                    current->leaf = false;              // mark that current node is not a particle
                }
                else{                                   // current node's child(NW) is a particle node(two particle are in the same region in this level)
                    TreeNode *tempNode = TwoParticleSubtree(current->NW, newPtc, tempMinX,tempMinY,tempMinZ,tempMaxX,tempMaxY,tempMaxZ)->NW; // creat a subtree contains both particles and connect it to the original node
                    // delete  current->NW->parent;
                    current->NW = tempNode;
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
                    TreeNode *tempNode = TwoParticleSubtree(current->NE, newPtc, tempMinX,tempMinY,tempMinZ,tempMaxX,tempMaxY,tempMaxZ)->NE; // creat a subtree contains both particles and connect it to the original node
                    // delete  current->NE->parent;
                    current->NE = tempNode;
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
                    TreeNode *tempNode = TwoParticleSubtree(current->SW, newPtc, tempMinX,tempMinY,tempMinZ,tempMaxX,tempMaxY,tempMaxZ)->SW; // creat a subtree contains both particles and connect it to the original node
                    // delete  current->SW->parent;
                    current->SW = tempNode;
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
                    TreeNode *tempNode = TwoParticleSubtree(current->SE, newPtc, tempMinX,tempMinY,tempMinZ,tempMaxX,tempMaxY,tempMaxZ)->SE; // creat a subtree contains both particles and connect it to the original node
                    // delete  current->SE->parent;
                    current->SE = tempNode;
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
        int parentChildNum = 0;
        if(Node->parent->NE != NULL){parentChildNum +=1; };
        if(Node->parent->NW != NULL){parentChildNum +=1; };
        if(Node->parent->SE != NULL){parentChildNum +=1; };
        if(Node->parent->SW != NULL){parentChildNum +=1; };
        if(parentChildNum == 0){Node->parent->leaf = true;}
    }
    if (Node->leaf){            // if the node is a particle node, disconnect to particle and delete self
       
        // delete Node->ptclPtr;
        Node->ptclPtr = NULL;
        delete Node->monople;
        Node->monople = NULL;
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
double *QuadrupleTree::Monople(TreeNode *Node){
    if (Node->monople != NULL)
    {   
        return Node->monople;
    }
    double *monoParaPtr = new double;     // {mass, x, y, z}
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
            double *monoNEptr = Monople(Node->NE);
            *monoParaPtr += *monoNEptr;
            *(monoParaPtr + 1) += *(monoNEptr + 1) * (*monoNEptr); // add center of mass times mass  
            *(monoParaPtr + 2) += *(monoNEptr + 2) * (*monoNEptr); // add center of mass times mass  
            // *(monoParaPtr + 3) += *(monoNEptr + 3) * (*monoNEptr); // add center of mass times mass  
            // delete[] monoNEptr;
        }
        if(Node->NW != NULL){
            double *monoNWptr = Monople(Node->NW);
            *monoParaPtr += *monoNWptr;
            *(monoParaPtr + 1) += *(monoNWptr + 1) * (*monoNWptr); // add center of mass times mass  
            *(monoParaPtr + 2) += *(monoNWptr + 2) * (*monoNWptr); // add center of mass times mass  
            // *(monoParaPtr + 3) += *(monoNWptr + 3) * (*monoNWptr); // add center of mass times mass  
            // delete[] monoNWptr;
        }
        if(Node->SE != NULL){
            double *monoSEptr = Monople(Node->SE);
            *monoParaPtr += *monoSEptr;
            *(monoParaPtr + 1) += *(monoSEptr + 1) * (*monoSEptr); // add center of mass times mass  
            *(monoParaPtr + 2) += *(monoSEptr + 2) * (*monoSEptr); // add center of mass times mass  
            // *(monoParaPtr + 3) += *(monoSEptr + 3) * (*monoSEptr); // add center of mass times mass  
            // delete[] monoSEptr;
        }
        if(Node->SW != NULL){
            double *monoSWptr = Monople(Node->SW);
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
std::vector<double> TreeNode::CalculateForce(TreeNode *Node, Particle &tarPtc){

    double r{0};double a{0}; 
    std::vector<double> acc{0., 0.};
    // r is the distance between the Node and the particle
    r = sqrt((tarPtc.posi[0] - *(Node->monople +1))*(tarPtc.posi[0] - *(Node->monople +1))+(tarPtc.posi[1] - *(Node->monople +2))*(tarPtc.posi[1] - *(Node->monople +2)));
    // a is the acceleration of the particle using the Newton's law
    // a is divided by one more r because the next step is divided by one less r
    a = G_CONST * *(Node->monople) / (r * r * r);
    acc[0] = a * (*(Node->monople +1) - tarPtc.posi[0]); // the acceleration of x component ( a * x/r )
    acc[1] = a * (*(Node->monople +2) - tarPtc.posi[1]); // the acceleration of y component ( a * y/r )

    return acc; // include x and y axis
}

//calculate the total force acting on the given particle
void QuadrupleTree::TotalForce(Particle &Ptc){
    std::queue<TreeNode*> q;
    TreeNode *current = root;
    std::vector<double> accSum {0.0, 0.0}; // the sum of acceleration, including x and y components
    std::vector<TreeNode*>section; // change NE, NW, SE, SW into a vector to use for-loop
    std::vector<double> acc; // to store the computing acceleration
    double r{0};double d{0};
    while (current){
        if (current->leaf==true and current->ptclPtr == &Ptc){
            if (q.empty() == false){
            current = q.front();
            q.pop();
            }
            else{break;}
            continue;
        }
        r = sqrt((Ptc.posi[0] - *(current->monople +1))*(Ptc.posi[0] - *(current->monople +1))+(Ptc.posi[1] - *(current->monople +2))*(Ptc.posi[1] - *(current->monople +2)));
        d = (maxX  - minX) / pow(2.0, (current->level) * 1.0);
        // this is the node where the particle itself exists
        if (r <= 0){
            break; 
        }
        // test if Multipole-Acceptance-Criterion(MAC) can be used in this node
        if (current->leaf == true){
            acc = current->CalculateForce(current, Ptc); // directly calculate the force(acceleration)
            accSum[0] += acc[0]; // x component
            accSum[1] += acc[1]; // y component
        }
        else{
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
            }
        }
        if (q.empty() == false){
            current = q.front();
            q.pop();
        }
        else{break;}
    }
    Ptc.acceleration[0] = accSum[0];
    Ptc.acceleration[1] = accSum[1];
    return;
}
void QuadrupleTree::TreeForce(){
    for (int i = 0; i < PtcVectorPtr->size(); i++){           
        TotalForce((*PtcVectorPtr)[i]);
    }
};



