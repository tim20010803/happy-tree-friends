// C++ code
#include <iostream>
#include <sstream>
#include <queue>
#include <vector>
#include <string>
struct Particle
{
    std::vector<float> posi;
    std::vector<float> velocity;
    std::vector<float> acceleration;
    float mass;
};

class QuadrupleTree;
class TreeNode{
//private:
public://for testing
    TreeNode *parent;
    TreeNode *NW;//initially set it point to nullptr and use it to determin if this node is leaf
    TreeNode *NE;
    TreeNode *SW;
    TreeNode *SE;
    Particle *monoplePtr;//store the total mass and center of mass in this subtree
    int level;//zero for root
    bool leaf;//if leaf == true then monoplePtr are the leaf of Particle

    TreeNode():NW(0),NE(0),SW(0),SE(0),parent(0),level(0),leaf(false),monoplePtr{NULL}{};
    TreeNode(Particle *newPtcl):NW(0),NE(0),SW(0),SE(0),parent(0),level(0),leaf(true),monoplePtr(newPtcl){};
    void printNode(){
        std::cout << "parent: " <<  parent <<  "\n";
        std::cout << "children: " <<"NW:"<< NW <<" NE:"<< NE <<" SW:"<< SW <<" SE:"<< SE << "\n";
        std::cout << "level: " << level << "\n";
        if (leaf == true)
        {
            std::cout << "this node is a Particle\n";
            std::cout << "position: (" << monoplePtr->posi[0] << "," << monoplePtr->posi[1] << "," << monoplePtr->posi[2] << ")\n";
            std::cout << "velocity: (" << monoplePtr->velocity[0] << "," << monoplePtr->velocity[1] << ","<<monoplePtr->velocity[2] << ")\n";
            std::cout << "mass: " << monoplePtr->mass << "\n\n";
        }
        else{std::cout << "\n\n";}
        
        
        };
        friend class QuadrupleTree;
};

class QuadrupleTree{
private:
    
    float maxX{1};
    float maxY{1};
    float maxZ{1};
    float minX{0};
    float minY{0};
    float minZ{0};
public:
    TreeNode *root;
    QuadrupleTree():root(0){};
    QuadrupleTree(Particle &firstPtc,float mX,float mY,float mZ,float MX,float MY, float MZ);//mX is minX, MX is maxX(boundry);

    void Construct(std::vector<Particle> Particles);
    void Insert(Particle& newPtc);
    
    TreeNode* leftmost(TreeNode *current);
    TreeNode* InorderSuccessor(TreeNode *current);
    void Inorder_by_parent();
};
QuadrupleTree::QuadrupleTree(Particle &firstPtc,float mX,float mY,float mZ,float MX,float MY, float MZ){//mX is minX, MX is maxX(boundry);
    root = new TreeNode;           // allocate memory for root
    maxX = MX; maxY = MY; maxZ = MZ; 
    minX = mX; minY = mY; minZ = mZ; 
    Insert(firstPtc);

}

// C++ code
// TreeNode* QuadrupleTree::leftmost(TreeNode *current){
//     while (current->NW != NULL){
//         current = current->NW;
//     }
//     return current;
// }
// void QuadrupleTree::Inorder_by_parent(){
//     TreeNode *current = new TreeNode;
//     current = leftmost(root);

//     while(current){
//         std::cout << current->monoplePtr << " ";
//         current = InorderSuccessor(current);
//     }
// }
// TreeNode* QuadrupleTree::InorderSuccessor(TreeNode *current){
//     if (current->NE != NULL){
//         return leftmost(current->NE);
//     }

//     // 利用兩個pointer: successor與current做traversal 

//     TreeNode *successor = current->parent;   
//     while (successor != NULL && current == successor->NE) {
//         current = successor;
//         successor = successor->parent;
//     }
//     return successor;
// }

// void QuadrupleTree::Construct(std::vector<Particle> Particles){
//     std::queue<TreeNode*> q;         // create a queue to handle level-roder rule
//     TreeNode *current = root;        // point *current to root
//     while (ss >> monoplePtr) {
//         if (monoplePtr >= 65 && monoplePtr <= 90){                // 處理current的NW
//             TreeNode *new_node = new TreeNode(monoplePtr);  // call constructor TreeNode(char s)
//             new_node->parent = current;
//             current->NW = new_node;
//             q.push(new_node);
//         }
//         if (!(ss >> monoplePtr)){                           // 有可能char array含有奇數筆資料
//             break;                                    // 所以在這裡結束迴圈
//         }
//         if (monoplePtr >= 65 && monoplePtr <= 90){                // 處理current的NE
//             TreeNode *new_node = new TreeNode;        // call constructor TreeNode()
//             new_node->parent = current;
//             current->NE = new_node;
//             new_node->monoplePtr = monoplePtr;                    // assign monoplePtr to new_node
//             q.push(new_node);
//         }
//         current = q.front();                          // 從queue中更新current
//         q.pop();                                      // 更新queue
//     }
// }
// TreeNode* TwoParticleTree(TreeNode* ptcT1,Particle ptcT2){
    
// }

TreeNode* TwoParticleSubtree(TreeNode* ptcTree1, Particle& ptc2,float tempminX,float tempminY,float tempminZ,float tempmaxX,float tempmaxY, float tempmaxZ){
    float mX = tempminX;float mY = tempminY;float mZ = tempminZ;//rename new boundary and cut it into four subregion
    float MX = tempmaxX;float MY = tempmaxY;float MZ = tempmaxZ;
    float midX = (mX + MX)/2.;float midY = (mY + MY)/2.;float midZ = (mZ + MZ)/2.;

    bool firstPtlN = ptcTree1->monoplePtr->posi[1] > midY;//determin each Particle in which part of subregion
    bool firstPtlE = ptcTree1->monoplePtr->posi[0] > midX;
    bool secPtlN = ptc2.posi[1] > midY;
    bool secPtlE = ptc2.posi[0] > midX;

    TreeNode* output = new TreeNode;
    output->level = ptcTree1->level - 1;
    
    ptcTree1->level += 1;
    TreeNode* tempOut = output;
    while (true)
    {
        if ( (firstPtlE == secPtlE) and (firstPtlN == secPtlN)){
            // std::cout <<"in the same category\n"; //print for debug
            tempOut->leaf = false;
            tempOut->monoplePtr = NULL;

            if (firstPtlE and firstPtlN){
                mX = midX;
                mY = midY;
                tempOut->NE = new TreeNode;
                tempOut->NE->parent = tempOut;
                tempOut->NE->level = tempOut->level + 1;
                tempOut = tempOut->NE;
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
            midX = (mX + MX)/2.; midY = (mY + MY)/2.; midZ = (mZ + MZ)/2.;
            firstPtlN = ptcTree1->monoplePtr->posi[1] > midY;//determin each Particle in which part of subregion
            firstPtlE = ptcTree1->monoplePtr->posi[0] > midX;
            secPtlN = ptc2.posi[1] > midY;
            secPtlE = ptc2.posi[0] > midX;

        }
        else{//two Particle are not in the same category
            //connect ptcTree1 to tempOut
            // std::cout <<"not in the same category\n";//print for debug
            tempOut->leaf = false;
            tempOut->monoplePtr = NULL;
            if (firstPtlE and firstPtlN){
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
            }
            //creat node for second paticle and connect it to tempOut
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

    std::queue<TreeNode*> q;
    TreeNode *current = root;
    float tempMaxX{maxX};float tempMaxY{maxY};float tempMaxZ{maxZ};
    float tempMinX{minX};float tempMinY{minY};float tempMinZ{minZ};

    while (current) {
        if (newPtc.posi[0] < ((tempMaxX + tempMinX) / 2.) and newPtc.posi[1] > ((tempMaxY + tempMinY) / 2.))
        {
            if (current->NW != NULL ){               // current的NW有分支
                if (current->NW->leaf == false)      //該分支不是particle
                {
                    q.push(current->NW);                // 將其推進queue中
                    current->leaf = false;
                }
                else{                                   //current的NW已經有particle,做一棵subtree接上去
                    current->NW = TwoParticleSubtree(current->NW, newPtc, tempMinX,tempMinY,tempMinZ,tempMaxX,tempMaxY,tempMaxZ)->NW;
                    delete [] current->NW->parent;
                    current->NW->parent = current;
                    break;
                }
                
            }
            else{         // current的NW沒有分支
                TreeNode *new_node = new TreeNode(&newPtc);   // 建立新的node, 將Particle放在這裡
                new_node->parent = current;
                new_node->level = current->level + 1;
                new_node->leaf = true;
                current->leaf = false;
                current->NW = new_node;
                break;                         
            }
            tempMaxX = ((tempMaxX + tempMinX) / 2.);
            tempMinY = ((tempMaxY + tempMinY) / 2.);
        }
        else if (newPtc.posi[0] >= ((tempMaxX + tempMinX) / 2.) and newPtc.posi[1] > ((tempMaxY + tempMinY) / 2.))
        {
            if (current->NE != NULL ){               // current的NE有分支
                if (current->NE->leaf == false)      //該分支不是particle
                {
                    q.push(current->NE);                // 將其推進queue中
                    current->leaf = false;
                }
                else{                                   //current的NE已經有particle,做一棵subtree接上去
                    current->NE = TwoParticleSubtree(current->NE, newPtc, tempMinX,tempMinY,tempMinZ,tempMaxX,tempMaxY,tempMaxZ)->NE;
                    delete [] current->NE->parent;
                    current->NE->parent = current;
                    break;
                }
            }
            else{                                          // current的NE沒有分支
                TreeNode *new_node = new TreeNode(&newPtc);   // 建立新的node, 將Particle放在這裡
                new_node->parent = current;
                new_node->level = current->level + 1;
                new_node->leaf = true;
                current->leaf = false;
                current->NE = new_node;
                break;
            }
            tempMinX = ((tempMaxX + tempMinX) / 2.);
            tempMinY = ((tempMaxY + tempMinY) / 2.);
        }
        else if (newPtc.posi[0] < ((tempMaxX + tempMinX) / 2.) and newPtc.posi[1] <= ((tempMaxY + tempMinY) / 2.))
        {
            if (current->SW != NULL ){               // current的SW有分支
                if (current->SW->leaf == false)      //該分支不是particle
                {
                    q.push(current->SW);                // 將其推進queue中
                    current->leaf = false;
                }
                else{                                   //current的SW已經有particle,做一棵subtree接上去
                    current->SW = TwoParticleSubtree(current->SW, newPtc, tempMinX,tempMinY,tempMinZ,tempMaxX,tempMaxY,tempMaxZ)->SW;
                    delete [] current->SW->parent;
                    current->SW->parent = current;
                    break;
                }
            }
            else{                                          // current的SW沒有分支
                TreeNode *new_node = new TreeNode(&newPtc);   // 建立新的node, 將Particle放在這裡
                new_node->parent = current;
                new_node->level = current->level + 1;
                new_node->leaf = true;
                current->leaf = false;
                current->SW = new_node;
                break;                
            }
            tempMaxX = ((tempMaxX + tempMinX) / 2.);
            tempMaxY = ((tempMaxY + tempMinY) / 2.);
        }
        else if (newPtc.posi[0] >= ((tempMaxX + tempMinX) / 2.) and newPtc.posi[1] <= ((tempMaxY + tempMinY) / 2.))
        {
            if (current->SE != NULL ){               // current的SE有分支且該分支不是particle
                if (current->SE->leaf == false)      //該分支不是particle
                {
                    q.push(current->SE);                // 將其推進queue中
                    current->leaf = false;
                }
                else{                                   //current的SE已經有particle,做一棵subtree接上去
                    current->SE = TwoParticleSubtree(current->SE, newPtc, tempMinX,tempMinY,tempMinZ,tempMaxX,tempMaxY,tempMaxZ)->SE;
                    delete [] current->SE->parent;
                    current->SE->parent = current;
                    break;
                }
            }
            else{                                          // current的SE沒有分支
                TreeNode *new_node = new TreeNode(&newPtc);   // 建立新的node, 將Particle放在這裡
                new_node->parent = current;
                new_node->level = current->level + 1;
                new_node->leaf = true;
                current->leaf = false;
                current->SE = new_node;
                break;                         
            }
            tempMinX = ((tempMaxX + tempMinX) / 2.);
            tempMaxY = ((tempMaxY + tempMinY) / 2.);
        }
        current = q.front();
        q.pop();
    }
}








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
    // AnodePtr->printNode();



    // TreeNode* BinA = TwoParticleSubtree(AnodePtr,b,0.,0.,0.,10.,10.,10.);
    // std::cout<<"\n\n" ;
    // if (BinA != NULL)
    // {
    // BinA->printNode();
    //     if (BinA->SE != NULL)
    //     {
    //         std::cout <<"SE:\n";
    //         BinA->SE->printNode();
    //     }
    //     if (BinA->NE != NULL)
    //     {
    //         std::cout <<"NE:\n";
    //         BinA->NE->printNode();
    //     }
    //     if (BinA->SW != NULL)
    //     {
    //         std::cout <<"SW:\n";
    //         BinA->SW->printNode();
    //     }
    //     if (BinA->NW != NULL)
    //     {
    //         std::cout <<"NW:\n";
    //         BinA->NW->printNode();
    //         BinA->NW->SW->printNode();
    //         BinA->NW->SW->NE->printNode();
    //         BinA->NW->SW->SW->printNode();
    //     }

    // }
    
    
    QuadrupleTree T(a,0.,0.,0.,10.,10.,10.); 


    T.Insert(b);
    T.Insert(c);
    T.Insert(d);
    T.root->printNode();
    T.root->SW->printNode();
    T.root->NW->printNode();
    T.root->NW->SE->printNode();
    T.root->NW->SW->printNode();
    T.root->NW->SW->SW->printNode();
    T.root->NW->SW->NE->printNode();
    // T.root->NE->printNode();
    // T.root->SW->printNode();
    
    // T.root->printNode();
    // T.root->NE->NE->printNode();

    //for the first test I found level update error and deletion of old Particle
    return 0;
}
