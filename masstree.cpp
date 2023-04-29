// C++ code
#include <iostream>
#include <sstream>
#include <queue>
#include <vector>
#include <string>
struct particle
{
    std::vector<float> posi;
    std::vector<float> velocity;
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
    particle monople;//store the total mass and center of mass in this subtree
    int level;//zero for root
    bool leaf;//if leaf == true then monople are the leaf of Particle

    TreeNode():NW(0),NE(0),SW(0),SE(0),parent(0),level(0),leaf(false),monople{{0.,0.},{0.,0.},0.}{};
    TreeNode(particle newPtcl):NW(0),NE(0),SW(0),SE(0),parent(0),level(0),leaf(true),monople(newPtcl){};
    void printNode(){
        std::cout << "parent: " <<  parent <<  "\n";
        std::cout << "children: " <<"NW:"<< NW <<"NE:"<< NE <<"SW:"<< SW <<"SE:"<< SE << "\n";
        std::cout << "If this is a leaf: " << leaf << "\n";
        std::cout << "level: " << level << "\n";
        std::cout << "monople position: (" << monople.posi[0] << "," << monople.posi[1] << "," << monople.posi[2] << ")\n";
        std::cout << "monople velocity: (" << monople.velocity[0] << "," << monople.velocity[1] << ","<<monople.velocity[2] << ")\n";
        std::cout << "monople mass: " << monople.mass << "\n\n";
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
    QuadrupleTree(particle firstPtc,float mX,float mY,float mZ,float MX,float MY, float MZ);//mX is minX, MX is maxX(boundry);

    void Construct(std::vector<particle> Particles);
    void Insert(particle newPtc);
    
    TreeNode* leftmost(TreeNode *current);
    TreeNode* InorderSuccessor(TreeNode *current);
    void Inorder_by_parent();
};
QuadrupleTree::QuadrupleTree(particle firstPtc,float mX,float mY,float mZ,float MX,float MY, float MZ){//mX is minX, MX is maxX(boundry);
    root = new TreeNode;           // allocate memory for root
    root->monople = firstPtc;      // assign character to root
    root->leaf = true;             // label this node as leaf (without subnode in it)
    float maxX{MX};float maxY{MY};float maxZ{MZ};
    float minX{mX};float minY{mY};float minZ{mZ};
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
//         std::cout << current->monople << " ";
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
//     while (ss >> monople) {
//         if (monople >= 65 && monople <= 90){                // 處理current的NW
//             TreeNode *new_node = new TreeNode(monople);  // call constructor TreeNode(char s)
//             new_node->parent = current;
//             current->NW = new_node;
//             q.push(new_node);
//         }
//         if (!(ss >> monople)){                           // 有可能char array含有奇數筆資料
//             break;                                    // 所以在這裡結束迴圈
//         }
//         if (monople >= 65 && monople <= 90){                // 處理current的NE
//             TreeNode *new_node = new TreeNode;        // call constructor TreeNode()
//             new_node->parent = current;
//             current->NE = new_node;
//             new_node->monople = monople;                    // assign monople to new_node
//             q.push(new_node);
//         }
//         current = q.front();                          // 從queue中更新current
//         q.pop();                                      // 更新queue
//     }
// }
// TreeNode* TwoParticleTree(TreeNode* ptcT1,Particle ptcT2){
    
// }

TreeNode* TwoParticleSubtree(TreeNode* ptcTree1, particle ptc2,float tempminX,float tempminY,float tempminZ,float tempmaxX,float tempmaxY, float tempmaxZ){
    float mX = tempminX;float mY = tempminY;float mZ = tempminZ;//rename new boundary and cut it into four subregion
    float MX = tempmaxX;float MY = tempmaxY;float MZ = tempmaxZ;
    float midX = (mX + MX)/2.;float midY = (mY + MY)/2.;float midZ = (mZ + MZ)/2.;

    bool firstPtlN = ptcTree1->monople.posi[1] > midY;//determin each particle in which part of subregion
    bool firstPtlE = ptcTree1->monople.posi[0] > midX;
    bool secPtlN = ptc2.posi[1] > midY;
    bool secPtlE = ptc2.posi[0] > midX;

    TreeNode* output = new TreeNode;
    output->monople.mass = ptcTree1->monople.mass + ptc2.mass;
    output->level = ptcTree1->level;
    output->monople.velocity={0.,0.,0.};
    for (int i = 0; i < 3; i++)
    {
        output->monople.posi[i] = ptcTree1->monople.posi[i] * ptcTree1->monople.mass + ptc2.posi[i] * ptc2.mass;
        output->monople.posi[i] /= output->monople.mass;
    }
    ptcTree1->level += 1;
    TreeNode* tempOut = output;
    while (true)
    {
        if ( (firstPtlE == secPtlE) and (firstPtlN == secPtlN)){
            std::cout <<"in the same category\n";
            tempOut->printNode();
            tempOut->leaf = false;

            if (firstPtlE and firstPtlN){
                mX = midX;
                mY = midY;
                tempOut->NE = new TreeNode(tempOut->monople);
                tempOut->NE->parent = tempOut;
                tempOut->NE->level = tempOut->level + 1;
                tempOut = tempOut->NE;
            }
            else if (firstPtlE == false and firstPtlN){
                MX = midX;
                mY = midY;
                tempOut->NW = new TreeNode(tempOut->monople);
                tempOut->NW->parent = tempOut;
                tempOut->NW->level = tempOut->level + 1;
                tempOut = tempOut->NW;
            }
            else if (firstPtlE and firstPtlN  == false){
                mX = midX;
                MY = midY;
                tempOut->SE = new TreeNode(tempOut->monople);
                tempOut->SE->parent = tempOut;
                tempOut->SE->level = tempOut->level + 1;
                tempOut = tempOut->SE;
            }
            else if (firstPtlE == false and firstPtlN == false){
                MX = midX;
                MY = midY;
                tempOut->SW = new TreeNode(tempOut->monople);
                tempOut->SW->parent = tempOut;
                tempOut->SW->level = tempOut->level + 1;
                tempOut = tempOut->SW;
            }
            midX = (mX + MX)/2.; midY = (mY + MY)/2.; midZ = (mZ + MZ)/2.;
            firstPtlN = ptcTree1->monople.posi[1] > midY;//determin each particle in which part of subregion
            firstPtlE = ptcTree1->monople.posi[0] > midX;
            secPtlN = ptc2.posi[1] > midY;
            secPtlE = ptc2.posi[0] > midX;

        }
        else{//two particle are not in the same category
            //connect ptcTree1 to tempOut
            std::cout <<"not in the same category\n";
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
                tempOut->NE = new TreeNode(ptc2);
                tempOut->NE->parent = tempOut;
                tempOut->NE->level = tempOut->level + 1;
                tempOut->NE->leaf = true;
            }
            else if (secPtlE == false and secPtlN){
                tempOut->NW = new TreeNode(ptc2);
                tempOut->NW->parent = tempOut;
                tempOut->NW->level = tempOut->level + 1;
                tempOut->NW->leaf = true;
            }
            else if (secPtlE and secPtlN  == false){
                tempOut->SE = new TreeNode(ptc2);
                tempOut->SE->parent = tempOut;
                tempOut->SE->level = tempOut->level + 1;
                tempOut->SE->leaf = true;

            }
            else if (secPtlE == false and secPtlN == false){
                tempOut->SW = new TreeNode(ptc2);
                tempOut->SW->parent = tempOut;
                tempOut->SW->level = tempOut->level + 1;
                tempOut->SW->leaf = true;
            }
            tempOut->printNode();
            break;
        }
    }
    
    
    return output;
    
}

void QuadrupleTree::Insert(particle newPtc){    

    std::queue<TreeNode*> q;
    TreeNode *current = root;
    float tempMaxX{maxX};float tempMaxY{maxY};float tempMaxZ{maxZ};
    float tempMinX{minX};float tempMinY{minY};float tempMinZ{minZ};

    while (current) {
        if (newPtc.posi[0] < ((tempMaxX + tempMinX) / 2.) and newPtc.posi[1] > ((tempMaxY + tempMinY) / 2.))
        {
            if (current->NW != NULL){               // current的NW有分支
                q.push(current->NW);                // 將其推進queue中
                current->leaf = false;
            }
            else{         // current的NW沒有分支,且current分支不是leaf(不是particle)
                TreeNode *new_node = new TreeNode(newPtc);   // 建立新的node, 將particle放在這裡
                new_node->parent = current;
                new_node->level = current->level + 1;
                new_node->leaf = true;
                current->leaf = false;
                current->NW = new_node;
                
                for (int i = 0; i < 3; i++)
                {
                    current->monople.posi[i] = newPtc.mass * newPtc.posi[i] + current->monople.mass * current->monople.posi[i];
                    current->monople.posi[i] /= (current->monople.mass + newPtc.mass);
                }
                current->monople.mass += newPtc.mass;
                current->monople.velocity = {0.,0.,0.};
                break;                         
            }
            tempMaxX = ((tempMaxX + tempMinX) / 2.);
            tempMinY = ((tempMaxY + tempMinY) / 2.);
        }
        else if (newPtc.posi[0] >= ((tempMaxX + tempMinX) / 2.) and newPtc.posi[1] > ((tempMaxY + tempMinY) / 2.))
        {
            if (current->NE != NULL){               // current的NE有分支
                q.push(current->NE);                // 將其推進queue中
                current->leaf = false;
            }
            else{                                          // current的NE沒有分支
                TreeNode *new_node = new TreeNode(newPtc);   // 建立新的node, 將particle放在這裡
                new_node->parent = current;
                new_node->level = current->level +1;
                current->leaf = false;
                current->NE = new_node;
                
                for (int i = 0; i < 3; i++)
                {
                    current->monople.posi[i] = newPtc.mass * newPtc.posi[i] + current->monople.mass * current->monople.posi[i];
                    current->monople.posi[i] /= (current->monople.mass + newPtc.mass);
                }
                current->monople.mass += newPtc.mass;
                current->monople.velocity = {0.,0.,0.};
                break;                         
            }
            tempMinX = ((tempMaxX + tempMinX) / 2.);
            tempMinY = ((tempMaxY + tempMinY) / 2.);
        }
        else if (newPtc.posi[0] < ((tempMaxX + tempMinX) / 2.) and newPtc.posi[1] <= ((tempMaxY + tempMinY) / 2.))
        {
            if (current->SW != NULL){               // current的SW有分支
                q.push(current->SW);                // 將其推進queue中
                current->leaf = false;
            }
            else{                                          // current的SW沒有分支
                TreeNode *new_node = new TreeNode(newPtc);   // 建立新的node, 將particle放在這裡
                new_node->parent = current;
                new_node->level = current->level + 1;
                current->leaf = false;
                current->SW = new_node;
                
                for (int i = 0; i < 3; i++)
                {
                    current->monople.posi[i] = newPtc.mass * newPtc.posi[i] + current->monople.mass * current->monople.posi[i];
                    current->monople.posi[i] /= (current->monople.mass + newPtc.mass);
                }
                current->monople.mass += newPtc.mass;
                current->monople.velocity = {0.,0.,0.};
                break;                         
            }
            tempMaxX = ((tempMaxX + tempMinX) / 2.);
            tempMaxY = ((tempMaxY + tempMinY) / 2.);
        }
        else if (newPtc.posi[0] >= ((tempMaxX + tempMinX) / 2.) and newPtc.posi[1] <= ((tempMaxY + tempMinY) / 2.))
        {
            if (current->SE != NULL){               // current的SE有分支
                q.push(current->SE);                // 將其推進queue中
                current->leaf = false;
            }
            else{                                          // current的SE沒有分支
                TreeNode *new_node = new TreeNode(newPtc);   // 建立新的node, 將particle放在這裡
                new_node->parent = current;
                new_node->level = current->level + 1;
                current->leaf = false;
                current->SE = new_node;
                
                for (int i = 0; i < 3; i++)
                {
                    current->monople.posi[i] = newPtc.mass * newPtc.posi[i] + current->monople.mass * current->monople.posi[i];
                    current->monople.posi[i] /= (current->monople.mass + newPtc.mass);
                }
                current->monople.mass += newPtc.mass;
                current->monople.velocity = {0.,0.,0.};
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
    particle a;
    a.posi = {1.,6.,4.};
    a.velocity = {4.,3.,2.};
    a.mass = {12.};
    particle b;
    b.posi = {2.,7.,8.};
    b.velocity = {1.,6.,7.};
    b.mass = {23.};
    particle c;
    c.posi = {3.,8.,8.};
    c.velocity = {1.,6.,7.};
    c.mass = {212.};
    TreeNode* AnodePtr = new TreeNode (a);
    TreeNode* BnodePtr = new TreeNode (b);
    TreeNode* CnodePtr = new TreeNode (c);
    // AnodePtr->printNode();



    TreeNode* BinA = TwoParticleSubtree(AnodePtr,b,0.,0.,0.,10.,10.,10.);
    std::cout<<"\n\n" ;
    if (BinA != NULL)
    {
    BinA->printNode();
        if (BinA->SE != NULL)
        {
            std::cout <<"SE:\n";
            BinA->SE->printNode();
        }
        if (BinA->NE != NULL)
        {
            std::cout <<"NE:\n";
            BinA->NE->printNode();
        }
        if (BinA->SW != NULL)
        {
            std::cout <<"SW:\n";
            BinA->SW->printNode();
        }
        if (BinA->NW != NULL)
        {
            std::cout <<"NW:\n";
            BinA->NW->printNode();
            BinA->NW->SW->printNode();
            BinA->NW->SW->NE->printNode();
            BinA->NW->SW->SW->printNode();
        }

    }
    
    
    // QuadrupleTree T(a,0.,0.,0.,10.,10.,10.); 
    // T.root->printNode();
    // T.Insert(b);
    // T.root->printNode();
    // T.root->NW->printNode();

    //for the first test I found level update error and deletion of old Particle
    return 0;
}
