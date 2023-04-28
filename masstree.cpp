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
    TreeNode(particle newPtcl):NW(0),NE(0),SW(0),SE(0),parent(0),level(0),leaf(false),monople(newPtcl){};
    void printNode(){
        std::cout << "parent: " <<  parent <<  "\n";
        std::cout << "children: " << NW << NE << SW << SE << "\n";
        std::cout << "If this is a leaf: " << leaf << "\n";
        std::cout << "level: " << level << "\n";
        std::cout << "monople position: (" << monople.posi[0] << "," << monople.posi[1] << "," << monople.posi[2] << ")\n";
        std::cout << "monople velocity: (" << monople.velocity[0] << "," << monople.velocity[1] << ","<<monople.velocity[2] << ")\n";
        std::cout << "monople mass: " << monople.mass << "\n";
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
    
    TreeNode* TwoParticleSubtree(TreeNode* ptcTree1, TreeNode* ptcTree2);//use this function to constrcut classify subtree
    
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
                current->monople.velocity = {0,0,0};
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
                current->monople.velocity = {0,0,0};
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
                current->monople.velocity = {0,0,0};
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
                current->monople.velocity = {0,0,0};
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
    a.posi = {1,6,4};
    a.velocity = {4,3,2};
    a.mass = {12};
    particle b;
    b.posi = {2,7,8};
    b.velocity = {1,6,7};
    b.mass = {23};
    particle c;
    c.posi = {3,8,8};
    c.velocity = {1,6,7};
    c.mass = {212};

    std::vector<float>minposition = {0,0,0};
    std::vector<float>maxposition = {5,5,5};
    QuadrupleTree T(a,0,0,0,10,10,10);  
    T.Insert(b);
    T.Insert(c);
    T.root->printNode();
    T.root->NE->printNode();
    T.root->NE->NE->printNode();
    //for the first test I found level update error and deletion of old Particle
    return 0;
}
