

#ifndef THETA
#define RESOLUTION 1e-10;
#define THETA 1.0
#endif
class QuadrupleTree;
class TreeNode{
// private:
public:                                // set all of parameter be public for testing. you should only access all those parameters by function in QuadrupleTree.
    TreeNode *parent{0};                  // point to the parent node
    TreeNode *NW{0};                      // northwest child node
    TreeNode *NE{0};                      // northeast child node
    TreeNode *SW{0};                      // southwest child node
    TreeNode *SE{0};                      // southeast child node
    Particle *ptclPtr{0};                 // the pointer of single particle if this node is a leaf(without any child node), otherwise it's a null pointer.  
    int level;                         // the depth level of this node (which is zero for root node)
    bool leaf;                         // if leaf is true then this node is a leaf the point toward the particle is ptclPtr
    double *monople{0};                    // storege the monople of all subtree
    std::vector<double> CalculateForce(TreeNode *Node, Particle &tarPtc); //Compute the force(acceleration) acting from this node to a particle p

    ~TreeNode();
    TreeNode():NW(NULL),NE(NULL),SW(NULL),SE(NULL),parent(NULL),level(0),leaf(false),ptclPtr{NULL},monople{NULL}{};                   // constuctor of TreeNode (which creats a TreeNode object and initialize parameters) 
    TreeNode(Particle *newPtcl):NW(NULL),NE(NULL),SW(NULL),SE(NULL),parent(NULL),level(0),leaf(true),ptclPtr(newPtcl),monople{NULL}{};// constuctor of TreeNode which stores the particle's location into pointer and set leaf==true (since this node is a particle)
    void PrintNode();
        friend class QuadrupleTree;   // give the QuadrupleTree class access to those private members(such as private functions and parameters)
};



class QuadrupleTree{
private:
    double minX{0};double minY{0};double minZ{0}; // the boundary of space (minimum and maximum of x,y,z of all particles)
    double maxX{1};double maxY{1};double maxZ{1};
    std::vector<Particle>* PtcVectorPtr{0};
    TreeNode *TwoParticleSubtree(TreeNode *ptcTree1, Particle &ptc2,double tempminX,double tempminY,double tempminZ,double tempmaxX,double tempmaxY, double tempmaxZ); 
    // TwoParticleSubtree function creats a tree which connects two node particle,
    // if a particle in the same region of existing particle in a node 
    // in the original level of first particle
public:
    TreeNode *root;
    QuadrupleTree():root(NULL){};
    QuadrupleTree(Particle &firstPtc,double mX,double mY,double mZ,double MX,double MY, double MZ);                // take one particle and boundary of all particle to ininitaize the tree, and mX is minX(minimum x), MX is maxX(maxmum x);
    QuadrupleTree(std::vector<Particle> &Particles,double mX,double mY,double mZ,double MX,double MY, double MZ);  // take several particles(with type of std::vector) and boundary of all particle to ininitaize the tree, and mX is minX(minimum x), MX is maxX(maxmum x);
    ~QuadrupleTree();                  // desturctor (to destroy the whole tree and release the memory space it takes)
    void DeleteNode(TreeNode *Node);                      // delete the Node and its all descendent
    void Insert(Particle& newPtc);                        // insert one particle in the tree
    void Trim(TreeNode *Node);                            // delete all empty subnodes and itself if the input node turn out to be a empty (those node without any child and not a particle node)
    double *Monople(TreeNode *Node);                       // return total mass and center of mass of this subtree: (mass, x, y, z) and initialize the monople at each tree (it don't modify any monople in the subtree, it just read it.)
    void TotalForce(Particle &Ptc);                       // calculate the total force(acceleration) of a given particle
    void TreeForce();
};
