#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>//I/O
#include <cmath>//abs(),pow(),ceil()
#include <float.h>//DBL_MAX
#include <algorithm>//sort()
#include <stack>
#include <queue>

using namespace std;

class Rnode{
public:
	vector<double> position;
	vector<double> size;
	vector<Rnode*> childs;
	int dimensions;
	bool isleaf;

	Rnode(int dimensions_, bool isleaf_ );
    void insertNode(Rnode row);
	void insertNodePtr(Rnode* row);
	bool intersects(Rnode Nrow);
	double interseccion(Rnode Nrow);
    double crecimiento(Rnode Nrow);
    double overlapcrecimiento(int index,Rnode Nrow);
    void intentacrecer(Rnode Nrow);
	double area();
    void actMBR();

};
class funcSortByDim{
public:
    int dim;
    int maxdim;
    int base;
    funcSortByDim();
    bool operator()(Rnode* vec1,Rnode* vec2);
    void up();
};
class Rtree{
public:
		Rnode *root;
		int M,m;
		int dimensions;
		int OTlastcallLvl;
        funcSortByDim funct;
		Rtree(int M_,int m_,int dimensions_);	
        void tryprint();
        Rnode* chooseSubTree(Rnode* actualnode, Rnode Nrow, int &level);
        Rnode * Split(Rnode * node);
        void reinsert(Rnode* N);
		void overflowTreatment(Rnode* actualnode,int level);
        void insert(Rnode row);
        void ordernodes(vector<Rnode*> &dataset,int ini,int tam,int dimens,int ordim,funcSortByDim functord);
        void bulkloadingAux(vector<Rnode*>& allPtr,int ini,int datasetsize);
        void bulkloading(vector<Rnode*> allptr,stack<int> lvl);
        bool search(Rnode row);
};

