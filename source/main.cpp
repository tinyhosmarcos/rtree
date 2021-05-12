#include "Rtree.cpp"
using namespace std;
int main(){
    int M=10,n=5;
    cout<<"lectura..."<<endl;
    vector<vector<string > > data_file;
    int dimensions=getData2("data_points.txt",data_file);
    cout<<"lectura terminada"<<endl;
    cout<<"transformando..."<<endl;
    vector<Rnode*> dataset(data_file.size());
    string2Rnode(data_file,dataset);
    cout<<"transformacion terminada..."<<endl;
    cout<<"SIZE: "<<dataset.size()<<endl;
    cout<<"pre bulk..."<<endl;
    stack<int> lvl;
    int allptrsize=dataset.size();
    lvl.push(dataset.size());
    while(lvl.top()>M){
        int tmp=ceil(double(lvl.top())/double(M));
        lvl.push(tmp);
        allptrsize+=tmp;
    }
    cout<<"pre bulk finished"<<endl;
    cout<<"inicializacion..."<<endl;
    vector<Rnode*> allptr(allptrsize);
    for(int i=0;i<allptrsize;i++){
        allptr[i]=new Rnode(dimensions,false);
    }
    for(int i=0;i<dataset.size();i++){
        allptr[i]->position=dataset[i]->position;
    }
    cout<<"inicializacion terminada"<<endl;
    Rtree a(M,n,dimensions);
    cout<<"bulk aux"<<endl;
    a.bulkloadingAux(allptr,0,dataset.size());
    cout<<"bulk aux end"<<endl;

    cout<<"bulk load..."<<endl;
    a.bulkloading(allptr,lvl);
    cout<<"indexacion completa"<<endl;




    return 0;
}
