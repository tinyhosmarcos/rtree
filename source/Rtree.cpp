#include "Rtree.h"
using namespace std;
bool compare(vector<double> a,vector<double> b){
	if(a.size() != b.size()){return false;}
	for(int i=0;i<a.size();i++){
		if(a[i]!=b[i]){return false;}
	}
	return true;
}
Rnode::Rnode(int dimensions_, bool isleaf_ ){
		position.resize(dimensions_,0);
		size.resize(dimensions_,0);
		dimensions=dimensions_;
		isleaf=isleaf_;
	}
void Rnode::insertNode(Rnode row){
                Rnode* ptr;
                ptr=new Rnode(dimensions,false);
                ptr->position=row.position;
                ptr->size=row.size;
                childs.push_back(ptr);
	}
void Rnode::insertNodePtr(Rnode* row){
            childs.push_back(row);
        }
bool Rnode::intersects(Rnode Nrow){
                vector<double> row=Nrow.position;
		for(int i=0;i<dimensions;i++){
			if(row[i]<position[i] or row[i]>position[i]+size[i]){
				return false;
			}
		}
		return true;
	}
double Rnode::interseccion(Rnode Nrow){
                vector<double> pos=Nrow.position;
                vector<double> sz=Nrow.size;
		vector<double> postmp(dimensions);
		vector<double> sztmp(dimensions);
		for(int i=0;i<dimensions;i++){
			if(position[i]<pos[i]){
				if(pos[i]<position[i]+size[i]){
					postmp[i]=pos[i];
					if(position[i]+size[i]<pos[i]+sz[i]){
						sztmp[i]=position[i]+size[i]-pos[i];
					}
					else{
						sztmp[i]=sz[i];
					}
				}
				else{return 0;}
			}
			else{
				if(position[i]<pos[i]+sz[i]){
					postmp[i]=position[i];
					if(pos[i]+sz[i]<position[i]+size[i]){
						sztmp[i]=pos[i]+sz[i]-position[i];	
					}
					else{
						sztmp[i]=size[i];
					}
					
				}
				else{return 0;}
			}
		}
double rpta=1;
		for(int i=0;i<dimensions;i++){
			rpta= rpta*sztmp[i];
		}
		return rpta;

	}
double Rnode::crecimiento(Rnode Nrow){
                vector<double> row=Nrow.position;
		double rpta =1;
		for(int i=0;i<dimensions;i++){
			if(row[i]< position[i]){rpta=rpta * (size[i]+abs(position[i]-row[i]));}
			else{
				if(row[i]> position[i]+size[i]){rpta= rpta* abs(position[i]-row[i]);}
				else{rpta = rpta*size[i];}
			}
		}
		return rpta-area();
	}
double Rnode::overlapcrecimiento(int index,Rnode Nrow){
                vector<double> row=Nrow.position;
		vector<double> postmp=childs[index]->position;
		vector<double> sztmp=childs[index]->size;
		for(int i=0;i<dimensions;i++){//crecimiento
			if(row[i]< childs[index]->position[i]){postmp[i]=row[i];}
			else{
				if(row[i]> childs[index]->position[i]+childs[index]->size[i]){
					sztmp[i]=row[i]-childs[index]->position[i];}
			}
		}
		double rpta=0;
		for(int i=0;i<childs.size() and i!=index;i++){//overlap de crecido - overlap de notmal
                        Rnode tmp(dimensions,false);
                        tmp.position=postmp;
                        tmp.size=sztmp;
                        double a=childs[i]->interseccion(tmp);
                        double b=childs[i]->interseccion(*(childs[index]));
			rpta=rpta+a-b;
		}
		return rpta;
	}
void Rnode::intentacrecer(Rnode Nrow){
                vector<double> row=Nrow.position;
		for(int i=0;i<dimensions;i++){
			if(row[i]< position[i]){position[i]=row[i];}
			if(row[i]> position[i]+size[i]){size[i]=abs(position[i]-row[i]);}				
		}
	}
double Rnode::area(){
		double rpta=1;
		for(int i=0;i<dimensions;i++){
			rpta= rpta*size[i];
		}
		return rpta;
	}
        //not tested
void Rnode::actMBR(){
            double min=DBL_MIN,max=DBL_MAX;
            for(int i=0;i<dimensions;i++){
                min=DBL_MAX;max=DBL_MIN;
                for(int j=0;j<childs.size();j++){
                    if(childs[j]->position[i]<min){min=childs[j]->position[i];}
                    if(childs[j]->position[i]+childs[j]->size[i] > max){max=childs[j]->position[i]+childs[j]->size[i];}
                }
                position[i]=min;
                size[i]=max-min;
            }
        }
        


funcSortByDim::funcSortByDim(){
        dim=0;maxdim=0;base=0;
    }
    
bool funcSortByDim::operator()(Rnode* vec1,Rnode* vec2){
        return vec1->position[dim]+(vec1->size[dim]/2) < vec2->position[dim]+(vec2->size[dim]/2);
    }
void funcSortByDim::up(){
        dim=((dim+1)%maxdim)+base;
    }



Rtree::Rtree(int M_,int m_,int dimensions_){
                root = new Rnode(dimensions_,true);
                //root=0;
		M = M_;
		m = m_;
		dimensions=dimensions_;
                funct.maxdim=dimensions_;
                OTlastcallLvl=-1;//Overflow treatment last call level
	}

void Rtree::tryprint(){
            queue<Rnode*> lista;
            lista.push(root);
            int c=lista.size();
            int nds=0;
            while(!lista.empty()){

                for(int i=0;i<lista.front()->childs.size();i++){
                    lista.push(lista.front()->childs[i]);
                    nds++;
                }
                lista.pop();
                c--;
                if(c==0){
                    cout<<nds<<endl;
                    nds=0;
                    c=lista.size();
                }
            }
        }

Rnode* Rtree::chooseSubTree(Rnode* actualnode, Rnode Nrow, int &level){
                vector<double> row=Nrow.position;
		if(actualnode->isleaf){return actualnode;}
		else{
			if(actualnode->childs[0]->isleaf){
				//determina minimun overlap cost
				int index=0;
				double value=DBL_MAX;
				bool notoverlapzero=true;
				bool notcrecimientozero=true;
				double tmp;
				for(int i=0;i<actualnode->childs.size();i++){
                                        if(notoverlapzero){tmp=actualnode->overlapcrecimiento(i,Nrow);}
					if(!notoverlapzero or tmp==0){//resolver ataduras overlap=0 por menor crecimiento
						if(notoverlapzero){notoverlapzero=false;value=DBL_MAX;}
                                                if(notcrecimientozero){tmp=actualnode->childs[i]->crecimiento(Nrow);}
						if(!notcrecimientozero or tmp==0){//resolver ataduras crecimiento=0 por menor area
							if(notcrecimientozero){notcrecimientozero=false;value=DBL_MAX;}
							tmp=actualnode->childs[i]->area();
							if(tmp<value){
								index = i;
								value = tmp;
							}
						}
						if(notcrecimientozero and tmp<value){
							index=i;	
							value=tmp;
						}
					}
					if(notoverlapzero and tmp<value){
						index = i ;
						value = tmp;
					}
				}
                                if(notoverlapzero or notcrecimientozero){actualnode->intentacrecer(Nrow);}
				actualnode=actualnode->childs[index];
			}
			else{
				//determina minimun area cost
				int index=0;
				double value=DBL_MAX;
				bool notcrecimientozero=true;
				double tmp;
				for(int i=0;i<actualnode->childs.size();i++){						
                                        if(notcrecimientozero){tmp=actualnode->childs[i]->crecimiento(Nrow);}
					if(!notcrecimientozero or tmp==0){//resolver ataduras crecimiento=0 por menor area
						if(notcrecimientozero){notcrecimientozero=false;value=DBL_MAX;}
						tmp=actualnode->childs[i]->area();
						if(tmp<value){
							index = i;
							value = tmp;
						}
					}
					if(notcrecimientozero and tmp<value){
						index=i;	
						value=tmp;
					}
					
				}
                                if(notcrecimientozero){actualnode->intentacrecer(Nrow);}//si el MBR hijo tuvo que crecer talvez el tngaque crecer yo pambien
				actualnode=actualnode->childs[index];
			}
		}
		level++;
                return chooseSubTree(actualnode,Nrow,level);
	}
	
void Rtree::reinsert(Rnode* N){
//            RI1 For all M+1 entries of a node N, compute the distance
//            between the centers of their rectangles and the center of the
//            bounding rectangle of N
            vector<double> distances(dimensions);
            for(int i=0;i<N->childs.size();i++){

            }
//            RI2 Sort the entries in decreasing order of their distances
//            computed in RI1

//            RI3 Remove the first p entries from N and adjust the bounding
//            rectangle of N

//            RI4 In the sort, defined in RI2, starting with the maximum
//            distance (= far reinsert) or minimum distance (= close
//            reinsert), invoke Insert to reinsert the entries
        }

void Rtree::overflowTreatment(Rnode* actualnode,int level){

//		If the level is not the root level and this is the first call of
//		OverflowTreatment in the given level during the insertion of
//		one data rectangle, then
		if(level!=0 and level!=OTlastcallLvl){
			OTlastcallLvl= level;
			//reinsert
		}
		else{
			//split
		}

	}

void Rtree::insert(Rnode row){
		int level=0;
		Rnode* N = chooseSubTree(root,row,level);
                N->insertNode(row);
                if(N->childs.size()==M+1){
			overflowTreatment(N,level);
		}
	}
	
void Rtree::ordernodes(vector<Rnode*> &dataset,int ini,int tam,int dimens,int ordim,funcSortByDim functord){

            int P=ceil(double(tam)/double(M));

            int NpS=M*ceil(pow(double(P),double(double(dimens-1)/double(dimens))));//Nodes per Stab
            sort(dataset.begin()+ini,dataset.begin()+ini+tam,functord);
            if(functord.dim==dimensions-1){return;}//si ya ordenaste por todas las dimensiones
            functord.up();
            for(int i=0;i<tam;i+=NpS){//por cada Stab llamar recrusivamente
                if(i+NpS>tam){
                    ordernodes(dataset,i,tam%NpS,dimens-1,ordim+1,functord);
                }
                else{
                    ordernodes(dataset,i,NpS,dimens-1,ordim+1,functord);
                }
            }

        }

void Rtree::bulkloadingAux(vector<Rnode*>& allPtr,int ini,int datasetsize){

            ordernodes(allPtr,ini,datasetsize,dimensions,0,funct);//ordenar por dimension x_i
            int index=ini+datasetsize-1;
            int i;
            for(i=0;i<datasetsize;i++){//todos los elementos del dataset

                if(i%M==0){
                    if(i!=0){
                        allPtr[index]->actMBR();
                    }
                    index++;
                }

                allPtr[index]->insertNodePtr(allPtr[ini+i]);//se guardan en otros nodos
            }
            allPtr[index]->actMBR();
            if(index-ini-datasetsize<M){return;}
            bulkloadingAux(allPtr,ini+datasetsize,index-(ini+datasetsize-1));
            //llamar de nuevo con nuevos indices para el dataset
        }

void Rtree::bulkloading(vector<Rnode*> allptr,stack<int> lvl){
            int index=0;
            int tmp=0;
            queue<Rnode*> currins;
            currins.push(root);
            while(!lvl.empty()){
                for(int i=0;i<lvl.top();i+=M){
                    for(int j=0;j<M and i+j<lvl.top();i++){
                        currins.front()->insertNodePtr(allptr[index]);
                        currins.push(allptr[index]);
                        index++;
                    }
                    currins.pop();
                }
                cout<<index-tmp<<endl;
                tmp=index;
                lvl.pop();
            }

            root->actMBR();
            cout<<"total indeado: "<<index<<endl;
        }

bool Rtree::search(Rnode row){
            queue<Rnode*> lista;
            lista.push(root);
            while(!lista.empty()){
                if(lista.front()->childs.size()==0){return true;}
                if(compare(lista.front()->position,row.position)){return true;}
                for(int i=0;i<lista.front()->childs.size();i++){
                    if(lista.front()->childs[i]->intersects(row)){
                        lista.push(lista.front()->childs[i]);
                    }
                }
                lista.pop();
            }
            return false;
        }


int getData2(string filename,vector<vector<string > > &data){
    data.clear();
    ifstream file(filename);
    string info;
    int c=0;
    int dim=0;

    getline(file,info,'\n');
    stringstream sst(info);
    while(getline(sst, info, ',')){
        dim++;
    }
    vector<string> tmp(dim);
    while(file.good()){
        getline(file,info,'\n');
        stringstream ss(info);
        while(getline(ss, info, ',')){
            c=c%dim;
            int i,f;
            for(i=0;info[i]=='"' or info[i]=='\n';i++){}
            for(f=info.size();info[f-1]=='"' or info[f-1]=='\n';f--){}
            tmp[c]=info.substr(i,f-i);
            if(c==dim-1){data.push_back(tmp);}
            c++;
        }
    }
    return dim;
}

void string2Rnode(vector<vector<string > >& input,vector<Rnode*>& output){

    int dim=input[0].size();
    for(int i=0;i<input.size();i++){
        output[i]=new Rnode(dim,false);
        for(int j=0;j<dim;j++){
            //cout<<"trans:"<<stod(input[i][j])<<endl;
            output[i]->position[j]=stod(input[i][j]);
            //cout<<output[i]->position[j]<<endl;
        }
    }
}

