#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <list>
#include <time.h>
#include <iomanip> 
#include <ctime>
#include <fstream>
#include <sstream>
#include <map>
#include <cstdlib>

float smoothing = 1.0;
// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory
using namespace std;

// Our graph consists of a list of nodes where each node is represented as follows:
class Graph_Node{

private:
	string Node_Name;  // Variable name 
	vector<int> Children; // Children of a particular node - these are index of nodes in graph.
	vector<string> Parents; // Parents of a particular node- note these are names of parents
	int nvalues;  // Number of categories a variable represented by this node can take
	vector<string> values; // Categories of possible values
	vector<float> CPT; // conditional probability table as a 1-d array . Look for BIF format to understand its meaning
    map<string, int> value_idx; // map value to idx

public:
	// Constructor- a node is initialised with its name and its categories
    vector<float> CPT_x;
    // string useless; //
    // float denominator;

    Graph_Node(string& name,int n,vector<string>& vals)
	{
		Node_Name=name;
		nvalues=n;
		values=vals;
	}
	string get_name()
	{
		return Node_Name;
	}
	vector<int> get_children()
	{
		return Children;
	}
	vector<string> get_Parents()
	{
		return Parents;
	}
	vector<float> get_CPT()
	{
		return CPT;
	}
    float access_CPT(int index){
        if(index > CPT.size()){
            cout<<"Index "<<index<<" out of bounds\n";
            return 0;
        }
        return CPT[index];
    }
	int get_nvalues()
	{
		return nvalues;
	}
	vector<string> get_values()
	{
		return values;
	}
	void set_CPT(vector<float>& new_CPT)
	{
		CPT.clear();
		CPT=new_CPT;
	}
    void set_Parents(vector<string>& Parent_Nodes)
    {
        Parents.clear();
        Parents=Parent_Nodes;
    }
    // add another node in a graph as a child of this node
    int add_child(int new_child_index )
    {
        for(int i=0;i<Children.size();i++)
        {
            if(Children[i]==new_child_index)
                return 0;
        }
        Children.push_back(new_child_index);
        return 1;
    }
    // compute CPT
    void compute(){
        // CPT.clear();
        int z1 = CPT_x.size();
        int z2 = nvalues;
        int z3 = z1 / z2 ;
        for(int i=0;i<z3;i++){
            float sum = 0.0;
            for(int j=0;j<z2;j++){
                sum += CPT_x[j*z3+i];
            }
            float smoothing = 0.0001;
            for(int j=0;j<z2;j++){
                CPT_x[j*z3+i]+=smoothing;
                if(CPT_x[j*z3+i] < 0.001*sum){
                    CPT_x[j*z3+i] += 0.001*sum;
                    sum += 0.001*sum;
                }
            }
            sum += smoothing*z2;
            for(int j=0;j<z2;j++){
                CPT[j*z3+i]=CPT_x[j*z3+i]/sum;
            }
        }
        // cout<<z1<<' '<<z2<<' '<<z3<<'\n';
        // for(int i=0;i<z1;i++){
        //     cout<<CPT_x[i]<<' ';
        // }
        // cout<<'\n';

    }
    void base_CPT_normalise(){
        float alpha = CPT.size() / nvalues;
        for(int i=0;i<alpha;i++){
            float sum = 0;
            int nx = 0;
            for(int j=0;j<nvalues;j++){
                if(CPT[j*alpha+i]== -1){nx++;}
                else{
                    sum += CPT[j*alpha+i];
                }
            }
            if(nx==0){continue;}
            float beta = (1-sum)/nx;
            for(int j=0;j<nvalues;j++){
                if(CPT[j*alpha+i]== -1){
                    CPT[j*alpha+i] = beta;
                }
            }
        }
    }

    void constr_map(){
        for(int i=0;i<values.size();i++){
            value_idx[values[i]] = i+1;
        }
    }

    int get_val_idx(string& val){
        if(value_idx[val]==0){
            cout<<val<<'\n';
            cout<<"value not found\n";
            return -1;
        }
        return value_idx[val]-1;
    }
};


 // The whole network represted as a list of nodes
class network{

	
private:
    map<string,int> node_index;
    map<pair<int,string>,int> converter;
public:
    vector <Graph_Node> Pres_Graph;
    string intro;
    vector <string> positions;
	int addNode(Graph_Node node)
	{
		Pres_Graph.push_back(node);
        node_index[node.get_name()] = Pres_Graph.size();
        // cout<<"Added"<<node.get_name()<<'\n';
		return 0;
	}
	int netSize()
	{
		return Pres_Graph.size();
	}
    // get the index of node with a given name
    int get_index(string& val_name)
    {
        // cout<<val_name<<'\n';
        if(node_index[val_name] == 0){
            cout<<"node not found\n";
            return -1;
        }
        return node_index[val_name]-1;
    }

    // void print_map(){
    //     for(auto it = node_index.begin();it!=node_index.end();it++){
    //         cout<<it->first<<" "<<it->second<<'\n';
    //         vector<float> cppt = Pres_Graph[it->second-1].get_CPT();
    //         for(int i=0;i<cppt.size();i++){
    //             cout<<cppt[i]<<" ";
    //         }
    //         cout<<'\n';

    //     }
    // }
    void constr_converter(){
        for(int i=0;i<Pres_Graph.size();i++){
            vector<string> values = Pres_Graph[i].get_values();
            for(int j=0;j<values.size();j++){
                converter[make_pair(i,values[j])] = j;
            }
        }
    }
    int get_index_for_CPT(int y, vector<int>& paren_nvals, vector<int>& paren_indices, vector<int>& data){
        int z = paren_nvals.size();
        int fac = 1;
        int idx = 0;
        for(int i=z-1;i>=0;i--){
            idx += fac*data[paren_indices[i]];
            fac *= paren_nvals[i];
        }
        return idx;
    }
    float get_value_from_CPT_2(int y, vector<int>& paren_nvals, vector<int>& paren_idx){
        int z = paren_nvals.size();
        int fac = 1;
        int idx = 0;
        for(int i=z-1;i>=0;i--){
            idx += fac*paren_idx[i];
            fac *= paren_nvals[i];
        }
        // vector<float> cppt = Pres_Graph[y].get_CPT();
        float z1 =  Pres_Graph[y].access_CPT(idx);
        // if(idx < 0){
        //     for(auto it = paren_nvals.begin();it!=paren_nvals.end();it++){
        //         cout<<*it<<" ";
        //     }
        //     cout<<'\n';
        //     for(auto it = paren_idx.begin();it!=paren_idx.end();it++){
        //         cout<<*it<<" ";
        //     }
        //     cout<<'\n';
        // }
        return z1;
        // return cppt[idx];
    }
    // float get_value_from_CPT(string& val_name,string& value, vector<string>& upper){
    //     vector<int> paren_nvals;
    //     vector<int> paren_idx;
    //     paren_nvals.push_back(Pres_Graph[get_index(val_name)].get_nvalues());
    //     paren_idx.push_back(Pres_Graph[get_index(val_name)].get_val_idx(value));
    //     int y = get_index(val_name);
    //     vector<string> paren = Pres_Graph[y].get_Parents();
    //     int z = paren.size();
    //     for(int i=0;i<z;i++){
    //         int alpha = get_index(paren[i]);
    //         int beta = Pres_Graph[alpha].get_nvalues();
    //         paren_nvals.push_back(beta);
    //         int gamma = Pres_Graph[alpha].get_val_idx(upper[i]);
    //         paren_idx.push_back(gamma);
    //     }
    //     int fac = 1;
    //     int idx = 0;
    //     for(int i=z;i>=0;i--){
    //         idx += fac*paren_idx[i];
    //         fac *= paren_nvals[i];
    //     }
    //     return Pres_Graph[y].access_CPT(idx);
    // }
};

network read_network(ifstream& myfile)
{
	network Alarm;
	string line;
	int find=0;
  	// ifstream myfile("alarm.bif"); 
  	string temp;
  	string name;
  	vector<string> values;
    string intro;
  	
    if (myfile.is_open())
    {
    	while (! myfile.eof() )
    	{
    		stringstream ss;
      		getline (myfile,line);
      		ss.str(line);
     		ss>>temp;
     		if(temp.compare("variable")==0)
     		{
     				ss>>name;
     				getline (myfile,line);
     				stringstream ss2;
     				ss2.str(line);
     				for(int i=0;i<4;i++)
     				{
     					ss2>>temp;
     				}
     				values.clear();
     				while(temp.compare("};")!=0)
     				{
     					values.push_back(temp);
     					ss2>>temp;
    				}
     				Graph_Node new_node(name,values.size(),values);
     				int pos=Alarm.addNode(new_node);
     		}
     		else if(temp.compare("probability")==0)
     		{
     				ss>>temp;
     				ss>>temp;
                    // list<Graph_Node>::iterator listIt;
                    // list<Graph_Node>::iterator listIt1;
     				// listIt=Alarm.search_node(temp);
                    int index=Alarm.get_index(temp);
                    ss>>temp;
                    values.clear();
     				while(temp.compare(")")!=0)
     				{
                        int par = Alarm.get_index(temp);
                        Alarm.Pres_Graph[par].add_child(index);
     					values.push_back(temp);
     					ss>>temp;
    				}
                    Alarm.Pres_Graph[index].set_Parents(values);
    				getline (myfile,line);
     				stringstream ss2;
     				ss2.str(line);
     				ss2>> temp;
     				ss2>> temp;
     				vector<float> curr_CPT;
                    string::size_type sz;
     				while(temp.compare(";")!=0)
     				{
     					curr_CPT.push_back(atof(temp.c_str()));
     					ss2>>temp;
    				}
                    Alarm.Pres_Graph[index].set_CPT(curr_CPT);
                    Alarm.Pres_Graph[index].base_CPT_normalise();
                    // cout<<"1"<<endl;
                    Alarm.Pres_Graph[index].constr_map();
                    // cout<<"2"<<endl;

     		}
            else if(temp.compare("property")==0){
                Alarm.positions.push_back(line);
            }
            else if(temp.compare("//")==0){

            }
            else
            {
                continue;
            }
     		
    		
    	}
    	Alarm.constr_converter();
    	myfile.close();
  	}
  	
  	return Alarm;
}

vector<vector<int> > read_data(network& Alarm,ifstream& myfile){
    string line;
//    ifstream myfile("reco rds.dat");
    string temp;
    string name;
    int vars = Alarm.Pres_Graph.size();
    vector<vector<int> > data;
    // cout<<"Reading data\n";
    if(myfile.is_open()){
        while(!myfile.eof()){
            // cout<<"Reading data\n";
            stringstream ss;
            int miss = -1;
            getline(myfile,line);
            ss.str(line);
            vector<int> data_pt;
            for(int i=0;i<vars;i++){
                ss>>temp;
                if(temp.compare("\"?\"")==0){
                    data_pt.push_back(-1);
                    miss = i;
                }
                else{
                    int index = Alarm.Pres_Graph[i].get_val_idx(temp);
                    data_pt.push_back(index);
                }
            }
            data_pt.push_back(miss);
            data.push_back(data_pt);
        }
        myfile.close();
    }
    return data;
}

void Expectation_Maximisation(network &alarm, vector<vector<int> >& data,time_t& start_time){
    time_t end_time;
    float time_limit = 110.0;
    bool flag = false;
    // start_time = std::time(&start_time);
    // 1. first compute a weighted data set, with weights and no ?
    vector<pair<float,vector<int> > > new_data; // weight and data set
    vector<vector<int> > parent_indices; // parent indices
    vector<vector<int> > parent_nvals; // parent indices
    vector<vector<int> > children_indices; // children indices
    int vars = alarm.Pres_Graph.size();
    for(int i=0;i<vars;i++){
        vector<int> temp;
        vector<int> temp2;
        temp.push_back(i);
        temp2.push_back(alarm.Pres_Graph[i].get_nvalues());
        vector<string> parents = alarm.Pres_Graph[i].get_Parents();
        int x = parents.size();
        for(int j=0;j<x;j++){
            int idx = alarm.get_index(parents[j]);
            temp.push_back(idx);
            temp2.push_back(alarm.Pres_Graph[idx].get_nvalues());
        }
        parent_indices.push_back(temp);
        parent_nvals.push_back(temp2);
        vector<int> temp3;
        temp3 = alarm.Pres_Graph[i].get_children();
        children_indices.push_back(temp3);
        // cout<<i<<" "<<'\n';
        // for(auto it = temp.begin();it!=temp.end();it++){
        //     cout<<*it<<" ";
        // }
        // cout<<'\n';
        // for(auto it = temp2.begin();it!=temp2.end();it++){
        //     cout<<*it<<" ";
        // }
        // cout<<'\n';
    }
    int b1 = 0;
    for(int r=0;r<500;r++){
    end_time = std::time(&end_time);
    new_data.clear();
    cout<<"TIME : "<<end_time - start_time<<'\n';
    if(end_time - start_time > time_limit){
        cout<<"Iteration number "<<r<<'\n';
        break;
    }
    int x = data.size();
    for(int i=0; i<x; i++){
        // Now we need to fill the ? with the all values and compute their weights
        // vector<int> temp = data[i];
        int index = data[i][37];
        // for(int j=0; j<vars; j++){
        //     if(data[i][j]==-1){
        //         index = j;
        //         break;
        //     }
        // }
        if(index==-1){
            new_data.push_back(make_pair(1,data[i]));
            // cout<<"Hallelujah"<<'\n';
        }
        else{
            // continue;
            // vector<string> values = alarm.Pres_Graph[index].get_values();
            int y = alarm.Pres_Graph[index].get_nvalues();
            vector<int> parents = parent_indices[index];
            vector<int> parent_vals;
            // parent_vals.push_back(temp[index]);
            for(int j=0; j<parents.size(); j++){
                parent_vals.push_back(data[i][parents[j]]);
            }
            // vector<int> children = children_indices[index];
            vector<float> weights;
            float wt1 = 1.0;
            for(int j=0; j<y; j++){
                // We need to find the weight using the Markov Blanket
                // 1. Find the parents
                data[i][index] = j;
                parent_vals[0] = j;
                float wt = alarm.get_value_from_CPT_2(index,parent_nvals[index],parent_vals); // this is the weight due to my parents
                // cout<<wt<<endl;
                // 2. Now I need to multiply this with the weights from my childrens parents
                for(int k=0;k<children_indices[index].size();k++){
                    // vector<int> child_parents = parent_indices[children_indices[index][k]];
                    vector<int> child_parent_vals;
                    // child_parent_vals.push_back(temp[children[k]]);
                    for(int l=0;l<parent_indices[children_indices[index][k]].size();l++){
                        if(parent_indices[children_indices[index][k]][l]==index){
                            child_parent_vals.push_back(j);
                        }
                        else{
                            child_parent_vals.push_back(data[i][parent_indices[children_indices[index][k]][l]]);
                        }
                    }
                    // cout<<children[k]<<'\n';
                    // for(auto z : parent_nvals[children[k]]){
                    //     cout<<z<<" ";
                    // }
                    // cout<<'\n';
                    // for(auto z : child_parent_vals){
                    //     cout<<z<<" ";
                    // }
                    // cout<<'\n';
                    wt1 = alarm.get_value_from_CPT_2(children_indices[index][k],parent_nvals[children_indices[index][k]],child_parent_vals);
                    // cout<<wt1<<'\n';
                    wt *= wt1;
                }
                // cout<<wt<<'\n';
                weights.push_back(wt);
            }
            data[i][index] = -1;
            // now normalise the weights;
            float sum = 0;
            for(int j=0;j<y;j++){
                sum += weights[j];
            }
            for(int j=0;j<y;j++){
                weights[j] /= sum;
            }
            // now add this to the new data array
            for(int j=0;j<y;j++){
                data[i][index] = j;
                new_data.push_back(make_pair(weights[j],data[i]));
                // cout<<weights[j]<<" ";
                // cout<<alarm.get_index_for_CPT(index,parent_nvals[index],parent_indices[index],data[i])<<' ';
                // for(auto it = data[i].begin();it!=data[i].end();it++){
                    // cout<<*it<<" ";
                // }
                // cout<<'\n';
            }
            data[i][index] = -1;
        }
    }
    // 2. Use that dataset to compute a new CPT_x
    int b2 = new_data.size();
    if(!flag){
    for(int i=0;i<vars;i++){
        vector<float> alpha = alarm.Pres_Graph[i].get_CPT();
        int beta = alpha.size();
        int gamma = alarm.Pres_Graph[i].get_nvalues();
        for(int j=0;j<beta;j++){
            alarm.Pres_Graph[i].CPT_x.push_back(0.05/gamma);
        }
        // alarm.Pres_Graph[i].denominator = 0.1*beta;
    }
    flag = true;
    }
    for(int i=0;i<b2;i++){
        for(int j=0;j<vars;j++){
            int idx = alarm.get_index_for_CPT(j,parent_nvals[j],parent_indices[j],new_data[i].second);
            alarm.Pres_Graph[j].CPT_x[idx] += new_data[i].first;
            // cout<<new_data[i].first<<" "<<'\n';
            // alarm.Pres_Graph[j].denominator += new_data[i].first;
        }
    }
    // b1 = b2;
    // 3. using the CPT_x compute and update the CPT and delete the new dataset
    for(int i=0;i<vars;i++){
        alarm.Pres_Graph[i].compute();
        // alarm.Pres_Graph[i].CPT_x.clear();
    }

    // 4. Repeat until convergence
    }
    //output the new CPT's
    // for(int i=0;i<37;i++){
    //     vector<float> cppt = alarm.Pres_Graph[i].get_CPT();
    //     for(int j=0;j<cppt.size();j++){
    //         cout<<cppt[j]<<" ";
    //     }
    //     cout<<'\n';
    // }
}
// void smooth_cpt(vector<float>& CPT, int nvals){
//     int z1 = CPT.size();
//     int z2 = nvals;
//     int z3 = z1/z2;
//     for(int i=0;i<z3;i++){
//         // float sum = 0.0;
//         // for(int j=0;j<z2;j++){
//         //     sum += CPT[j*z3+i];
//         // }
//         // float smoother = 0.001*sum;
//         // for(int j=0;j<z2;j++){
//         //     CPT[j*z3+i]+=smoother;
//         // }
//         // sum += smoother*z2;
//         for(int j=0;j<z2;j++){
//             // CPT[j*z3+i]=CPT[j*z3+i]/sum;
//             if()
//         }
//     }
// }


void write_data(network& alarm){
    ofstream outfile("solved_alarm.bif");
    string output;
    outfile<<"// Bayesian Network in the Interchange Format\n";
    outfile<<"// Produced by BayesianNetworks package in JavaBayes\n";
    outfile<<"// Output created Sun Nov 02 17:58:15 GMT+00:00 1997\n";
    outfile<<"// Bayesian network \n";
    outfile<<"network \"Alarm\" { //37 variables and 37 probability distributions\n}\n";
    // write the variables
    for(int i=0; i<37; i++){
        outfile<<"variable  ";
        outfile<<alarm.Pres_Graph[i].get_name();
        outfile<<" ";
        outfile<<"{ //";
        outfile<<to_string(alarm.Pres_Graph[i].get_nvalues());
        outfile<<" values\n";
        outfile<<"\ttype discrete[";
        outfile<<to_string(alarm.Pres_Graph[i].get_nvalues());
        outfile<<"] { ";
        vector<string> values = alarm.Pres_Graph[i].get_values();
        for(int j=0; j<values.size(); j++){
            outfile<<" ";
            outfile<<values[j];
            outfile<<" ";
        }
        outfile<<"};\n";
        outfile<<alarm.positions[i];
        outfile<<"\n";
        outfile<<"}\n";

    }
    // write the probability
    for(int i=0;i<37;i++){
        outfile<<"probability (  ";
        outfile<<alarm.Pres_Graph[i].get_name();
        vector<string> Parents = alarm.Pres_Graph[i].get_Parents();
        for(int j=0;j<Parents.size();j++){
            outfile<<"  ";
            outfile<<Parents[j];
            // output+=" ";
        }
        outfile<<" ) { //";
        vector<string> values = alarm.Pres_Graph[i].get_values();
        int y1 = Parents.size();
        outfile<<to_string(y1+1);
        vector<float> cppt = alarm.Pres_Graph[i].get_CPT();
        // smooth_cpt(cppt,alarm.Pres_Graph[i].get_nvalues());
        // smooth_cpt(cppt,alarm.Pres_Graph[i].get_nvalues());
        // smooth_cpt(cppt,alarm.Pres_Graph[i].get_nvalues());
        outfile<<" variable(s) and ";
        int x1 = cppt.size();
        outfile<<to_string(x1);
        outfile<<" values";
        outfile<<"\n\ttable ";
        for(int j=0;j<cppt.size();j++){
            outfile<<" ";
            float y = cppt[j];
            if(y < 0.01){
                y = 0.01;
            }
            if(y > 0.99){
                y = 0.99;
            }
            outfile<<fixed<<setprecision(4)<<y;
            outfile<<" ";
        }
        outfile<<";\n}\n";
    }

}

int main(int argc, char** argv)
{
    time_t start_time;
    start_time = std::time(&start_time);
	network Alarm;
    std::ifstream infile(std::string(argv[1]) + "");
    std::ifstream infile2(std::string(argv[2]) + "");
	Alarm=read_network(infile);
	cout<<"Perfect! Hurrah! \n";
    vector<vector<int> > data;
    data = read_data(Alarm,infile2);
    Expectation_Maximisation(Alarm,data,start_time);
    time_t end_time = std::time(&end_time);
    cout<<"TIME : "<<end_time - start_time<<'\n';
    write_data(Alarm);
    // for(int i=0;i<37;i++){
    //     cout<<data[1][i]<<" ";
    // }
    // cout<<'\n';
    // cout<<data.size()<<endl;
    // // Alarm.print_map();
    // string x = "\"HRBP\"";
    // // string t = "\"High\"";
    // vector<int> upper;
    // upper.push_back(1);
    // upper.push_back(2);
    // // cout<<Alarm.get_value_from_CPT(x,t,upper)<<'\n';
    // cout<<Alarm.get_value_from_CPT_2(x,2,upper)<<'\n';
}