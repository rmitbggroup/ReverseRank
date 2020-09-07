#include <iostream>
#include <limits.h>
#include <string.h>
#include <stdlib.h>           /* C standard library                   */
#include <iostream>
#include <fstream>
#include <cstdio>
#include <ctime>
#include <queue>
#include <algorithm>
#include <map>
#include <cmath>
using namespace std;

// Number of vertices in given graph

#define N 500
const int M = 14000;
#define V N+2
const int K=100;

double T[N][N];
double e[N][N];
double G[N][N];
double standG[N][N];
double  graph[V][V];
double x[N],y[N], w[N], nw[N];
using namespace std;
double total_time=0.0;
double ct_lp =0.0;

string path ="C:/Users/liang/PycharmProjects/whynotquery/Wrap/Netflix/";
string xp =path+"x"+to_string(N)+".txt" ;
string yp =path+"y"+to_string(N)+".txt" ;
string wp =path+"w"+to_string(N)+".txt" ;
string tp = path+to_string(N)+"T_3miu0.5.txt";



/* Returns true if there is a path from source 's' to sink 't' in
  residual graph. Also fills parent[] to store the path */
bool bfs(int s, int t, int parent[])
{
    // Create a visited array and mark all vertices as not visited
    bool visited[V];
    memset(visited, 0, sizeof(visited));

    // Create a queue, enqueue source vertex and mark source vertex
    // as visited
    queue <double> q;
    q.push(s);
    visited[s] = true;
    parent[s] = -1;

    // Standard BFS Loop
    while (!q.empty())
    {
        int u = q.front();
        q.pop();

        for (int v=0; v<V; v++)
        {
            if (visited[v]==false && graph[u][v] > 0)
            {
                q.push(v);
                parent[v] = u;
                visited[v] = true;
            }
        }
    }

    // If we reached sink in BFS starting from source, then return
    // true, else false
    return (visited[t] == true);
}

// Returns the maximum flow from s to t in the given graph
double fordFulkerson( int s, int t)
{
    int u, v;

    // Create a residual graph and fill the residual graph with
    // given capacities in the original graph as residual capacities
    // in residual graph
//    double rGraph[V][V]; // Residual graph where rGraph[i][j] indicates
    // residual capacity of edge from i to j (if there
    // is an edge. If rGraph[i][j] is 0, then there is not)
//    for (u = 0; u < V; u++)
//        for (v = 0; v < V; v++)
//            rGraph[u][v] = graph[u][v];

    int parent[V];  // This array is filled by BFS and to store path

    double max_flow = 0;  // There is no flow initially

    // Augment the flow while tere is path from source to sink
    while (bfs(s, t, parent))
    {
//        cout<<s<<" "<<t<<endl;
        // Find minimum residual capacity of the edges along the
        // path filled by BFS. Or we can say find the maximum flow
        // through the path found.
        double path_flow = 1e10;
        for (v=t; v!=s; v=parent[v])
        {
            u = parent[v];
            path_flow = min(path_flow, graph[u][v]);
        }

        // update residual capacities of the edges and reverse edges
        // along the path
        for (v=t; v != s; v=parent[v])
        {
            u = parent[v];
            e[u][v]+=path_flow;
            graph[u][v] -= path_flow;
            graph[v][u] += path_flow;
        }

        // Add path flow to overall flow
        max_flow += path_flow;
    }

    // Return the overall flow
    return max_flow;
}



double gettarget()
{
    double sum=0.0;
    for(int i=0;i<N;i++)
    {
        sum+=1.0*w[i]*(x[i]-y[i]);
    }
    return sum;
}


void readT()
{
    int n;
    ifstream infile;

    infile.open (tp);
    infile >> n ;
//    cout<<n<<" ";
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            T[i][j]=0.5;
//            T[i][j] = 1.0;
//            cerr<<T[i][j];
        }
    }
}

void readxyw()
{
    ifstream infilex,infiley,infilew;
//    cerr<<wp<<endl;
    infilex.open (xp);
    infiley.open (yp);
    infilew.open (wp);

    for(int i=0;i<N;i++)
    {
        infiley>>y[i];
        infilex>>x[i];
        infilew>>w[i];
        nw[i]=w[i];
//        cerr<<w[i]<<endl;
    }

}
//C:\Users\liang\CLionProjects\mwa_greedy\mid_x\x_2
double item[14000][N];
double wg[N], aw[N];
double getorank()
{
    double sum=0.0;
    for(int i=0;i<N;i++)
    {
        sum+=1.0*w[i]*(x[i]-y[i]);
    }
    return sum;
}

double getnrank()
{
    double sum=0.0;
    for(int i=0;i<N;i++)
    {
        sum+=1.0*wg[i]*(x[i]-y[i]);
    }
    return sum;
}

void updatenw()
{
    for(int i=0;i<N;i++)
    {
        nw[i]=w[i];
    }
}

void readembedding()
{
    readT();

}

struct dpair{
    int x,y;
    double diff;
};

bool compare(dpair a, dpair b)
{
    return a.diff<b.diff;
}

double cost=0.0;
dpair P[N*N];
int ct=0;
dpair ans[N*N];
void ConstructFlow()
{
    clock_t start;
    double duration;

    start =  clock();
//    readxyw();
//    readT();
    updatenw();

    if(gettarget()<=0)
    {
        for(int i=0;i<N;i++)
        {
            swap(x[i], y[i]);
        }
    }

    if(gettarget()<=0)
    {
//        cerr<<gettarget()<<endl;
        return;
    }
    ct_lp+=1.0;
//    cerr<<"flow"<<endl;
    memset(graph, 0, sizeof(graph));
    memset(e, 0, sizeof(e));
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            graph[i][j] = T[i][j];
        }
    }



    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            dpair t;
            t.diff = -x[i] + y[i] + x[j]  - y[j];
            t.x = i;
            t.y = j;
            P[ct] = t;
            ct++;
        }
    }
    sort(P, P+ct, compare);
    double target = gettarget();


    for(int i=0;i<N*N;i++)
    {
//        cout<<P[i].diff<<endl;

        int source = P[i].x;
        int sink = P[i].y;
        double eff= P[i].diff;
        if(eff>=0)
            break;
//        cout<<source<<" "<<sink<<" "<<eff<<endl;
        graph[N][source] = 1+nw[source];
        graph[sink][N+1] = 1- nw[sink];
        double flow = fordFulkerson( N, N+1);
        graph[N][source] = 0;
        graph[sink][N+1]=0;

//        cout << "The maximum possible flow is " << flow <<endl;
        if(flow<1e-10)
            continue;
        if(abs(flow*eff) < target)
        {
//            cout<<source<<" "<<sink<<" "<<flow<<endl;
            cost += flow;
            target -= abs(flow*eff);
            nw[source] -= flow;
            nw[sink] += flow;
            G[source][sink]+=flow;
            standG[source][sink]+=flow*flow;
        }
        else
        {

            flow = abs(target/eff);
//            cout<<source<<" "<<sink<<" "<<flow<<endl;
            cost += flow;
            target -= abs(flow*eff);
            nw[source] -= flow;
            nw[sink] += flow;
            G[source][sink]+=flow;
            standG[source][sink]+=flow*flow;
            break;
        }
//        if(flow>1e-10)
//        cout<<source<<" "<<sink<<" cost:"<<cost<<"eff"<<eff<<endl;
    }

//    cout<<cost<<" "<<getntarget()<<endl;
//    cout<<CostFunct()<<endl;
    duration = (  clock() - start ) / (double) CLOCKS_PER_SEC;
//    cout<<"Time: "<< duration <<endl;
    total_time+=duration;

}


struct something{
    int id;
    double score;
}rg[M],ri[M];


bool com(something t1, something t2)
{
    return t1.score>t2.score;
}
map<int,int> topset;
double Measure(int k)
{
    double ans=0;
    for(int i=0;i<M;i++)
    {
        rg[i].id=i;
        rg[i].score=0;
        ri[i].id=i;
        ri[i].score=0;
        for(int j=0;j<N;j++)
        {
            rg[i].score+=item[i][j]*wg[j];
            ri[i].score+=item[i][j]*w[j];
        }
    }
    sort(rg,rg+M, com);
    sort(ri,ri+M, com);


    topset.clear();
    for(int i=0;i<k;i++)
    {
        topset[rg[i].id]=1;
    }

    for(int i=0;i<k;i++)
    {
        if(topset[ri[i].id]==1)
        {
            ans+=1.0;
        }
    }
    return ans/k;

}

void updatew()
{
    for(int i=0;i<N;i++)
    {
        w[i]=nw[i];
    }
}


int a[101],b[101];


void readembed(int u_id)
{
//    cerr<<u_id<<endl;
//    string respath = path+"/res/r"+to_string(u_id);
//    string pairspath=path+"pairs";
    ofstream myfile, hitsp, hamfile;

    string respath=path+"/mf/mfwre"+to_string(u_id);
    myfile.open(respath);
    ifstream infileaw,infileitem,infilew;
    tp=path+"G";
    string awp=path+"gl";
//    string itemp=path+"/top1000/i"+to_string(u_id);
    string itemp= path+"item";
//    string wgp = path+"/top1000/u"+to_string(u_id);
    string wgp=path+"100user";
//    tp = path+"g0.5";
    infileaw.open (awp);
    infileitem.open (itemp);
    infilew.open (wgp);
    hamfile.open(path+"/mf/mfham"+to_string(u_id));
//    inpair.open(pairspath);
    readT();
    hitsp.open(path+"/mf/mfhits"+to_string(u_id));

//    for(int i=0;i<100;i++)
//    {
//        inpair>>a[i];
//        inpair>>b[i];
////        cerr<<a[i]<<" "<<b[i]<<endl;
//    }
    for(int i=0;i<N;i++)
    {
        infileaw>>aw[i];
            w[i]=aw[i];
//        w[i]=1.0/N;
    }

    for(int j=0;j<=u_id;j++)
    {
        for(int i = 0; i<N;i++)
        {
            infilew>>wg[i];
        }
    }



    for(int i=0;i<M;i++)
    {
        for(int j=0;j<N;j++)
        {
            infileitem >> item[i][j];
        }
    }
    srand(time(0));
    int cx[N], cy[N], candidate=0;


    for(int step=1;step<=2000;step++) {
        cout<<"Round:"<<step<<endl;
        int ct = 0;
        candidate = 0;
        for (int loop = 1; loop <= 1; loop++) {


            ct++;
            while (1) {
                int ta = rand() % M;
                if(topset[ta]==1)
                    continue;
                int tb = rand() % K;
                for (int i = 0; i < N; i++) {
                    x[i] = item[ta][i];
                    y[i] = item[rg[tb].id][i];
                }

                if (getorank() * getnrank() < 0) {
//                    cout<<a<<" "<<b<<endl;
                    break;
                }
            }
//            cerr<<"find a,b"<<endl;
            if (getorank() < 0) {
                for (int i = 0; i < N; i++) {
                    swap(x[i], y[i]);
                }
            }
//            cerr<<getorank()<<endl;
//            cerr<<getnrank()<<endl;
            ConstructFlow();
//            printw();
//            cout << "cost:" << CostFunct() << endl;
            int ham=0;
            for(int i=0;i<N;i++)
            {
                if(abs(w[i] - nw[i])>1e-6)
                    ham++;
            }
            updatew();
            hamfile<<ham<<endl;
        }
//        cerr<<step<<endl;
//        cout<<"count: "<<ct<<endl<<endl;
//        cerr<<getorank()<<endl;
        hitsp<<Measure(K)<<endl;
        for(int i=0;i<N;i++)
        {
            myfile<<w[i]<<" ";
        }
        myfile<<endl;
    }
    hitsp.close();
    hamfile.close();
    myfile.close();
}

double total_cost=0;
void save_e(int ob, int t_size, string xy, string cap)
{
    string path = "C:/Users/liang/CLionProjects/mwa_greedy/e/";
    path+=to_string(ob)+"_T"+to_string(t_size)+xy+cap;
    ofstream my;
    my.open(path);
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            my<<e[i][j]<<" ";
        }
        my<<endl;
    }
}
double MDACost()
{
    double s=0;
    for(int i=0;i<N;i++)
    {
        double t= abs(w[i]-nw[i]);
        if(t>1e-9)
            s+=1.0;
    }
    return s;
}
//void Scale()
//{
//    string path = "C:/Users/liang/CLionProjects/mwa_greedy/data/";
//    xp = path+"x"+to_string(N)+".txt";
//    yp = path+"y"+to_string(N)+".txt";
//    wp = path+"w"+to_string(N)+".txt";
//    tp=path+to_string(N)+"T_2miu0.1.txt";
//    cerr<<wp<<endl;
//    readT();
//    readxyw();
//    ConstructFlow();
////    save_e(ob, t_size, xy, cap);
//    total_cost+=CostFunct();
//    cerr << total_time / ct_lp << endl;
//    cerr<<total_cost/ct_lp<<endl;
//}
void Tests(int ob, int t_size, string xy, string cap)
{
    string path = "C:/Users/liang/CLionProjects/mwa_greedy/high_y";
    xp ="C:/Users/liang/CLionProjects/mwa_greedy/"+xy+"_x/x_"+to_string(ob);
    yp = "C:/Users/liang/CLionProjects/mwa_greedy/"+xy+"_y/y_"+to_string(ob);
    wp = "C:/Users/liang/CLionProjects/mwa_greedy/w/w";
    tp="C:/Users/liang/CLionProjects/mwa_greedy/T/"+to_string(N)+"T_"+to_string(t_size)+"miu"+cap+".txt";
    readT();
    readxyw();
    ConstructFlow();
    save_e(ob, t_size, xy, cap);
    total_cost+=MDACost();

}

void writedata(int ts, string rank, string cap) {
    ofstream myfile, yourfile;
    string p = "C:/Users/liang/PycharmProjects/whynotquery/Experiment Results/TS=" + to_string(ts) + "/cap=" + cap + "/" +
            rank + "/mda_flow_time";
    string costp = "C:/Users/liang/PycharmProjects/whynotquery/Experiment Results/TS=" + to_string(ts) + "/cap=" + cap + "/" +
               rank + "/mda_flow_cost";
    myfile.open(p);
    yourfile.open(costp);
//    cerr << p << endl;
    ct_lp = 0.0;
    total_time = 0.0;
    total_cost=0;
    for (int i = 0; i < 50; i++) {
        if (i == 43)
            continue;
        Tests(i, ts, rank, cap);
    }
    cerr << total_time / ct_lp << endl;
    myfile << total_time / ct_lp;
    myfile.close();
    cerr << total_cost / ct_lp << endl;
    yourfile<<total_cost/ct_lp<<endl;
    yourfile.close();

}



void GenerateG(string dataset)
{
    memset(G,0, sizeof(G));
    memset(standG, 0, sizeof(standG));
    ofstream   outG, outS;

    path="C:/Users/liang/PycharmProjects/whynotquery/Wrap/500/";
    path+=dataset+"/";
    outG.open(path+"G.txt");
    outS.open(path+"S.txt");
    for(int u=0;u<1;u++)
    {
        readembed(u);
    }

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            G[i][j]/=2000;
            standG[i][j]/=2000;
            outG<<G[i][j]<<" ";
            outS<<standG[i][j]<<" ";
        }
        outG<<endl;
        outS<<endl;
    }
    outG.close();
    outS.close();

}

void RandG(string dataset, int bign)
{
    ifstream Gfile;
    ofstream   outG, outS;

    path="C:/Users/liang/PycharmProjects/whynotquery/BPR/500/";
    path+=dataset+"/";
    Gfile.open("G_t.txt");

    outG.open(path+to_string(bign)+"G.txt");

    for(int i=0;i<100;i++)
    {
        for(int j=0;j<100;j++)
        {
            Gfile>>G[i][j];
        }
    }
    for(int i=0;i<bign;i++)
    {
        for(int j=0;j<bign;j++)
        {
            G[i][j]=G[i%100][j%100]*1.5;
            outG<<G[i][j]<<" ";


        }
        outG<<endl;

    }
    outG.close();

}
void BPR()
{
    path="C:/Users/liang/PycharmProjects/whynotquery/Wrap/Netflix/";
    for(int u=0;u<2;u++)
    {
        readembed(u);
    }

}
// Driver program to test above functions
int main()
{
    // Let us create a graph shown in the above example
//    ConstructFlow();
//    Gettop100();
//    BPR();
    string names[5]={"Amazon", "Foursquare", "Movielens", "Netflix","RateBeer"};


        RandG(names[0], 50);
    RandG(names[0], 100);
    RandG(names[0], 200);
    RandG(names[0], 500);
    RandG(names[0], 800);
    RandG(names[0], 1000);

//    readembed(100);
//    string caps[3]={"0.05","0.1","0.2"};
//    Scale();
//    for(int i=1;i<=3;i++)
//    {
//        for(int c=0;c<3;c++)
//        {
//            writedata(i, "low", caps[c]);
//            writedata(i, "mid", caps[c]);
//            writedata(i, "high",caps[c]);
//        }
//
//    }

    return 0;

}