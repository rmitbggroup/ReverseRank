#include <stdio.h>            /* C input/output                       */
#include <stdlib.h>           /* C standard library                   */
#include "glpk.h"     /* GNU GLPK linear/mixed integer solver */
#include <iostream>
#include <map>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <cstdio>
#include <ctime>
#include <vector>
#include<queue>
#include <random>


//const int N = 100;
const int N = 6; //dimension
int M = 14000; //items
const int K = 100; //users
const int R = 20;
double G[N][N];
double tg[N][N];
double T[N][N];
double x[N],y[N], w[N], nw[N];
double e[N][N];
using namespace std;
double total_time=0.0;
double ct_lp =0.0;
int lp_ham=0;
string xp  ;
string yp  ;
string wp  ;
double getorank()
{
    double sum=0.0;
    for(int i=0;i<N;i++)
    {
        sum+=1.0*w[i]*(x[i]-y[i]);
    }
    return sum;
}


//string xp="C:/Users/liang/CLionProjects/whynot/x"+to_string(N)+".txt";
//string yp = "C:/Users/liang/CLionProjects/whynot/y"+to_string(N)+".txt";
//
//string wp = "C:/Users/liang/CLionProjects/whynot/w"+to_string(N)+".txt";

string tp;
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
            infile >> T[i][j];
            T[i][j]*=1.2;
//            cout<<T[i][j];
        }
    }
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

int getxy(int a,int b)
{
    return (a-1)*N+b;
}

int ia[1+26000000], ja[1+26000000];
double ar[1+26000000], z, x1, x2, x3;

double total_cost=0;
//static void *tls = NULL;
/* NOTE: in a re-entrant version of the package this variable should be
 * placed in the Thread Local Storage (TLS) */
void savee()
{
    ofstream oe;
    oe.open("e");
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            oe<<e[i][j]<<" ";
        }
        oe <<endl;
    }
}
double u[N];

int limit[N];
bool LP(int pre)
{

//    cerr<<"In LP"<<endl;
    //int ob, int t_size, string xy, string cap
//    string path = "C:/Users/liang/CLionProjects/whynot/";
//     xp ="C:/Users/liang/CLionProjects/whynot/"+xy+"_x/x_"+to_string(ob);
//     yp = "C:/Users/liang/CLionProjects/whynot/"+xy+"_y/y_"+to_string(ob);
//     wp = "C:/Users/liang/CLionProjects/whynot/w/w";
//     tp="C:/Users/liang/CLionProjects/whynot/T/"+to_string(N)+"T_"+to_string(t_size)+"miu"+cap+".txt";
//    xp=path+"data/x"+to_string(N)+".txt";
//    yp=path+"data/y"+to_string(N)+".txt";
//    wp=path+"data/w"+to_string(N)+".txt";
//    tp = path+"data/"+to_string(N)+"T_2miu0.1.txt";


//    cerr<<tp<<endl;
//    readT();
//    readxyw();
//    printxyw();
//    string outp="C:/Users/liang/CLionProjects/whynot/res/"+to_string(t_size)+xy+cap+"out.txt";
//    freopen(outp.c_str(),"w",stdout);
    lp_ham=0;
    clock_t start;
    memset(e,0,sizeof(e));
    double duration;

    start = std::clock();

    /* Your algorithm here */


    /* declare variables */

    /* create problem */
    glp_prob *lp;
    s1: lp = glp_create_prob();
    s2: glp_set_prob_name(lp, "LP");
    s3: glp_set_obj_dir(lp, GLP_MIN);

    s4: glp_add_rows(lp, 3*N+1);

    // set delta_w[i] * N
    for(int i=1;i<=N;i++)
    {
        glp_set_row_name(lp, i, ("dw"+to_string(i)).c_str() );
        if(limit[i-1]==1)
        {
            glp_set_row_bnds(lp, i, GLP_DB, 0, 0);
        }
        else
        {
            glp_set_row_bnds(lp, i, GLP_DB, -1.0 - w[i-1],1.0 - w[i-1]);
        }

    }

    // set target
    glp_set_row_name(lp, N+1, "target");
    glp_set_row_bnds(lp, N+1, GLP_UP, 0.0, -1.0*gettarget());

//    cerr<<" finish set w and target"<<endl;

    //Set u * 2N
    for(int i=1;i<=N;i++)
    {
        glp_set_row_name(lp, N+1+i, ("u+"+to_string(i)).c_str() );
        glp_set_row_bnds(lp, N+1+i, GLP_UP,0.0,0.0);
    }
    for(int i=1;i<=N;i++)
    {
        glp_set_row_name(lp, 2*N+1+i, ("u-"+to_string(i)).c_str() );
        glp_set_row_bnds(lp, 2*N+1+i, GLP_UP,0.0,0.0);
    }
//    cerr<<" finish set u"<<endl;
//
//    s5: glp_set_row_name(lp, 1, "p");
//    s6: glp_set_row_bnds(lp, 1, GLP_UP, 0.0, 100.0);
//    s7: glp_set_row_name(lp, 2, "q");
//    s8: glp_set_row_bnds(lp, 2, GLP_UP, 0.0, 600.0);
//    s9: glp_set_row_name(lp, 3, "r");
//    s10: glp_set_row_bnds(lp, 3, GLP_UP, 0.0, 300.0);

    glp_add_cols(lp, N*N+N);
//    s11: glp_add_cols(lp, 3);

    //set  0<=e_i,j<T_i,j
    for(int i=1;i<=N;i++)
    {
        for(int j=1;j<=N;j++)
        {
            glp_set_col_name(lp, getxy(i,j), ("e"+to_string(i)+","+to_string(j)).c_str());
            glp_set_col_bnds(lp, getxy(i,j), GLP_DB, 0.0, T[i-1][j-1]+1e-10);
        }
    }

    //add u >=0
    for(int i=1;i<=N;i++)
    {
        glp_set_col_name(lp, N*N+i, ("u"+to_string(i)).c_str());
        glp_set_col_bnds(lp, N*N+i, GLP_LO, 0.0, 0.0);

    }


    // set the  e ij no use
    for(int i=1;i<=N*N;i++)
    {
        glp_set_obj_coef(lp, i, 0.0);
    }
    // set the sum_u min
    for(int i=1;i<=N;i++)
    {
        glp_set_obj_coef(lp, N*N+i, 1.0);
    }


    // set -1-wi<= wi' - wi <=1-wi
    int ct=1;
    for(int i=1;i<=N;i++)
    {
        // wi

        for(int j=1;j<=N;j++)
        { //e [j,k]
            for(int k=1;k<=N;k++)
            {

                ia[ct] = i;
                ja[ct] = getxy(j,k);
                if(j==k)
                    ar[ct]=0;
                else if(k==i) //inflow
                    ar[ct] = 1.0; //+ej,i
                else if(j==i) //outflow
                    ar[ct] = -1.0;// - ei,k
                else
                    ar[ct] = 0.0;
                ct++;
            }
        }
        // u is not used
        for(int j=1;j<=N;j++)
        {
            ia[ct] = i;
            ja[ct] = N*N+j;
            ar[ct] = 0.0;
            ct++;
        }
    }

    // set target N+1
    for(int i=1;i<=N;i++)
    {
        for(int j=1;j<=N;j++)
        {

            ia[ct] = N+1;
            ja[ct] = getxy(i,j);
            ar[ct] = (x[j-1] - y[j-1]) -  (x[i-1] - y[i-1]);
            ct++;
        }

    }
    //set u , not used
    for(int j=1;j<=N;j++)
    {

        ia[ct] = N+1;
        ja[ct] = N*N+j;
        ar[ct] = 0;
        ct++;
    }

    //set u
    for(int i=1;i<=N;i++)
    {
        for(int j=1;j<=N;j++)
        {
            for(int k=1;k<=N;k++)
            {
                ia[ct] = i+N+1;
                ja[ct] = getxy(j,k);
                if(j==k)
                    ar[ct]=0.0;
                else if(k==i)
                    ar[ct] = 1.0;
                else if(j==i)
                    ar[ct] = -1.0;
                else
                    ar[ct] = 0.0;
                ct++;
            }
        }

        for(int j=1;j<=N;j++)
        {
            ia[ct] = i+N+1;
            ja[ct] = N*N+j;
            if(j==i)
                ar[ct] = -1.0;
            else
                ar[ct] = 0.0;
            ct++;
        }
    }

    //set anather u
    for(int i=1;i<=N;i++)
    {
        for(int j=1;j<=N;j++)
        {
            for(int k=1;k<=N;k++)
            {
                ia[ct] = i+2*N+1;
                ja[ct] = getxy(j,k);
                if(j==k)
                    ar[ct] = 0.0;
                else if(k==i)
                    ar[ct] = -1.0;
                else if(j==i)
                    ar[ct] = 1.0;
                else
                    ar[ct] = 0.0;
                ct++;
            }
        }
        for(int j=1;j<=N;j++)
        {
            ia[ct] = i+2*N+1;
            ja[ct] = N*N+j;
            if(j==i)
                ar[ct] = -1.0;
            else
                ar[ct] = 0.0;
            ct++;
        }
    }
    cout<<ct<<endl;
    s30: glp_load_matrix(lp, ct - 1, ia, ja, ar);
    s31: glp_simplex(lp, NULL);
    int state = glp_get_status(lp);
    if(state==GLP_NOFEAS)
    {
        return false;
    }
//    s32: z = glp_get_obj_val(lp);
//    s33: x1 = glp_get_col_prim(lp, 1);
//    s34: x2 = glp_get_col_prim(lp, 2);
//    s35: x3 = glp_get_col_prim(lp, 3);


    double minchange=1000.0;
    int index_min = 0;
    for(int i=1;i<=N;i++)
    {
        u[i-1] = glp_get_col_prim(lp, N*N+i);
        limit[i-1] = 1;
        if(u[i-1]>1e-10) {
            printf("u[%d] = %lf\n", i, u[i-1]);
            total_cost+=u[i-1];
            lp_ham+=1;
            limit[i-1]=0;
            if(u[i-1] <minchange)
            {
                minchange = u[i-1];
                index_min= i-1;
            }

        }

    }
    cerr<<lp_ham<<endl;

    limit[index_min] = 1;


    for(int i=1;i<=N;i++)
    {
//        nw[i-1] = w[i-1];
        for(int j=1;j<=N;j++) {
            double tt = glp_get_col_prim(lp, getxy(i, j));
            if (tt > 1e-8) {
                e[i-1][j-1]=tt;
                G[i-1][j-1]+=tt;
//                printf("e[%d][%d] = %lf\n ", i, j, tt);
                if(pre==0)
                {
                    w[i - 1] -= tt;
                    w[j - 1] += tt;
                }

            }
        }
    }
//    cerr<<glp_get_row_name(lp,N+1)<<endl;
//    cerr<<glp_get_row_stat(lp,N+1)<<endl;
    //    printxyw();

//    s36: printf("\nz = %g; x1 = %g; x2 = %g; x3 = %g\n", z, x1, x2, x3);

    /* housekeeping */

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
//    cerr<<duration<<endl;
    total_time+=duration;
    ct_lp+=1.0;
    cout<<"Time: "<< duration <<'\n';
    glp_delete_prob(lp);
    glp_free_env();

    cerr<<"Finished Deleting"<<endl;
    return true;

}


void writedata(int ts, string rank, string cap)
{
    ofstream myfile, yourfile;
    string p ="C:/Users/liang/PycharmProjects/whynotquery/Experiment Results/TS="+to_string(ts)+"/cap="+cap+"/"+rank+"/time";
    string costp ="C:/Users/liang/PycharmProjects/whynotquery/Experiment Results/TS="+to_string(ts)+"/cap="+cap+"/"+rank+"/cost";

    myfile.open (p);
    yourfile.open(costp);
    cerr<<p<<endl;
    ct_lp=0.0;
    total_time=0.0;
    total_cost=0;
    for(int i=0;i<50;i++)
    {
        if(i==43)
            continue;
//        LP(i,ts, rank,  cap);
    }
    cerr<<total_time/ct_lp<<endl;
    cerr<<total_cost/ct_lp<<endl;
    myfile << total_time/ct_lp;
    yourfile<<total_cost/ct_lp;
    myfile.close();
    yourfile.close();
}

double item[15000][N];
double wg[N], aw[N];

double getnrank()
{
    double sum=0.0;
    for(int i=0;i<N;i++)
    {
        sum+=1.0*wg[i]*(x[i]-y[i]);
    }
    return sum;
}
void updatew()
{
    for(int i=0;i<N;i++)
    {
        w[i]=nw[i];

    }
}


struct something{
    int id;
    double score;
}rg[15000],ri[15000];

bool compare(something t1, something t2)
{
    return t1.score>t2.score;
}

map<int,int> topset;
double  Measure(int k)
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
    sort(rg,rg+M, compare);
    sort(ri,ri+M, compare);

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

double  TopkDiff(int k)
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
            rg[i].score+=item[i][j]*aw[j];
            ri[i].score+=item[i][j]*w[j];
        }
    }
    sort(rg,rg+M, compare);
    sort(ri,ri+M, compare);

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



void readembed(int user_id)
{
    ifstream infileaw,infileitem,infilew, inpair;
    string path="C:/Users/liang/PycharmProjects/whynotquery/Wrap/Netflix/";
    tp = path+"G_t.txt";
    string awp=path+"gl";
    string itemp=path+"item";
    string wgp = path+"100user";
//    inpair.open(path+"pairs");
    infileaw.open (awp);
    infileitem.open (itemp);
    infilew.open (wgp);
    ofstream myfile, hitsfile, hamfile,oute, difffile;
    myfile.open (path+"/lp/lpwre"+to_string(user_id));
    hitsfile.open(path+"/lp/lphits"+to_string(user_id));
    difffile.open(path+"/lp/lpdiff"+to_string(user_id));
    hamfile.open(path+"/lp/lpham"+to_string(user_id));
    readT();

    memset(G, 0, sizeof(G));
    oute.open("e.txt");


    for(int i=0;i<N;i++)
    {
        infileaw>>aw[i];
//        w[i]=-0.5;
        w[i]=aw[i];

    }
    for(int i=0;i<=user_id;i++)
    {
        for(int j=0;j<N;j++)
        {
            infilew>>wg[j];
        }
    }
    for(int i=0;i<M;i++)
    {
        for(int j=0;j<N;j++)
        {
            infileitem >> item[i][j];
        }
    }
//    srand(time(0));
    int cx[N], cy[N], candidate=0;
    for(int step=1;step<=R;step++) {
        cerr<<"Step:"<<step<<endl;
        int ct = 0;
        candidate = 0;
        for (int loop = 1; loop <= 1; loop++) {

//                for(int i=0;i<N;i++)
//                {
//                    x[i]=item[a[ct]][i];
//                    y[i]=item[b[ct]][i];
//                }

            ct++;
//            cerr<<getorank() * getnrank()<<endl;
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
            if (getorank() < 0) {
                for (int i = 0; i < N; i++) {
                    swap(x[i], y[i]);
                }
            }
//            cerr<<gettarget()<<endl;
            cerr<<getorank()<<endl;
            cerr<<getnrank()<<endl;
               if(!LP(0))
               {
                   loop--;
               }


        }
//        cerr<<gettarget()<<endl;
        cerr<<getorank()<<endl;
        hitsfile<<Measure(K)<<endl;
        difffile<<TopkDiff(K)<<endl;
        hamfile<<lp_ham<<endl;
        // /        cout<<"count: "<<ct<<endl<<endl;
        for(int i=0;i<N;i++)
        {
            myfile<<w[i]<<" ";
        }

        myfile<<endl;
    }
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            G[i][j]/=100.0;
            oute<<G[i][j]<<" ";
        }
        oute<<endl;
    }

}


int inset[N*2], outset[N*2];



struct Node
{
    double cap, flow, cost;
};

Node edge[N+2][N+2];

double sum;
int size = 7;
int s = 0, t = 6;
int pre[N*2];
int q[N*4];
int vis[N*2];
double dis[N*2];
double flowdis[N*2];
double flow_sum = 0;
double outsum[N*2];
double insum[N*2];
double inweight[N*2];
double outweight[N*2];


void init(int source, int des, int size_n)
{
    size=size_n;

    int i, j;
    for(i = 0; i < size; i++)
        for(j = 0; j < size; j++)
        {
            edge[i][j].cost = 100000;
            edge[i][j].flow = 0;
            edge[i][j].cap = 0;
        }

    memset(insum,0,sizeof(insum));
    memset(outsum,0,sizeof(outsum));
    memset(inweight,0,sizeof(inweight));
    memset(outweight,0,sizeof(outweight));


    s=source, t=des;
    flow_sum=0;
}

double spfa()
{
    int i, temp;
    for(i = 0; i < size; i++)
    {
        vis[i] = -1;
        dis[i] = 10000000;
    }

    int front = 0, tail= 1;
    vis[s] = 1; dis[s] = 0;
    q[front] = s;
    while(front < tail)
    {

        temp = q[front];
        vis[temp] = -1;
        front++;
        for(i = 0; i < size; i++)
        {
            double res = edge[temp][i].cap - edge[temp][i].flow;
            if(res > 0 && dis[i] > dis[temp] + edge[temp][i].cost )
            {
                dis[i] = dis[temp] + edge[temp][i].cost;
                pre[i] = temp;
                if(vis[i] == -1)
                {
                    vis[i] = 1;
                    q[tail++] = i;
                }
            }
        }
    }
    return dis[t];
}


double goal = 0;
double MCMF()
{
    int temp;
    double  flow, sum = 0, cost, val;


    while(spfa() < 8000)
    {
        temp = t;
        flow = 100000; cost = 0;
        int tar_di=-1, sou_di=-1;
        while(temp != s)
        {
            val = max(edge[pre[temp]][temp].cap - edge[pre[temp]][temp].flow, 0.0);
            cost += edge[pre[temp]][temp].cost;
            if(val < flow)
                flow = val;
            temp = pre[temp];

            if(pre[temp]==s)
                tar_di = temp;
        }

        temp = t;
//        cout << flow << " " << cost << endl;
        sum += flow * cost;
        flow_sum += flow;
        while(temp != s)
        {
            if(temp==t)
                sou_di = pre[temp];
            edge[pre[temp]][temp].flow += flow;
            edge[temp][pre[temp]].flow -= flow;
//			cout << pre[temp] << " " << temp << " ";
//			cout << edge[pre[temp]][temp].flow << endl;
            temp = pre[temp];
        }
//        cerr<<sou_di<<"->"<<tar_di<<":"<<flow*cost - flow*2*2<<endl;
        outsum[sou_di] += flow*cost - flow*2*2;
        insum[tar_di] += flow*cost - flow*2*2;

        inweight[tar_di] += flow;
        outweight[sou_di] += flow;

        tg[sou_di][tar_di]+=abs(flow*cost - flow*2*2);
        e[sou_di][tar_di]+=1.0;
    }
    cerr<<"MCMF"<<endl;
    return sum;

}


double HeuristicMCMF()
{
    init(N,N+1,N+2);
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            if(i!=j)
            {
                edge[i][j].cap = T[i][j];
//                cout<<"T[i][j]:"<<T[i][j]<<endl;
                edge[i][j].cost = 0;
            }
        }
    }

//    cerr<<"Inset:";
    for(int i=0;i<N;i++)
    {
        if(inset[i])
        {
//            cerr<<i<<" ";
            //range -1, 1
            //edge[N][i].cap = max(min(1 + w[i],2.0),0.0);
            //range 0, 1
            edge[N][i].cap = max(min( w[i], 1.0),0.0);
//            cout<<edge[N][i].cap<<endl;
            //x[i] > y[i]
            edge[N][i].cost= (y[i] - x[i] + 2);
//            cout<<edge[N][i].cost<<endl;
        }

    }
//    cerr<<endl;

//    cerr<<"Outset: ";
    for(int i=0;i<N;i++)
    {
        if(outset[i])
        {
//            cerr<<i<<" ";
//            cout<<"out:"<<i<<" ";
            //range -1, 1
//            edge[i][N+1].cap = max(min(1 - w[i], 2.0),0.0);
            // range 0, 1
            edge[i][N+1].cap = max(min(1- w[i], 1.0),0.0);
//            cout<<edge[i][N+1].cap<<endl;
            //y[i] > x[i]
            edge[i][N+1].cost= (x[i] - y[i] + 2);
//            cout<<edge[i][N+1].cost<<endl;
        }
    }
//    cerr<<endl;

   double cla = MCMF();
    int flag = 1;
    double res = cla - 2*2*flow_sum;
//    cerr<<cla<<" "<<flow_sum<<endl;
//    cerr<<"Res Before:"<<res<<endl;
    if(res<0)
    {
//        cerr<<"Insum"<<endl;

        for(int i=0;i<N;i++)
        {

            if(abs(insum[i])<1e-10&&inset[i])
            {
                inset[i]=0;
                insum[i]=0;
//                cout<<"remove Inset"<<i<<endl;
            }

            if(goal+res-insum[i]<=1e-10&&inset[i])
            {

//                cerr<<"Goal:"<<goal<<" & res "<<res<<endl;
//                cerr<<i<<":"<<insum[i];
//                cerr<<" "<<inset[i]<<endl;

                flag=0;
                res -= insum[i];
                inset[i]=0;
                insum[i]=0;

//                cout<<"remove Inset"<<i<<endl;
            }



        }

        if(flag)
        {
            for(int i=0;i<N;i++)
            {

                if(abs(outsum[i])<1e-10&&outset[i])
                {
                    outset[i]=0;
                    outsum[i]=0;
//                cout<<"remove Outset"<<i<<endl;
                }

                if(goal+res-outsum[i]<=1e-10&&outset[i])
                {
//                    cerr<<i<<":"<<outsum[i]<<" ";
//                    cerr<<outset[i]<<endl;

                    flag=0;
                    res-=outsum[i];
                    outset[i]=0;
                    outsum[i]=0;

//                cout<<"remove Outset"<<i<<endl;
                }

            }
        }

    }
//    cerr<<res<<endl;


    if(flag==0)
    {
        return HeuristicMCMF();
    }
//    double ingoal=0;
//    for(int i=0;i<N;i++)
//    {
//        if(inset[i]||outset[i])
//        {
//            cout<<i<<" "<<insum[i]<<" "<<outsum[i]<<endl;
//            ingoal+=insum[i];
//        }
//    }
//    cout<<"Goal:"<<goal<<endl;
//    cout<<"Ingoal"<<ingoal<<endl;
    cout<<"Res: "<<res<<endl;
    return res;

}


struct ccc
{
    int from, to;
    double weights;
};

bool comp(ccc a, ccc b)
{
    return a.weights>b.weights;
}
void naive()
{
    ccc ar[N*N];
    int counter=0;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            ar[counter].from=i;
            ar[counter].to=j;
            double mxf=w[i]+1.0;
            double mxt=1.0 - w[j];
            ar[counter].weights=(x[i]-y[i] + y[j] - x[j])*min(mxf, mxt);
            counter++;
        }
    }
    sort(ar, ar+N*N, comp);
    double goal=0.0;
    for(int i=0;i<N;i++)
    {
        goal+=w[i]*(x[i] - y[i]);
    }

    counter=0;
    while(goal>0)
    {
        int from=ar[counter].from;
        int to=ar[counter].to;
        double wt=ar[counter].weights;
        double eff= (x[from]-y[from] + y[to] - x[to]);
        if(wt<goal)
        {
            double flow= wt/(eff);
            w[from]-=flow;
            w[to]+=flow;
            goal-=wt;
        }
        else
        {
            double flow = goal/eff;
            w[from]-=flow;
            w[to]+=flow;
            goal=0.0;
        }

    }

}

int testmc(int pre_lp)
{
//    loaddata();
    flow_sum=0;

//    for(int i=0;i<N;i++)
//    {
////        x[i]=item[1][i];
////        y[i]=item[2][i];
//        goal += (x[i]-y[i])*w[i];
//    }

    clock_t start;
    double duration;

    start = std::clock();
    memset(inset,0,sizeof(inset));
    memset(outset,0,sizeof(outset));



    if(pre_lp==0)
    {
        for(int i=0;i<N;i++)
        {
            if(x[i]>y[i])
            {
                inset[i]=1;
            }
            else
            {
                outset[i]=1;
            }
        }
    }
    else
    {
        for(int i=0;i<N;i++)
        {
            if(x[i]>y[i]&&u[i]>1e-8)
            {
                inset[i]=1;
            }
            else if(x[i]<y[i]&&u[i]>1e-8)
            {
                outset[i]=1;
            }
        }
    }




    double cla = HeuristicMCMF();
//    cout<<cla<<endl;
    int Ham=0;
    double rate=1.0;
    double sum_xx=0;
    for(int i=0;i<N;i++)
    {
        if(inset[i]||outset[i])
        {
            Ham++;
            sum_xx+=abs(insum[i]);
//            cout<<i<<" "<<insum[i]<<" "<<outsum[i]<<endl;
        }
    }
    rate= abs(goal)/abs(sum_xx);
//    cerr<<"Rate:"<<rate<<endl;
    rate=min(rate, 1.0);
    for(int i=0;i<N;i++)
    {
        if(inset[i]||outset[i])
        {

            cerr<<i<<": -"<<inweight[i]*rate<<"| +"<<outweight[i]*rate<<endl;

            w[i] -= inweight[i]*rate;
            w[i] += outweight[i]*rate;
            //change it to above zero
            w[i] = max(0.0, min(w[i],1.0));

//            cout<<i<<": "<<w[i]<<endl;
//            w[i] = max(-1.0,min(w[i],1.0));
//            cout<<i<<" "<<insum[i]<<" "<<outsum[i]<<endl;
        }
    }
    duration = std::clock() - start;

//    cout<<"Time: " <<duration/CLOCKS_PER_SEC<<endl;
    return Ham;
}

double wt[105][N];

int pa,pb;
string path;
ofstream curtimefile;

void VaryTest(int g,int user_id, int rounds, int method, int g_size, int ue)
{
    ifstream infileaw,infileitem,infilew,pairspath;
    string name[5]={"lp", "mcmf", "lpmc","naive","iterlp"};



    tp = path+"G.txt";

    string awp=path+"gl";
    string itemp=path+"item";
    string wgp = path+"100user";
//    std::random_device rd;
    mt19937 mt_rand(time(0));


    ofstream hamcompare;
//    inpair.open(path+"pairs");
    infileaw.open (awp);
    infileitem.open (itemp);
    infilew.open (wgp);
    ofstream myfile, hitsfile, hamfile, timefile, difffile;
    timefile.open(path+name[method-1]+"/time"+to_string(user_id));
    myfile.open (path+name[method-1]+"/wres"+to_string(user_id));
    if(ue==1)
    {
        hitsfile.open(path+name[method-1]+"/hits"+to_string(user_id));
        difffile.open(path+name[method-1]+"/diff"+to_string(user_id));
    }

    hamfile.open(path+name[method-1]+"/"+to_string(g)+"ham"+to_string(user_id));
    hamcompare.open(path+name[method-1]+"/"+to_string(g)+"fea"+to_string(g_size)+to_string(user_id));
    readT();

    for(int i=0;i<N;i++)
    {
        infileaw>>aw[i];

        w[i]=aw[i];

    }
    for(int i=0;i<=user_id;i++)
    {
        for(int j=0;j<N;j++)
        {
            infilew>>wg[j];
        }
    }

    double sum_gl=0.0;
    for(int i=0;i<N;i++)
    {
        sum_gl+=wg[i];
    }
    for(int i=0;i<N;i++)
    {
        aw[i]=sum_gl/N;
        w[i]=aw[i];
    }
    Measure(K);
    for(int i=0;i<M;i++)
    {
        for(int j=0;j<N;j++)
        {
            infileitem >> item[i][j];
        }
    }


    int cx[N], cy[N], candidate=0;
    int ham=0;
    int nofes=0;
    double ttham=0;
    double duration=0;
    clock_t start = clock();
    for(int step=1;step<=rounds;step++) {
        cerr<<"Step:"<<step<<endl;
        int ct = 0;
        candidate = 0;
        start = clock();
        for (int loop = 1; loop <= 1; loop++) {


            ct++;
            while (1) {
                int ta = mt_rand() % M;
                if(topset[ta]==1)
                    continue;
                int tb =mt_rand() % K;
                for (int i = 0; i < N; i++) {
                    x[i] = item[ta][i];
                    y[i] = item[rg[tb].id][i];
                }
                pb = rg[tb].id;
                pa=ta;
                if (getorank() * getnrank() < 0) {
                    break;
                }


            }
            if (getorank() < 0) {
                for (int i = 0; i < N; i++) {
                    swap(x[i], y[i]);
                }
                swap(pa,pb);
            }
            cerr<<pa<<" "<<pb<<endl;
//            cerr<<gettarget()<<endl;
            cerr<<getorank()<<endl;
//            cerr<<getnrank()<<endl;
            goal = getorank();

            for(int i=0;i<N;i++)
            {
                aw[i]=w[i];
                limit[i]=0;
            }
            if(method==1)
            { //linear programming
                LP(0);
                ham = lp_ham;
            }
            else if (method==2)
            { //mincostmaxflow  itrssp
                ham = testmc(0);
            }
            else if(method==3)
            { //linear programming + itrssp
                LP(1);
                cerr<<lp_ham<<endl;
                ham = testmc(1);
                cerr<<ham<<endl;
            }else  if(method==4)
            {  // greedy algorithm
                naive();
            }
            else
            {  //updated lp
                for (int i = 0; i < N; i++)
                {
                    limit[i] = 0;
                    nw[i] = w[i];
                }
                while(1)
                {
                    updatew();
                    bool feasibile = LP(0);
                    if(!feasibile)
                        break;
                }


            }

            if(getorank()>1e-5)
                nofes++;
            else
            {
                ttham+=ham;
                hamfile<<ham<<endl;

                curtimefile<<"method"<<method<<"gcap"<<g<<"gsize"<<g_size<<"time"<<(clock() -start)/CLOCKS_PER_SEC<<endl;
            }

        }
        double ingoal=0;
        for(int i=0;i<N;i++)
        {
            if(inset[i]||outset[i])
            {
                cout<<i<<" "<<insum[i]<<" "<<outsum[i]<<endl;
                ingoal+=insum[i];
            }
        }
        cout<<"Goal:"<<goal<<endl;
        cout<<"Ingoal"<<ingoal<<endl;
//        cerr<<gettarget()<<endl;
        cerr<<getorank()<<endl;
        if(ue==1)
        {
            hitsfile<<Measure(K)<<endl;
            difffile<<TopkDiff(K)<<endl;
        }
//       } cout<<"count: "<<ct<<endl<<endl;


        for(int i=0;i<N;i++)
        {
            myfile<<w[i]<<" ";
        }
        myfile<<endl;
    }

    hamcompare<<"G-cap"<<g<<"Ham"<<name[method-1]<<" "<<ttham/(rounds-nofes+1)<<endl;
    cerr<<"Ham"<<name[method-1]<<" "<<ttham/(rounds-nofes+1)<<endl;
    hamcompare<<"G-cap"<<g<<"Name"<<name[method-1]<< g_size<<" no-fea rate"<<nofes*1.0/rounds<<endl;
    cerr<<"No_feas_rate"<<nofes*1.0/rounds<<endl;
    timefile<<(clock() -start)/CLOCKS_PER_SEC<<endl;
}


void mcreadembed(int user_id)
{
    ifstream infileaw,infileitem,infilew, filewt;
    string path="C:/Users/liang/PycharmProjects/whynotquery/Wrap/Netflix/";
    tp = path+"1G";
    string awp=path+"gl";
    string itemp=path+"item";
    string wgp = path+"100user";
    string costp=path+"/mcmf/mdaham";
    filewt.open(wgp);
//    inpair.open(path+"pairs");
    infileaw.open (awp);
    infileitem.open (itemp);
    infilew.open (wgp);
    ofstream myfile, hitsfile, hamfile, difffile;
    myfile.open (path+"/mcmf/mdawre"+to_string(user_id));
    hitsfile.open(path+"/mcmf/mdahits"+to_string(user_id));
    difffile.open(path+"/mcmf/mdadiff"+to_string(user_id));
    hamfile.open(costp+to_string(user_id));
    readT();


    for(int i=0;i<N;i++)
    {
        infileaw>>aw[i];
//        w[i]=-0.5;
        w[i]=aw[i];

    }
    for(int i=0;i<100;i++)
    {
        for(int j=0;j<N;j++)
        {
            filewt>>wt[i][j];
        }
    }
    for(int i=0;i<=user_id;i++)
    {
        for(int j=0;j<N;j++)
        {
            infilew>>wg[j];
        }
    }
    for(int i=0;i<M;i++)
    {
        for(int j=0;j<N;j++)
        {
            infileitem >> item[i][j];
        }
    }

//    srand(time(0));
    int cx[N], cy[N], candidate=0;
    int ham=0;
    for(int step=1;step<=R;step++) {
        cerr<<"Step:"<<step<<endl;
        int ct = 0;
        candidate = 0;
        for (int loop = 1; loop <= 1; loop++) {

//                for(int i=0;i<N;i++)
//                {
//                    x[i]=item[a[ct]][i];
//                    y[i]=item[b[ct]][i];
//                }

            ct++;
//            cerr<<getorank() * getnrank()<<endl;
            while (1) {
//                pause(1);

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
            if (getorank() < 0) {
                for (int i = 0; i < N; i++) {
                    swap(x[i], y[i]);
                }
            }
//            cerr<<gettarget()<<endl;
            cerr<<getorank()<<endl;
//            cerr<<getnrank()<<endl;
            ham = testmc(0);
            if(getorank()>1e-5)
                loop--;

        }
        double ingoal=0;
        for(int i=0;i<N;i++)
        {
            if(inset[i]||outset[i])
            {
                cout<<i<<" "<<insum[i]<<" "<<outsum[i]<<endl;
                ingoal+=insum[i];
            }
        }
        cout<<"Goal:"<<goal<<endl;
        cout<<"Ingoal"<<ingoal<<endl;
//        cerr<<gettarget()<<endl;
        cerr<<getorank()<<endl;
        hitsfile<<Measure(K)<<endl;
        difffile<<TopkDiff(K)<<endl;
//        cout<<"count: "<<ct<<endl<<endl;

        hamfile<<ham<<endl;
        for(int i=0;i<N;i++)
        {
            myfile<<w[i]<<" ";
        }
        myfile<<endl;
    }
}


void testframe(string cur)
{

    string dataset[5]={"Amazon","Netflix", "Foursquare","Movielens", "RateBeer"};
    int msize[5]={14000, 14000,14000,3600,2700};
    for(int d = 0;d<3;d++)
    {
//        if(cur.size()==0)
//            path = "C:/Users/liang/PycharmProjects/whynotquery/Wrap/500/"+dataset[d]+"/";
//        else
//        {
//            path=cur;
//            curtimefile.open("timesaver.txt");
//            d=10;
//        }

//        path="C:/Users/liang/PycharmProjects/whynotquery/Wrap/"+dataset[d]+"/";
        M = msize[d];
        path="C:/Users/liang/PycharmProjects/whynotquery/BPR/"+dataset[d]+"/";

        for(int i=25;i<=30;i++)
        {
            int user_id=i;
            VaryTest(2, user_id, 50, 5, 2, 1);
//            VaryTest(2,user_id, 50,1,2,1);
//            VaryTest(2,user_id, 50,2,2,1);
//            VaryTest(2,user_id, 50,3,2,1);
//            VaryTest(2,user_id, 50,4,2,1);
//            VaryTest(2, user_id, 50, 5, 2,1);
        }

//        for(int vary_cap=0;vary_cap<6;vary_cap++)
//        {
//            VaryTest(vary_cap,user_id,100,1,2,0);
//            VaryTest(vary_cap,user_id,100,2,2,0);
//            VaryTest(vary_cap,user_id,100,3,2,0);
//        }
//        for(int vary_size=0;vary_size<5;vary_size++)
//        {
//            VaryTest(1,user_id,100,1,vary_size,0);
//            VaryTest(1,user_id,100,2,vary_size,0);
//            VaryTest(1,user_id,100,3,vary_size,0);
//        }




    }

//
//    for(int s=0;s<5;s++)
//    {
//        for(int g=0;g<7;g++)
//        {
////            VaryTest(g,0,100, 1,s);
////            VaryTest(g,0,100, 2,s);
////            VaryTest(g,0, 100,3, s);
//        }
//    }

}



string curpath= "C:/Users/liang/Downloads/housevalue.txt";
string logpath = "C:/Users/liang/Downloads/logs/log.txt";
string gpath = "C:/Users/liang/Downloads/logs/genG.txt";
vector<vector<double> > itemvec, rankitems;
ofstream logfile;

double uplimit[6] ={2101.0, 11.0, 1605.0, 128.0, 4309.0, 71.0};
vector<int> invert(6);
void printrank(int flag)
{
    if(flag>1)
    {
//        for(int i=0;i<6;i++)
        {
            cout<<"Price BedRooms  Landsize Education DistoTrain Transportation"<<endl;
        }
        logfile<<"10 Properties Rank"<<endl;
    }

    vector<pair<double, int>>rankscore;
    for(int i=0;i<10;i++)
    {
        rankscore.push_back({0,0});
        for(int j=0;j<6;j++)
        {
            rankscore[i].first += rankitems[i][j]*w[j];
        }

        rankscore[i].second  = i;


    }
    sort(rankscore.begin(), rankscore.end());
    for(int i=0;i<10;i++)
    {
        int id = rankscore[9-i].second;
        invert[i] = id;
        if(flag>1)
            cout<<i<<" ";
        if(flag%2)
        {
            logfile<<i<<" ";
        }
        for(int j=0;j<6;j++)
        {
            if(flag%2)
            {
//                logfile<<rankitems[id][j]<<" ";
                if(j==0)
                    logfile<<2001-rankitems[id][j]*uplimit[j]<<" ";
                else if(j==4)
                    logfile<<4500 - rankitems[id][j]*uplimit[j]<<" ";
                else if(j==5)
                    logfile<<80-rankitems[id][j]*uplimit[j]<<" ";
                else logfile<< rankitems[id][j]*uplimit[j]<<" ";
            }
            if(flag>1)
            {
                if(j==0)
                    cout<<2001-rankitems[id][j]*uplimit[j]<<" ";
                else if(j==4)
                    cout<<4500 - rankitems[id][j]*uplimit[j]<<" ";
                else if(j==5)
                    cout<<80-rankitems[id][j]*uplimit[j]<<" ";
                else cout<< rankitems[id][j]*uplimit[j]<<" ";
            }

        }
        if(flag%2)
            logfile<<endl;
        if(flag>1)
            cout<<endl;
    }
    if(flag>1)
    {
        for(int i=0;i<6;i++)
        {
            cout<<w[i]<<" ";
        }
        cout<<endl;
    }
    if(flag%2)
    {
        for(int i=0;i<6;i++)
        {
            logfile<<w[i]<<" ";
        }
        logfile<<endl;
    }


}
void inputall(int numItems)
{
    int preitem[10]={248, 337 ,271 ,357 ,516 ,191, 117, 338, 522, 82 };
    logfile.open(logpath,ofstream::app);
    logfile<<"#"<<endl;
    ifstream inputitems;
    inputitems.open(curpath);
    mt19937 mt_rand(time(0));

    for(int ct=0; ct<numItems; ct++)
    {
        vector<double> tempvec;
        for(int i=0;i<6;i++)
        {
            double tmp;
            inputitems >> tmp;
            tempvec.push_back(tmp);
//            cout<<tmp<<" ";
        }
//        cout<<endl;
        itemvec.push_back(tempvec);
    }
    for(int i=0;i<10;i++)
    {
//        int temp = mt_rand()%numItems;
        int temp = preitem[i];
        vector<double> tmp;
        logfile<<temp<<" ";
        for(int j=0;j<6;j++)
        {
            tmp.push_back(itemvec[temp][j]);
        }
        rankitems.push_back(tmp);
    }
    logfile<<endl;

    cout<<" Selected"<<endl;
    for(int i=0;i<10;i++)
    {
        cout<<i<<": ";
        for(int j=0;j<6;j++)
        {
            logfile<<rankitems[i][j]<<" ";
            cout<<rankitems[i][j]<<" ";
        }
        logfile<<endl;
        cout<<endl;
    }


}
void initW()
{
    for(int i=0;i<6;i++)
    {
        w[i]= 1.0/6;
    }
}

void loadG()
{
    ifstream inputG;
    inputG.open(gpath);
    logfile<<"The G"<<endl;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            inputG>>T[i][j];
            logfile<<T[i][j]<<" ";
        }
        logfile<<endl;
    }

}
void BuildG()
{
    // Price, Bedrooms, landsize, education, dis, transportation.
    logfile<<"The G"<<endl;
    for(int i=0; i<6; i++)
    {
        for(int j=0; j<6; j++)
        {

            T[i][j]=0.2;
            logfile<<T[i][j]<<" ";
            if(i==0||j==0)
                T[i][j] = 0.05;
//            if(i==5)
//                T[i][j] = 0.05;
////            if(i==1)
////                T[i][j] = 0.08;
//
//            if(j==5)
//                T[i][j] = 0.1;

        }
        logfile<<endl;
    }

}


void simulation()
{

    inputall(663);
    initW();
//    BuildG();
    loadG();
    while(1)
    {
        printrank(3);
        cout<<"Which pair do you want to reverse? (0 0 for break)"<<endl;
        int tmpx, tmpy;
        cin>>tmpx>>tmpy;
        logfile<<tmpx<<" "<<tmpy<<endl;
        if(tmpx==0&&tmpy==0)
            break;

        for(int i=0; i<6; i++)
        {
            int idx, idy;
            idx = invert[tmpx];
            idy = invert[tmpy];
            x[i] = rankitems[idx][i];
            y[i] = rankitems[idy][i];
        }
        cout<<2001 - x[0]*uplimit[0]<<endl;
        cout<<2001 - y[0]*uplimit[0]<<endl;

        if (getorank() < 0) {
            for (int i = 0; i < 6; i++) {
                swap(x[i], y[i]);
            }
        }
        //set goal before use MinCostMaxFlow
        goal = getorank();
        cout<<"Ham Dis: "<<testmc(0)<<endl;
    }
}

void GenerateG(int repTimes)
{
    const int numItems = 663;
    ifstream inputitems;
    inputitems.open(curpath);
    mt19937 mt_rand(time(0));
    ofstream Gpath;
    Gpath.open(gpath);


    for(int ct=0; ct<numItems; ct++)
    {
        vector<double> tempvec;
        for(int i=0;i<N;i++)
        {
            double tmp;
            inputitems >> tmp;
            tempvec.push_back(tmp);
        }

        itemvec.push_back(tempvec);
    }
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            T[i][j] = 1.0;
            tg[i][j] = 0.0;
            e[i][j] = 0.0;
        }
    }


    for(int c=1;c<=repTimes; c++)
    {
        int tx = mt_rand()%numItems;
        int ty = mt_rand()%numItems;
        initW();
        for(int i=0;i<N;i++)
        {
            x[i] = itemvec[tx][i];
            cerr<<y[i]<<" ";
            y[i] = itemvec[ty][i];
        }
        cerr<<endl;


        if (getorank() < 0) {
            for (int i = 0; i < N; i++) {
                swap(x[i], y[i]);
            }
        }
        //set goal before use MinCostMaxFlow
        goal = getorank();
        cerr<<"Goal:"<<goal<<endl;
        cout<<"Ham Dis: "<<testmc(0)<<endl;
    }

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            tg[i][j]/=e[i][j];

            if(e[i][j]>1e-6)
                Gpath<<tg[i][j]<<" ";
            else
                Gpath<<0.0<<" ";
        }
        Gpath<<endl;
    }





}
int main(int argc, char** argv)
{

    //testframe(path); 
    // this is for experiments.

    //    GenerateG(10000); 
    // this is to generate the transition graph

    //    simulation();
    // this is for simulation of the real world data
    //    logfile<<"#"<<endl;


    return 0;
}