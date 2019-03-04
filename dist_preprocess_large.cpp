#include <bits/stdc++.h>
using namespace std;

#define _CRT_SECURE_NO_DEPRECATE

typedef long long ll;
typedef vector<int> vi;
typedef pair<int, int> ii;
typedef pair<int, ii> iii;
typedef vector<ii> vii;
typedef set<int> si;
typedef map<string, int> msi;

#define REP(i, a, b) \
for (int i = int(a); i <= int(b); i++)
#define IREP(i, a, b)\
for(int i = int(a); i >= int(b); i--)
#define TRvi(c, it) \
for (vi::iterator it = (c).begin(); it != (c).end(); it++)
#define TRvii(c, it) \
for (vii::iterator it = (c).begin(); it != (c).end(); it++)
#define TRmsi(c, it) \
for (msi::iterator it = (c).begin(); it != (c).end(); it++)

#define INF 1000000000
#define CLR(a) memset(a, 0, sizeof(a))
#define RESET(a) memset(a, -1, sizeof(a))
#define MP(a,b) make_pair((int)a, (int)b)
#define MP2(a,b,c) make_pair((int)a, make_pair((int)b, (int)c))
#define pb(a) push_back((ii)a)
#define EPS 1e-6
#define char_to_int(c)((int)c-(int)'A')
#define MAX 500010

bool contracted[MAX];
int delNeighbors[MAX], orderPos[MAX], d_dist[MAX], hops[MAX];
vii AdjList[MAX], revAdjList[MAX];
int n, m, s, t;
priority_queue<ii, vector<ii>, greater<ii>> d_pq, imp_pq;

void dijkstra(int src, int Max){
    memset(d_dist, -1, sizeof(d_dist)); memset(hops, 0, sizeof(hops));
    d_dist[src] = 0; while(!d_pq.empty())d_pq.pop();
    d_pq.push(MP(0, src));
    int u, w;
    while(!d_pq.empty()){
        u = d_pq.top().second; w = d_pq.top().first; d_pq.pop();
        if(w > d_dist[u])continue;
        if(hops[u] > 3)continue;
        if(d_dist[u] > Max)return;
        TRvii(AdjList[u], v){
            if(contracted[v->first])continue;
            if((d_dist[v->first] > d_dist[u] + v->second) || d_dist[v->first] == -1){
                d_dist[v->first] = d_dist[u] + v->second;
                hops[v->first] = hops[u]+1;
                d_pq.push(MP(d_dist[v->first], v->first));
            }
        }
    }
}

void contractNode(int u){
    contracted[u] = 1;
    TRvii(AdjList[u], v){
        delNeighbors[v->first]++;
    }
    TRvii(revAdjList[u], v){
        delNeighbors[v->first]++;
    }
    int in_Max = 0, out_Max = 0;
    TRvii(AdjList[u], v){
        out_Max = max(out_Max, v->second);
    }
    TRvii(revAdjList[u], v){
        in_Max = max(in_Max, v->second);
    }
    int Max = in_Max+out_Max, totCost;
    TRvii(revAdjList[u], v){
        if(contracted[v->first])continue;
        dijkstra(v->first, Max);
        TRvii(AdjList[u], w){
            if(contracted[w->first])continue;
            totCost = v->second+w->second;
            if(d_dist[w->first] > totCost || d_dist[w->first] == -1){
                AdjList[v->first].pb(MP(w->first, totCost)); revAdjList[w->first].pb(MP(v->first, totCost));
            }
        }
    }
}

int bidijkstra(){
    priority_queue<ii, vii, greater<ii>> f_pq, r_pq;
    vi f_dist(n+1, INF), r_dist(n+1, INF); int estimate = INF;
    f_dist[s] = 0; r_dist[t] = 0;
    f_pq.push(MP(0, s)); r_pq.push(MP(0, t));
    int u,w;
    while(!f_pq.empty() || !r_pq.empty()){
        if(!f_pq.empty()){
            u = f_pq.top().second; w = f_pq.top().first; f_pq.pop();
            if(w > f_dist[u])continue;
            if(f_dist[u] <= estimate){
                TRvii(AdjList[u], v){
                    if((f_dist[v->first] > f_dist[u] + v->second) && (orderPos[v->first] > orderPos[u])){
                        f_dist[v->first] = f_dist[u] + v->second; f_pq.push(MP(f_dist[v->first], v->first));
                    }
                }
                estimate = min(estimate, f_dist[u]+r_dist[u]);
            }
        }
        if(!r_pq.empty()){
            u = r_pq.top().second; w = r_pq.top().first; r_pq.pop();
            if(w > r_dist[u])continue;
            if(r_dist[u] <= estimate){
                TRvii(revAdjList[u], v){
                    if((r_dist[v->first] > r_dist[u] + v->second) && (orderPos[v->first] > orderPos[u])){
                        r_dist[v->first] = r_dist[u] + v->second; r_pq.push(MP(r_dist[v->first], v->first));
                    }
                }
                estimate = min(estimate, f_dist[u]+r_dist[u]);
            }
        }
    }
    return estimate;
}

inline int reCompute(int i){
    return (14*((AdjList[i].size()*revAdjList[i].size())-AdjList[i].size()-revAdjList[i].size()) + 25*(AdjList[i].size()+revAdjList[i].size()) + 10*delNeighbors[i]);
}

void preprocess(){
    int imp;
    for(int i = 1; i <= n; i++){
        imp = 14*((AdjList[i].size()*revAdjList[i].size())-AdjList[i].size()-revAdjList[i].size()) + 25*(AdjList[i].size()+revAdjList[i].size());
        imp_pq.push(MP(imp, i));
    }
    int u, u_imp, t_imp, cnt = 0;
    while(!imp_pq.empty()){
        u_imp = imp_pq.top().first; u = imp_pq.top().second; imp_pq.pop();
        t_imp = reCompute(u);
        if(imp_pq.size() != 0 && t_imp > imp_pq.top().first){
            imp_pq.push(MP(t_imp, u)); continue;
        }
        orderPos[u] = cnt;
        contractNode(u);
        cnt++;
    }
}

void init(){
    memset(delNeighbors, 0, sizeof(delNeighbors)); memset(contracted, 0, sizeof(contracted)); memset(orderPos, 0, sizeof(orderPos));
    int a, b, w;
    REP(i,1,m){
        cin>>a>>b>>w;
        AdjList[a].pb(MP(b,w)); revAdjList[b].pb(MP(a,w));
    }
    preprocess();
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL); cout.tie(NULL);

    cin>>n>>m;
    init();
    cout<<"Ready"<<endl;
    int q;
    cin>>q;
    REP(i,1,q){
        cin>>s>>t;
        if(s== t){
            cout<<0<<endl; continue;
        }
        int ans = bidijkstra();
        if(ans < INF)cout<<ans<<endl;
        else cout<<-1<<endl;
    }
    return 0;
}
