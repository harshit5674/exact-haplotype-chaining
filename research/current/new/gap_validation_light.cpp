// Independent validation only, memory-light variant (2026-10-20) of gap_validation.cpp for graphs whose
// all-pairs distance table (V^2 words) does not fit, such as base-level windows
// (experiments/baselevel_feasibility). Only the index differs from gap_validation.cpp:
//  - last reachable vertices come from the definition: for each label, the path's vertices are taken
//    from last to first and each marks, by graph search, every vertex it reaches and no later path
//    vertex reached; a vertex reached from a later path vertex is skipped, with everything it reaches;
//  - remote distances come from one single-source DAG shortest-path pass per distinct last vertex,
//    stored only for the (vertex, label) pairs that use them: V*H words instead of V*V.
// The sweep, the checks and the C ABI are unchanged. No dependency on gap_native.cpp or its index.
// Graph-major sweep with Fenwick trees over READ endpoints (the engines sweep reads
// with trees over GRAPH endpoints). Never used for timings.
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
using I = int64_t;
static constexpr I NEG=-(I(1)<<61), INF=(I(1)<<61);
static thread_local std::string error;
static void require(bool b,const char* msg){if(!b)throw std::runtime_error(msg);}
struct Index {
 int V,H; std::vector<I> len,phi,dist; std::vector<int> last,rank;
 std::vector<std::vector<int>> labels,paths,adj;
 Index(int v,int h,int e,const I* l,const int* edges,const int* po,const int* ps):
 V(v),H(h),len(l,l+v),phi(size_t(v)*h,-1),dist(size_t(v)*h,INF),last(size_t(v)*h,-1),
 rank(size_t(v)*h,-1),labels(v),paths(h),adj(v){
  require(V>0&&H>0,"empty graph");
  for(int k=0;k<e;k++){int a=edges[2*k],b=edges[2*k+1];require(0<=a&&a<b&&b<V,"non-DAG edge");adj[a].push_back(b);}
  for(int g=0;g<H;g++){
   paths[g]=std::vector<int>(ps+po[g],ps+po[g+1]); I off=0; int prev=-1;
   require(!paths[g].empty(),"empty path");
   for(int x:paths[g]){
    require(x>=0&&x<V&&x>prev,"non-simple ordered path");
    if(prev>=0)require(std::find(adj[prev].begin(),adj[prev].end(),x)!=adj[prev].end(),"path edge absent");
    require(len[x]>0,"nonpositive length");phi[size_t(x)*H+g]=off;off+=len[x];
    rank[size_t(x)*H+g]=labels[x].size();labels[x].push_back(g);prev=x;
   }
  }
  // Definition-based lookup, independent of the production last/distance DP: last[x][g] is the last
  // vertex of P_g from which x is reachable. Later path vertices claim first; a claimed vertex was reached
  // from a later path vertex, and so was everything it reaches, so the search stops there.
  std::vector<int> stamp(V,-1),stack;
  for(int g=0;g<H;g++)for(auto it=paths[g].rbegin();it!=paths[g].rend();++it){
   int u=*it;if(stamp[u]==g)continue;
   stamp[u]=g;last[size_t(u)*H+g]=u;stack.assign(1,u);
   while(!stack.empty()){int a=stack.back();stack.pop_back();
    for(int b:adj[a])if(stamp[b]!=g){stamp[b]=g;last[size_t(b)*H+g]=u;stack.push_back(b);}}
  }
  // Remote distances D(last[x][g], x) for labels g missing at x: one single-source pass per distinct source.
  std::vector<std::vector<size_t>> need(V);
  for(int x=0;x<V;x++)for(int g=0;g<H;g++)
   if(rank[size_t(x)*H+g]<0&&last[size_t(x)*H+g]>=0)need[last[size_t(x)*H+g]].push_back(size_t(x)*H+g);
  std::vector<I> d(V,INF);
  for(int s=0;s<V;s++)if(!need[s].empty()){
   std::fill(d.begin()+s,d.end(),INF);d[s]=0;
   for(int a=s;a<V;a++)if(d[a]!=INF)for(int b:adj[a])d[b]=std::min(d[b],d[a]+len[a]);
   for(size_t k:need[s]){require(d[k/H]!=INF,"last vertex does not reach");dist[k]=d[k/H];}
  }
 }
};
struct Anchor {int v;I x,y,c,d,w;};
struct Event {int v;I p;int type,i,h;bool operator<(const Event& b)const{return std::tie(v,p,type,i,h)<std::tie(b.v,b.p,b.type,b.i,b.h);}};
struct Batch {
 Index& g;int n;I gamma;std::vector<Anchor>a;std::vector<int>off,owner;
 std::vector<I>C,B,R;
 Batch(Index& graph,int count,const I* flat,I penalty):g(graph),n(count),gamma(penalty),off(n+1),B(n),R(n,NEG){
  require(gamma>=0,"negative gamma");
  std::vector<I> ends;std::vector<Event> ev;
  for(int i=0;i<n;i++){
   const I* f=flat+6*i;Anchor z{int(f[0]),f[1],f[2],f[3],f[4],f[5]};
   require(z.v>=0&&z.v<g.V&&1<=z.x&&z.x<=z.y&&z.y<=g.len[z.v]&&1<=z.c&&z.c<=z.d,"invalid anchor");
   require(!g.labels[z.v].empty(),"anchor without biological labels");a.push_back(z);ends.push_back(z.d);
   off[i+1]=off[i]+g.labels[z.v].size();
   ev.push_back({z.v,z.x,1,i,-1});ev.push_back({z.v,z.x,2,i,-1});ev.push_back({z.v,z.y,3,i,-1});
   for(int h=0;h<g.H;h++)if(g.rank[size_t(z.v)*g.H+h]<0){
    int u=g.last[size_t(z.v)*g.H+h];if(u>=0)ev.push_back({u,g.len[u]+1,0,i,h});
   }
  }
  C.assign(off[n],NEG);owner.resize(off[n]);
  for(int i=0;i<n;i++)std::fill(owner.begin()+off[i],owner.begin()+off[i+1],i);
  std::sort(ends.begin(),ends.end());ends.erase(std::unique(ends.begin(),ends.end()),ends.end());
  std::sort(ev.begin(),ev.end());int stride=ends.size()+1;
  std::vector<I> trees(size_t(stride)*g.H,NEG);
  auto query=[&](int h,I c){int k=std::lower_bound(ends.begin(),ends.end(),c)-ends.begin();I best=NEG;
   for(;k;k-=k&-k)best=std::max(best,trees[size_t(h)*stride+k]);return best;};
  auto update=[&](int h,I d,I val){int k=std::lower_bound(ends.begin(),ends.end(),d)-ends.begin()+1;
   for(;k<stride;k+=k&-k)trees[size_t(h)*stride+k]=std::max(trees[size_t(h)*stride+k],val);};
  for(const auto& e:ev){int i=e.i;const auto& z=a[i];I alpha=z.w-z.x-z.c+2;
   if(e.type==0){I q=query(e.h,z.c);if(q!=NEG){I val=q+alpha-g.phi[size_t(e.v)*g.H+e.h]-g.dist[size_t(z.v)*g.H+e.h]-gamma;R[i]=std::max(R[i],val);}}
   else if(e.type==1){for(int h:g.labels[z.v]){I q=query(h,z.c);I val=z.w;
    if(q!=NEG)val=std::max(val,q+alpha-g.phi[size_t(z.v)*g.H+h]);C[off[i]+g.rank[size_t(z.v)*g.H+h]]=val;}}
   else if(e.type==2){I best=z.w;for(int k=off[i];k<off[i+1];k++){C[k]=std::max(C[k],R[i]);best=std::max(best,C[k]);}
    B[i]=best;for(int k=off[i];k<off[i+1];k++)C[k]=std::max(C[k],best-gamma);}
   else for(int h:g.labels[z.v])update(h,z.d,C[off[i]+g.rank[size_t(z.v)*g.H+h]]+z.y+z.d+g.phi[size_t(z.v)*g.H+h]);
  }
 }
 void check(const I* values,const I* best,const I* remote,const I* back){
  for(size_t k=0;k<C.size();k++)require(back[k]>=-1&&back[k]<(I)C.size(),"pointer out of range");
  for(int i=0;i<n;i++){
   require(best[i]==B[i],"B mismatch");require(remote[i]==R[i],"R mismatch");const auto& z=a[i];
   for(int k=off[i];k<off[i+1];k++){
    require(values[k]==C[k],"C mismatch");I q=back[k];int h=g.labels[z.v][k-off[i]];
    require(q>=-1&&q<(I)C.size(),"pointer out of range");
    if(q==-1){require(values[k]==z.w,"restart score mismatch");continue;}
    int j=owner[q];const auto& u=a[j];int source_h=g.labels[u.v][q-off[j]];
    if(j==i){require(source_h!=h,"self switch");require(values[q]==B[i]&&values[k]==values[q]-gamma,"same switch score mismatch");
     require(back[q]==-1||owner[back[q]]!=i,"consecutive same-anchor switches");}
    else {
     require(u.d<z.c,"non-strict read precedence");I gap;
     if(source_h==h){require(u.v<z.v||(u.v==z.v&&u.y<z.x),"non-strict local graph precedence");
      gap=g.phi[size_t(z.v)*g.H+h]+z.x-g.phi[size_t(u.v)*g.H+h]-u.y+z.c-u.d-2;
      require(values[k]==values[q]+z.w-gap,"local pointer score mismatch");}
     else {require(g.rank[size_t(z.v)*g.H+source_h]<0,"remote source present at target");
      int w=g.last[size_t(z.v)*g.H+source_h];require(w>=0&&u.v<=w,"illegal remote exit");
      gap=g.phi[size_t(w)*g.H+source_h]+g.dist[size_t(z.v)*g.H+source_h]+z.x-g.phi[size_t(u.v)*g.H+source_h]-u.y+z.c-u.d-2;
      require(values[k]==values[q]+z.w-gap-gamma,"remote pointer score mismatch");}
     require(gap>=0,"negative gap");
    }
   }
  }
  std::vector<unsigned char> color(C.size());
  for(I start=0;start<(I)C.size();start++)if(!color[start]){
   I k=start;while(k>=0&&!color[k]){color[k]=1;k=back[k];}
   require(k<0||color[k]!=1,"traceback cycle");k=start;
   while(k>=0&&color[k]==1){color[k]=2;k=back[k];}
  }
 }
};
extern "C" {
const char* validation_error(){return error.c_str();}
void* validation_index(int v,int h,int e,const I*l,const int*edges,const int*po,const int*ps){try{return new Index(v,h,e,l,edges,po,ps);}catch(const std::exception&x){error=x.what();return nullptr;}}
void validation_index_free(void*p){delete static_cast<Index*>(p);}
void* validation_batch(void*p,int n,const I*a,I gamma){try{return new Batch(*static_cast<Index*>(p),n,a,gamma);}catch(const std::exception&x){error=x.what();return nullptr;}}
void validation_batch_free(void*p){delete static_cast<Batch*>(p);}
int validation_check(void*p,const I*C,const I*B,const I*R,const I*back){try{static_cast<Batch*>(p)->check(C,B,R,back);return 0;}catch(const std::exception&x){error=x.what();return 1;}}
}
