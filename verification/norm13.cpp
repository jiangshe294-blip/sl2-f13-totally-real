// Independent of the norm13 producer: integer pseudo-remainders, q13 Newton
// traces over E, and characteristic-coefficient Newton recurrence over M.
#define main old_exact_verifier_main
#include "exact_V.hpp"
#undef main
struct E{Vec n=Vec(14);Z d=1;};struct M{E a,b;};
void en(E&x){Z c=x.d;for(Z&z:x.n)mpz_gcd(c.get_mpz_t(),c.get_mpz_t(),z.get_mpz_t());if(x.d<0)c=-c;if(c!=1){x.d/=c;for(Z&z:x.n)z/=c;}}
E es(Q q){E x;x.n[0]=q.get_num();x.d=q.get_den();return x;}
E ea(const E&x,const E&y,int sign=1){E z;mpz_lcm(z.d.get_mpz_t(),x.d.get_mpz_t(),y.d.get_mpz_t());Z sx=z.d/x.d,sy=sign*(z.d/y.d);for(int i=0;i<14;i++)z.n[i]=sx*x.n[i]+sy*y.n[i];en(z);return z;}
E esc(E x,Q q){for(Z&z:x.n)z*=q.get_num();x.d*=q.get_den();en(x);return x;}
Z G13;long ep=0;
E em(const E&x,const E&y){Vec p(27);for(int i=0;i<14;i++)if(x.n[i]!=0)for(int j=0;j<14;j++)if(y.n[j]!=0)mpz_addmul(p[i+j].get_mpz_t(),x.n[i].get_mpz_t(),y.n[j].get_mpz_t());for(int k=26;k>=14;k--){Z top=p[k];for(int i=0;i<k;i++)p[i]*=G;for(int j=0;j<14;j++)p[k-14+j]-=top*gi[j];p[k]=0;}E z;z.d=x.d*y.d*G13;for(int i=0;i<14;i++)z.n[i]=p[i];en(z);ep++;cap();return z;}
E readE(){std::vector<Q>x(14);E r;for(Q&q:x){q=rq();mpz_lcm(r.d.get_mpz_t(),r.d.get_mpz_t(),q.get_den().get_mpz_t());}for(int i=0;i<14;i++)r.n[i]=x[i].get_num()*(r.d/x[i].get_den());en(r);return r;}
E getcol(const V&x,int j){E e;e.d=x.d;for(int i=0;i<14;i++)e.n[i]=x.a[i][j];en(e);return e;}
E onlyE(const V&x){for(int i=0;i<14;i++)for(int j=1;j<13;j++)if(x.a[i][j]!=0)throw std::runtime_error("expected E coordinate");return getcol(x,0);}
bool equalE(const E&a,const E&b){for(int i=0;i<14;i++)if(a.n[i]*b.d!=b.n[i]*a.d)return false;return true;}
void eqV(const V&a,const V&b){for(int i=0;i<14;i++)for(int j=0;j<13;j++)if(a.a[i][j]*b.d!=b.a[i][j]*a.d)throw std::runtime_error("V mismatch");}
V eV(const E&e,const V&v){std::array<E,13> cols;V w;for(int j=0;j<13;j++){cols[j]=em(e,getcol(v,j));mpz_lcm(w.d.get_mpz_t(),w.d.get_mpz_t(),cols[j].d.get_mpz_t());}for(int j=0;j<13;j++)for(int i=0;i<14;i++)w.a[i][j]=cols[j].n[i]*(w.d/cols[j].d);norm(w);return w;}
M ma(const M&x,const M&y,int sign=1){return{ea(x.a,y.a,sign),ea(x.b,y.b,sign)};}
M msc(const M&x,Q q){return{esc(x.a,q),esc(x.b,q)};}
M mm(const M&x,const M&y,const E&ac){return{ea(em(x.a,y.a),em(ac,em(x.b,y.b))),ea(em(x.a,y.b),em(x.b,y.a))};}
M readM(){V a=readv(),b=readv();return{onlyE(a),onlyE(b)};}
void outV(const V&v){std::cout<<"{\"denominator\":\""<<v.d<<"\",\"numerator_rows_x_then_y\":[";for(int i=0;i<14;i++){if(i)std::cout<<',';std::cout<<'[';for(int j=0;j<13;j++){if(j)std::cout<<',';std::cout<<'"'<<v.a[i][j]<<'"';}std::cout<<']';}std::cout<<"]}";}
V embed(const E&e){V v;v.d=e.d;for(int i=0;i<14;i++)v.a[i][0]=e.n[i];return v;}
void outM(const M&m){std::cout<<'[';outV(embed(m.a));std::cout<<',';outV(embed(m.b));std::cout<<']';}
int main(){alarm(35);try{
 int k;std::cin>>k;if(k<1||k>13)throw std::runtime_error("stage range");
 std::vector<Q>g(15);G=1;for(Q&q:g){q=rq();mpz_lcm(G.get_mpz_t(),G.get_mpz_t(),q.get_den().get_mpz_t());}if(g[14]!=1)throw std::runtime_error("nonmonic g");for(int i=0;i<15;i++)gi[i]=g[i].get_num()*(G/g[i].get_den());mpz_pow_ui(G38.get_mpz_t(),G.get_mpz_t(),38);mpz_pow_ui(G13.get_mpz_t(),G.get_mpz_t(),13);E ac=readE();
 V u0=readv(),u1=readv(),v0=readv(),v1=readv();M expected=readM();
 std::vector<M>p(k+1),e(k+1);e[0].a=es(1);for(int i=1;i<k;i++){p[i]=readM();e[i]=readM();}
 V power0,power1;if(k==1){power0=u0;power1=u1;}else{V bb=mul(v1,u1);power0=add(mul(v0,u0),eV(ac,bb));power1=add(mul(v0,u1),mul(v1,u0));}
 // Derive q13(a,Y)=g(Y)/(Y-a) and its first twelve Newton traces.
 std::array<E,14>q;for(int j=0;j<=13;j++){std::vector<Q>vals(14);for(int l=0;l<=13-j;l++)vals[l]=g[j+l+1];E c;for(Q&z:vals)mpz_lcm(c.d.get_mpz_t(),c.d.get_mpz_t(),z.get_den().get_mpz_t());for(int l=0;l<14;l++)c.n[l]=vals[l].get_num()*(c.d/vals[l].get_den());en(c);q[j]=c;}
 std::array<E,13>tr;tr[0]=es(13);for(int j=1;j<=12;j++){tr[j]=esc(q[13-j],-j);for(int i=1;i<j;i++)tr[j]=ea(tr[j],em(q[13-i],tr[j-i]),-1);}
 auto trace=[&](const V&v){E total;for(int j=0;j<13;j++)total=ea(total,em(getcol(v,j),tr[j]));return total;};p[k]={trace(power0),trace(power1)};
 M sum;for(int i=1;i<=k;i++)sum=ma(sum,msc(mm(e[k-i],p[i],ac),Q(i%2?1:-1)));e[k]=msc(sum,Q(1,k));
 if(k==13&&(!equalE(e[k].a,expected.a)||!equalE(e[k].b,expected.b)))throw std::runtime_error("recomputed norm13 differs from published beta");
 std::cout<<"{\"status\":\""<<(k==13?"PASS_RECOMPUTED_NORM13_EQUALS_PUBLISHED_BETA":"PASS_RECOMPUTED_NORM13_STAGE")<<"\",\"stage\":"<<k<<",\"power_u_k\":[";outV(power0);std::cout<<',';outV(power1);std::cout<<"],\"trace_u_k\":";outM(p[k]);std::cout<<",\"elementary_e_k\":";outM(e[k]);std::cout<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}return 0;}
