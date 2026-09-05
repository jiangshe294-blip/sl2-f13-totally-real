#include <gmpxx.h>
#include <array>
#include <vector>
#include <iostream>
#include <chrono>
#include <stdexcept>
#include <unistd.h>
// New univariate E arithmetic. No producer includes or reduction tables.
using Z=mpz_class;using Q=mpq_class;
using Vec=std::vector<Z>;
struct E{Vec c=Vec(14);Z d=1;};struct M{E x,y;};
Z G,G13;Vec gi(15);E rad;
const auto start=std::chrono::steady_clock::now();int multiplications=0;
void ck(bool b,const char*s){if(!b)throw std::runtime_error(s);}
Z rz(){std::string s;ck(bool(std::cin>>s),"integer input");return Z(s);}
Q rq(){std::string s;ck(bool(std::cin>>s),"rational input");Q q(s);ck(q.get_den()!=0,"rational denominator");q.canonicalize();return q;}
bool zero(const E&a){for(const Z&x:a.c)if(x!=0)return false;return true;}
void norm(E&a){ck(a.d>0,"denominator must be positive");Z h=a.d;for(const Z&x:a.c){mpz_gcd(h.get_mpz_t(),h.get_mpz_t(),x.get_mpz_t());if(h==1)break;}if(h!=1){a.d/=h;for(Z&x:a.c)x/=h;}}
E readE(){E a;a.d=rz();for(Z&x:a.c)x=rz();norm(a);return a;}
E scal(const Q&q){E a;a.c[0]=q.get_num();a.d=q.get_den();return a;}
E add(const E&a,const E&b,int sign=1){E z;mpz_lcm(z.d.get_mpz_t(),a.d.get_mpz_t(),b.d.get_mpz_t());Z u=z.d/a.d,v=sign*(z.d/b.d);for(int i=0;i<14;i++)z.c[i]=u*a.c[i]+v*b.c[i];norm(z);return z;}
E mul(const E&a,const E&b){multiplications++;if(zero(a)||zero(b))return E();Vec p(27);for(int i=0;i<14;i++)if(a.c[i]!=0)for(int j=0;j<14;j++)if(b.c[j]!=0)mpz_addmul(p[i+j].get_mpz_t(),a.c[i].get_mpz_t(),b.c[j].get_mpz_t());
 // Each step replaces p by G*p - top*a^(k-14)*G*g.
 // Its rational meaning is preserved by multiplying the denominator by G.
 for(int k=26;k>=14;k--){Z top=p[k];for(Z&x:p)x*=G;for(int i=0;i<15;i++)p[k-14+i]-=top*gi[i];}
 E z;z.d=a.d*b.d*G13;for(int i=0;i<14;i++)z.c[i]=p[i];for(int i=14;i<27;i++)ck(p[i]==0,"unreduced high power");norm(z);return z;}
M prod(const M&a,const M&b){
 E ac=mul(a.x,b.x),bd=mul(a.y,b.y);
 E cross=add(add(mul(add(a.x,a.y),add(b.x,b.y)),ac,-1),bd,-1);
 return {add(ac,mul(rad,bd)),cross};
}
void outE(const E&a){std::cout<<"{\"denominator\":\""<<a.d<<"\",\"numerator_coefficients_a_ascending\":[";for(int i=0;i<14;i++){if(i)std::cout<<',';std::cout<<'"'<<a.c[i]<<'"';}std::cout<<"]}";}
void init(){std::vector<Q>g(15),r(14);G=1;for(Q&q:g){q=rq();mpz_lcm(G.get_mpz_t(),G.get_mpz_t(),q.get_den().get_mpz_t());}ck(g[14]==1,"g not monic");for(int i=0;i<15;i++)gi[i]=g[i].get_num()*(G/g[i].get_den());mpz_pow_ui(G13.get_mpz_t(),G.get_mpz_t(),13);
 rad.d=1;for(Q&q:r){q=rq();mpz_lcm(rad.d.get_mpz_t(),rad.d.get_mpz_t(),q.get_den().get_mpz_t());}for(int i=0;i<14;i++)rad.c[i]=r[i].get_num()*(rad.d/r[i].get_den());norm(rad);
}
int main(){alarm(35);try{std::string mode;std::cin>>mode;
 if(mode=="horner"){
  int k;std::cin>>k;ck(0<=k&&k<28,"coefficient index");init();M beta{readE(),readE()},previous{readE(),readE()};Q coefficient=rq();M now=prod(previous,beta);now.x=add(now.x,scal(coefficient));
  if(k==0)ck(zero(now.x)&&zero(now.y),"P28(beta) IS NONZERO");
  std::cout<<"{\"status\":\""<<(k==0?"PASS_EXACT_P28_BETA_ZERO":"PASS_EXACT_HORNER_STEP")<<"\",\"coefficient_index\":"<<k<<",\"residual_pair\":[";outE(now.x);std::cout<<',';outE(now.y);std::cout<<"],\"E_multiplications\":"<<multiplications<<",\"zero_coordinates\":"<<(k==0?28:0)<<",\"native_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
 }else if(mode=="encoding"){
  Z D=rz();ck(D>0,"positive scale");int begin,end;std::cin>>begin>>end;ck(0<=begin&&begin<end&&end<=57,"encoding range");Z content=0;int checked=0;
  for(int i=begin;i<end;i++){Q c=rq();Z z=rz();ck(c.get_num()*D==z*c.get_den(),"primitive coefficient differs from scaled rational coefficient");if(i%2)ck(c==0&&z==0,"odd coefficient nonzero");if(i==56)ck(c==1&&z==D,"monic leading coefficient mismatch");mpz_gcd(content.get_mpz_t(),content.get_mpz_t(),z.get_mpz_t());checked++;}
  std::cout<<"{\"status\":\"PASS_EXACT_POLYNOMIAL_ENCODING_SLICE\",\"start\":"<<begin<<",\"end\":"<<end<<",\"coefficients_checked\":"<<checked<<",\"slice_integer_content\":\""<<content<<"\",\"native_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
 }else if(mode=="content"){
  int n;std::cin>>n;Z h=0;for(int i=0;i<n;i++){Z z=rz();mpz_gcd(h.get_mpz_t(),h.get_mpz_t(),z.get_mpz_t());}ck(h==1,"integer polynomial not primitive");std::cout<<"{\"status\":\"PASS_GLOBAL_PRIMITIVE_CONTENT\",\"content\":1}\n";
 }else throw std::runtime_error("unknown mode");
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}return 0;}
