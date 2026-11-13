#include <Rcpp.h>
#include <vecLib/cblas.h>
#include <vecLib/clapack.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <limits>
#include <cstring>
#include <cstdint>
using namespace Rcpp;
using Vec = std::vector<float>;

// Owning R raw buffers retain FP32 bits across R calls and serialization.
// Numeric inputs remain accepted for boundary and equivalence tests.
// Column-major float matrix view; input is an R matrix or a packed raw buffer.
struct Matrix32 {
  RObject owner;
  Vec owned;
  const float* data;
  int rows, cols;
  explicit Matrix32(SEXP value) : owner(value) {
    if (Rf_inherits(value, "smm_fp32_matrix")) {
      List x(value);
      rows=as<int>(x["nrow"]); cols=as<int>(x["ncol"]);
      RawVector raw=x["data"];
      if(rows<0 || cols<0 || raw.size()!=static_cast<R_xlen_t>(rows)*cols*sizeof(float))
        stop("Invalid FP32 buffer dimensions");
      if(reinterpret_cast<uintptr_t>(RAW(raw)) % alignof(float)) stop("Unaligned FP32 buffer");
      data=reinterpret_cast<const float*>(RAW(raw));
    } else {
      NumericMatrix x(value); rows=x.nrow(); cols=x.ncol();
      owned.assign(x.begin(),x.end()); data=owned.data();
    }
  }
  int nrow() const { return rows; }
  int ncol() const { return cols; }
  size_t size() const { return static_cast<size_t>(rows)*cols; }
  const float* begin() const { return data; }
  const float* end() const { return data+size(); }
  float operator[](size_t i) const { return data[i]; }
  float operator()(int row,int col) const { return data[row+rows*col]; }
};

// Pack column-major FP32 values x with dimensions rows and cols for the R boundary.
List packed_matrix(const Vec& x,int rows,int cols) {
  RawVector raw(x.size()*sizeof(float));
  if(!x.empty()) std::memcpy(RAW(raw),x.data(),raw.size());
  List result=List::create(_["data"]=raw,_["nrow"]=rows,_["ncol"]=cols);
  result.attr("class")="smm_fp32_matrix";
  return result;
}

namespace {
static_assert(sizeof(float) == 4, "FP32 requires 4-byte float");
// Compute column-major C = alpha op(A) op(B) + beta C. ar,ac,br,bc are matrix dimensions; ta,tb request transposition.
void product(const float* A, int ar, int ac, bool ta,
             const float* B, int br, int bc, bool tb, Vec& out) {
  const int m = ta ? ac : ar, n = tb ? br : bc, k = ta ? ar : ac;
  if (k != (tb ? bc : br) || out.size() != static_cast<size_t>(m*n))
    stop("Invalid native matrix product dimensions");
  cblas_sgemm(CblasColMajor, ta ? CblasTrans : CblasNoTrans,
    tb ? CblasTrans : CblasNoTrans, m, n, k, 1.0f, A, ar, B, br, 0.0f, out.data(), m);
}

// Sum squared entries of x in the half-open slice [begin, begin + length).
float squares(const Vec& x, size_t begin, size_t length) {
  float value = 0;
  for (size_t i = begin; i < begin + length; ++i) value += (float)x[i]*x[i];
  return static_cast<float>(value);
}
// Sum entries of x in the half-open slice [begin, begin + length).
float total(const Vec& x, size_t begin, size_t length) {
  float value = 0;
  for (size_t i = begin; i < begin + length; ++i) value += x[i];
  return static_cast<float>(value);
}
// Convert column-major FP32 values x, dimensions rows/cols, into an R double matrix.
NumericMatrix matrix_copy(const Vec& x, int rows, int cols) {
  NumericMatrix result(rows, cols);
  std::copy(x.begin(), x.end(), result.begin());
  return result;
}

// Apply a scalar, temporal or separable spatial quadratic inverse.
struct Inverse {
  bool isotropic, spatial;
  float scale;
  Vec diagonal, lower1, lower2, row_basis, column_basis, denominator;
  int p, q, row_blocks=1;
  Inverse() : isotropic(true), spatial(false), scale(1), p(0), q(0) {}
  Inverse(const List& input, int pp, int qq) : p(pp), q(qq) {
    row_blocks = input.containsElementNamed("row_blocks") ? as<int>(input["row_blocks"]) : 1;
    isotropic = as<bool>(input["isotropic"]);
    spatial = input.containsElementNamed("spatial") && as<bool>(input["spatial"]);
    if (isotropic) scale = as<float>(input["diagonal"]);
    else if (spatial) {
      row_basis = as<Vec>(input["row_basis"]);
      column_basis = as<Vec>(input["column_basis"]);
      denominator = as<Vec>(input["denominator"]);
      if ((int)row_basis.size()!=p*p || (int)column_basis.size()!=q*q ||
          (int)denominator.size()!=p*q) stop("Invalid spatial spectral inverse");
    } else {
      diagonal = as<Vec>(input["diagonal"]);
      lower1 = as<Vec>(input["lower1"]); lower2 = as<Vec>(input["lower2"]);
      if ((int)diagonal.size() != q || (int)lower1.size() != std::max(0,q-1) ||
          (int)lower2.size() != std::max(0,q-2)) stop("Invalid band factor");
    }
  }
  void apply(const Vec& input, Vec& output, int columns) const {
    if (isotropic) {
      for (size_t i = 0; i < input.size(); ++i) output[i] = input[i] / scale;
      return;
    }
    const int d = p*q;
    if (spatial) {
      Vec left(d), spectral(d), back(d), answer(d);
      for (int c=0;c<columns;++c) {
        product(row_basis.data(),p,p,true,input.data()+c*d,p,q,false,left);
        product(left.data(),p,q,false,column_basis.data(),q,q,false,spectral);
        for (int i=0;i<d;++i) spectral[i]/=denominator[i];
        product(row_basis.data(),p,p,false,spectral.data(),p,q,false,back);
        product(back.data(),p,q,false,column_basis.data(),q,q,true,answer);
        std::copy(answer.begin(),answer.end(),output.begin()+c*d);
      }
      return;
    }
    for (int c = 0; c < columns; ++c) for (int r = 0; r < p; ++r) {
      for (int t = 0; t < q; ++t) {
        const int i = c*d + r + p*t;
        float value = input[i];
        if (t > 0) value -= lower1[t-1]*output[i-p];
        if (t > 1) value -= lower2[t-2]*output[i-2*p];
        output[i] = value / diagonal[t];
      }
      for (int t = q-1; t >= 0; --t) {
        const int i = c*d + r + p*t;
        float value = output[i];
        if (t+1 < q) value -= lower1[t]*output[i+p];
        if (t+2 < q) value -= lower2[t]*output[i+2*p];
        output[i] = value / diagonal[t];
      }
    }
  }
};

// Primal/observation-space Cholesky solve with optional iterative refinement.
struct LinearSolve {
  int n, d, k, type, dimension, p, q;
  float diagonal;
  const float* X;
  Vec factor, V, sample, feature, intermediate;
  Inverse inverse;
  LinearSolve(int nn, int dd, int kk, int pp, int qq, const float* xx) :
    n(nn), d(dd), k(kk), p(pp), q(qq), X(xx), sample(n*k), feature(d*k), intermediate(d*k) {}
  void reset(const List& value) {
    type = as<int>(value["type"]);
    Matrix32 R(value["factor"]);
    dimension = R.nrow();
    if (R.ncol() != dimension || dimension != (type == 0 ? d : n))
      stop("Invalid Cholesky factor dimensions");
    factor.assign(R.begin(), R.end());
    if (type == 1) diagonal = as<float>(value["diagonal"]);
    if (type == 2) {
      Matrix32 packed_V(value["V"]); V.assign(packed_V.begin(),packed_V.end());
      inverse = Inverse(as<List>(value["inverse"]), p, q);
    }
  }
  void triangular(Vec& value) const {
    cblas_strsm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit,
      dimension, k, 1.0f, factor.data(), dimension, value.data(), dimension);
    cblas_strsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
      dimension, k, 1.0f, factor.data(), dimension, value.data(), dimension);
  }
  void solve(const Vec& rhs, Vec& answer) {
    if (type == 0) { answer = rhs; triangular(answer); return; }
    if (type == 1) {
      product(X,n,d,false,rhs.data(),d,k,false,sample);
      triangular(sample);
      product(X,n,d,true,sample.data(),n,k,false,feature);
      for (int i = 0; i < d*k; ++i)
        answer[i] = rhs[i]/diagonal - feature[i]/(diagonal*diagonal);
    } else {
      inverse.apply(rhs, intermediate, k);
      product(X,n,d,false,intermediate.data(),d,k,false,sample);
      triangular(sample);
      product(V.data(),d,n,false,sample.data(),n,k,false,feature);
      for (int i = 0; i < d*k; ++i) answer[i] = intermediate[i] - feature[i];
    }
  }
};

// Reusable SVD workspace for nuclear proximal maps and spectral projections.
struct Spectral {
  int p, q, rank;
  Vec input, singular, u, vt, work, values_work;
  std::vector<__CLPK_integer> iwork;
  Spectral(int pp, int qq) : p(pp), q(qq), rank(std::min(p,q)),
    input(p*q), singular(rank), u(p*rank), vt(rank*q), iwork(8*rank) {
    for (int mode = 0; mode < 2; ++mode) {
      char job = mode == 0 ? 'S' : 'N';
      __CLPK_integer m=p,n=q,lda=p,ldu= mode == 0 ? p : 1,
        ldvt=mode == 0 ? rank : 1,info=0,lwork=-1;
      float query=0;
      sgesdd_(&job,&m,&n,input.data(),&lda,singular.data(),u.data(),&ldu,
        vt.data(),&ldvt,&query,&lwork,iwork.data(),&info);
      if (info != 0) stop("SVD workspace query failed");
      (mode == 0 ? work : values_work).resize(static_cast<size_t>(query));
    }
  }
  void decompose(const float* x, bool vectors) {
    for (int i = 0; i < p*q; ++i)
      if (!std::isfinite(x[i]))
        stop("SVD input contains non-finite value at index %d (p=%d, q=%d)", i, p, q);
    std::copy(x,x+p*q,input.begin());
    char job = vectors ? 'S' : 'N';
    Vec& workspace = vectors ? work : values_work;
    __CLPK_integer m=p,n=q,lda=p,ldu=vectors ? p : 1,ldvt=vectors ? rank : 1,
      info=0,lwork=workspace.size();
    sgesdd_(&job,&m,&n,input.data(),&lda,singular.data(),u.data(),&ldu,
      vt.data(),&ldvt,workspace.data(),&lwork,iwork.data(),&info);
    if (info < 0)
      stop("SGESDD invalid argument %d (p=%d, q=%d)", -info, p, q);
    if (info > 0) {
      // Retry only this matrix decomposition; SGESDD overwrites its input.
      const int gesdd_info = info;
      char jobu = vectors ? 'S' : 'N', jobvt = jobu;
      __CLPK_integer fallback_lwork = -1;
      float query = 0;
      std::copy(x, x+p*q, input.begin());
      sgesvd_(&jobu,&jobvt,&m,&n,input.data(),&lda,singular.data(),
        u.data(),&ldu,vt.data(),&ldvt,&query,&fallback_lwork,&info);
      if (info != 0 || !std::isfinite(query) || query < 1)
        stop("SGESVD workspace query failed (info=%d, SGESDD info=%d)",
          static_cast<int>(info), gesdd_info);
      Vec fallback_work(static_cast<size_t>(std::ceil(query)));
      fallback_lwork = fallback_work.size();
      std::copy(x, x+p*q, input.begin());
      sgesvd_(&jobu,&jobvt,&m,&n,input.data(),&lda,singular.data(),
        u.data(),&ldu,vt.data(),&ldvt,fallback_work.data(),&fallback_lwork,&info);
      if (info != 0)
        stop("SGESVD fallback failed (info=%d, SGESDD info=%d, p=%d, q=%d)",
          static_cast<int>(info), gesdd_info, p, q);
      Rcpp::warning("SGESDD info=%d; SGESVD fallback succeeded (p=%d, q=%d)",
        gesdd_info, p, q);
    }
  }
  float norm(const float* x, bool nuclear) {
    decompose(x,false);
    return nuclear ? total(singular,0,rank) : singular[0];
  }
  void shrink(const Vec& x, Vec& answer, int k, float threshold) {
    if (threshold == 0) { answer = x; return; }
    const int d=p*q;
    for (int c=0; c<k; ++c) {
      decompose(x.data()+c*d,true);
      int kept=0;
      for (int j=0; j<rank; ++j) {
        float value=std::max(0.0f,singular[j]-threshold);
        if (value > 0) ++kept;
        for (int i=0; i<p; ++i) u[i+p*j] *= value;
      }
      std::fill(answer.begin()+c*d,answer.begin()+(c+1)*d,0.0f);
      if (kept > 0) cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
        p,q,kept,1.0f,u.data(),p,vt.data(),rank,0.0f,answer.data()+c*d,p);
    }
  }
};

// Compute direct multiclass hinge proximal rows: a is input, answer output, labels 1-based, costs loss weights, rho ADMM penalty, h score scale, n/k row/class counts.
void hinge(const Vec& a, Vec& answer, const IntegerVector& labels,
           const Matrix32& costs, int n, int k, float rho, float h) {
  // Exact binary proximal map in margin units. The interior branch avoids
  // subtracting the large scaled multiplier from itself in FP32.
  if (k == 2) {
    for (int i=0; i<n; ++i) {
      const int y=labels[i]-1, w=1-y;
      const float ay=a[i+n*y], aw=a[i+n*w];
      const float cap=costs(i,0)/(rho*h), difference=aw-ay;
      if (difference <= -h) {
        answer[i+n*y]=ay; answer[i+n*w]=aw;
      } else if (difference >= 2*cap-h) {
        answer[i+n*y]=ay+cap; answer[i+n*w]=aw-cap;
      } else {
        const float midpoint=0.5f*ay+0.5f*aw;
        answer[i+n*y]=midpoint+0.5f*h;
        answer[i+n*w]=midpoint-0.5f*h;
      }
    }
    return;
  }
  Vec c(k-1), sorted(k-1);
  for (int i=0; i<n; ++i) {
    const int y=labels[i]-1;
    int j=0;
    for (int col=0; col<k; ++col) if (col!=y)
      c[j++]=1+a[i+n*col]/h-a[i+n*y]/h;
    sorted=c; std::sort(sorted.begin(),sorted.end(),std::greater<float>());
    float cumulative=0,total=0;
    for (j=0; j<k-1; ++j) {
      cumulative+=sorted[j];
      if (sorted[j]>cumulative/(j+2)) total=std::max(0.0f,cumulative/(j+2));
    }
    const float cost=costs(i,0)/(rho*h*h);
    float threshold=total;
    if (total>cost) {
      cumulative=0;
      for (j=0; j<k-1; ++j) {
        cumulative+=sorted[j];
        float candidate=(cumulative-cost)/(j+1);
        if (sorted[j]>candidate) threshold=candidate;
      }
    }
    float sum=0; j=0;
    for (int col=0; col<k; ++col) if(col!=y) {
      float alpha=std::max(0.0f,c[j++]-threshold);
      answer[i+n*col]=h*(a[i+n*col]/h-alpha); sum+=alpha;
    }
    answer[i+n*y]=h*(a[i+n*y]/h+sum);
  }
}

// Compute FP32 certificate candidates; final acceptance is checked in R FP64.
struct Certificate {
  Vec objective, dual, gap, B, N, linear, inverse, alpha;
  Certificate(int n,int d,int k) : objective(k),dual(k),gap(k),B(n*k),
    N(d*k),linear(d*k),inverse(d*k),alpha(n*k) {}
  void ovr(const Matrix32& X,const Matrix32& Y,const Matrix32& C,
           const Vec& W,const Vec& b,const Vec& U,const Vec& V,const Vec& raw,
           float rho,float h,float tau,float gamma,int p,int q,int k,
           Spectral& spectral,const Inverse& H0) {
    const int n=X.nrow(),d=X.ncol();
    for(int col=0;col<k;++col) {
      float plus=0,minus=0;
      for(int i=0;i<n;++i) {
        int index=i+n*col;
        alpha[index]=std::min(C[index],std::max(0.0f,-rho*h*V[index]));
        (Y[index]>0 ? plus : minus)+=alpha[index];
      }
      for(int i=0;i<n;++i) {
        int index=i+n*col;
        if(plus>minus && Y[index]>0) alpha[index] *= (float)minus/(float)plus;
        if(minus>plus && Y[index]<0) alpha[index] *= (float)plus/(float)minus;
        B[index]=Y[index]*alpha[index];
      }
      for(int i=0;i<d;++i) N[i+d*col]=tau==0 ? 0 : rho*U[i+d*col];
      if(tau>0) {
        float norm=spectral.norm(N.data()+d*col,false);
        if(norm>tau) for(int i=0;i<d;++i) N[i+d*col] *= tau/norm;
      }
    }
    product(X.begin(),n,d,true,B.data(),n,k,false,linear);
    for(int i=0;i<d*k;++i) linear[i]-=N[i];
    H0.apply(linear,inverse,k);
    for(int col=0;col<k;++col) {
      float rough=0,loss=0,dot=0;
      if (H0.spatial) {
        for (int t=0;t<q;++t) for(int r=0;r<p;++r) {
          int i=d*col+r+p*t;
          if((r+1)%(p/H0.row_blocks)!=0) { float delta=W[i]-W[i+1]; rough+=delta*delta; }
          if(t+1<q) { float delta=W[i]-W[i+p]; rough+=delta*delta; }
        }
      } else for(int t=0;t<q-2;++t) for(int r=0;r<p;++r) {
        float value=W[d*col+r+p*t]-2*W[d*col+r+p*(t+1)]+W[d*col+r+p*(t+2)];
        rough+=(float)value*value;
      }
      for(int i=0;i<n;++i) {
        int index=i+n*col;
        loss+=C[index]*std::max(0.0f,1-Y[index]*(raw[index]+b[col]));
      }
      for(int i=0;i<d;++i) dot+=(float)linear[i+d*col]*inverse[i+d*col];
      float nuclear=tau==0 ? 0 : spectral.norm(W.data()+d*col,true);
      objective[col]=squares(W,d*col,d)/2+tau*nuclear+gamma*(float)rough+(float)loss;
      dual[col]=total(alpha,n*col,n)-(float)dot/2;
      gap[col]=std::abs(objective[col]-dual[col])/
        std::max(1.0f,std::max(std::abs(objective[col]),std::abs(dual[col])));
    }
  }
  void direct(const Matrix32& X,const Matrix32& C,const IntegerVector& labels,
              const Vec& W,const Vec& U,const Vec& V,const Vec& raw,
              float rho,float h,float tau,float gamma,int p,int q,int k,Spectral& spectral,const Inverse& H0) {
    const int n=X.nrow(),d=X.ncol();
    float loss=0,sum_alpha=0;
    for(int i=0;i<n;++i) {
      const int y=labels[i]-1;
      float worst=0; float mass=0;
      for(int col=0;col<k;++col) if(col!=y) {
        worst=std::max(worst,raw[i+n*col]-raw[i+n*y]+1);
        B[i+n*col]=std::max(0.0f,rho*h*V[i+n*col]); mass+=B[i+n*col];
      }
      loss+=C(i,0)*worst;
      float scale=mass>C(i,0) ? C(i,0)/(float)mass : 1;
      float adjusted=0;
      for(int col=0;col<k;++col) if(col!=y) {
        B[i+n*col]*=scale; adjusted+=B[i+n*col];
      }
      B[i+n*y]=-(float)adjusted; sum_alpha+=adjusted;
    }
    float nuclear=0;
    for(int col=0;col<k;++col) {
      for(int i=0;i<d;++i) N[i+d*col]=tau==0 ? 0 : rho*U[i+d*col];
      if(tau>0) {
        nuclear+=spectral.norm(W.data()+d*col,true);
        float norm=spectral.norm(N.data()+d*col,false);
        if(norm>tau) for(int i=0;i<d;++i) N[i+d*col]*=tau/norm;
      }
    }
    product(X.begin(),n,d,true,B.data(),n,k,false,linear);
    for(int i=0;i<d*k;++i) linear[i]+=N[i];
    float rough=0,quadratic_dual=0;
    for(int col=0;col<k;++col) {
      if(H0.spatial) {
        for(int t=0;t<q;++t) for(int row=0;row<p;++row) {
          int i=col*d+row+p*t;
          if((row+1)%(p/H0.row_blocks)!=0) { float delta=W[i]-W[i+1]; rough+=delta*delta; }
          if(t+1<q) { float delta=W[i]-W[i+p]; rough+=delta*delta; }
        }
      } else for(int t=0;t<q-2;++t) for(int row=0;row<p;++row) {
        int i=col*d+row+p*t;
        float delta=W[i]-2*W[i+p]+W[i+2*p]; rough+=delta*delta;
      }
    }
    H0.apply(linear,inverse,k);
    for(int i=0;i<d*k;++i) quadratic_dual+=linear[i]*inverse[i];
    objective[0]=squares(W,0,d*k)/2+tau*nuclear+gamma*rough+loss;
    dual[0]=sum_alpha-quadratic_dual/2;
    gap[0]=std::abs(objective[0]-dual[0])/
      std::max(1.0f,std::max(std::abs(objective[0]),std::abs(dual[0])));
  }
};
}

// Run FP32 ADMM. X_input/Xc_input are raw/centered designs; Y_input and costs_input are binary labels/weights for OvR; labels are 1-based direct labels. initial holds W,S,U,M,V,b; control specifies penalties, tolerances, scaling, iteration budget and FP64 acceptance callback. factor(rho) supplies the linear factor; inverse_H0 supplies the quadratic inverse.
// [[Rcpp::export]]
List smm_admm_fp32(SEXP X_input, SEXP Xc_input, SEXP Y_input,
                    SEXP costs_input, IntegerVector labels, List initial, List control,
                    Function factor, List inverse_H0) {
  Matrix32 X(X_input), Xc(Xc_input), Y(Y_input), costs(costs_input);
  const bool ovr=as<std::string>(control["family"])=="ovr";
  const int n=X.nrow(),d=X.ncol(),p=as<int>(control["p"]),q=as<int>(control["q"]);
  NumericMatrix startW=initial["W"];
  const int k=startW.ncol(),limit=as<int>(control["max_iter"]),interval=as<int>(control["interval"]);
  const float tau=as<float>(control["tau"]),gamma=as<float>(control["gamma"]),
    h=as<float>(control["h"]),abs_tol=as<float>(control["abs_tol"]),
    rel_tol=as<float>(control["rel_tol"]),gap_tol=as<float>(control["gap_tol"]);
  float rho=as<float>(control["rho"]),diagnostic_rho=rho;
  const bool adaptive=as<bool>(control["adaptive_rho"]);
  const int refinement=control.containsElementNamed("refinement") ? as<int>(control["refinement"]) : 0;
  const float alpha=control.containsElementNamed("relaxation") ? as<float>(control["relaxation"]) : 1.0f;
  if(!std::isfinite(alpha) || alpha<=0 || alpha>=2) stop("ADMM relaxation must lie in (0,2)");
  Vec center=as<Vec>(control["center"]);
  if(n<1 || p<1 || q<1 || d!=p*q || k<2 || startW.nrow()!=d ||
     Xc.nrow()!=n || Xc.ncol()!=d || limit<1 || interval<0 ||
     (ovr && (Y.nrow()!=n || Y.ncol()!=k || costs.nrow()!=n || costs.ncol()!=k)) ||
     (!ovr && (labels.size()!=n || costs.nrow()!=n || costs.ncol()!=1)))
    stop("Invalid native ADMM inputs");
  Vec W=as<Vec>(initial["W"]),S=as<Vec>(initial["S"]),U=as<Vec>(initial["U"]),
    M=as<Vec>(initial["M"]),V=as<Vec>(initial["V"]),b=as<Vec>(initial["b"]);
  Vec Sold(d*k),Mold(n*k),target(n*k),rhs(d*k),raw(n*k),margin(n*k),
    shrink_input(d*k),a(n*k),dz(n*k),cross(d*k),yv(n*k),adjoint(d*k),
    means(k),primal(ovr?k:1),dual(ovr?k:1),ptol(ovr?k:1),dtol(ovr?k:1);
  if(S.size()!=W.size() || U.size()!=W.size() || M.size()!=(size_t)n*k ||
     V.size()!=M.size() || b.size()!=(size_t)k || center.size()!=(size_t)d)
    stop("Invalid native ADMM initial state");
  LinearSolve solver(n,d,k,p,q,Xc.begin());
  solver.reset(as<List>(factor(rho)));
  Spectral spectral(p,q);
  Inverse H0(inverse_H0,p,q);
  Certificate certificate(n,d,k);
  Function accept_fp64=control["accept_fp64"];
  int changes=0,iteration=0,checks=0,acceptance_checks=0,rejected_stops=0;
  bool converged=false;
  Vec linear_sample(n*k), linear_cross(d*k), linear_residual(d*k), correction(d*k);
  std::vector<int> hist_iteration;
  Vec hist_objective,hist_dual,hist_gap,hist_primal,hist_residual,hist_rho;
  for(iteration=1;iteration<=limit;++iteration) {
    if(iteration%25==0) checkUserInterrupt();
    Sold=S; Mold=M;
    if(ovr) {
      for(int col=0;col<k;++col) {
        for(int i=0;i<n;++i) {
          int index=i+n*col; target[index]=Y[index]*(M[index]-V[index])/h;
        }
        means[col]=total(target,n*col,n)/n;
        for(int i=0;i<n;++i) target[i+n*col]-=means[col];
      }
      product(Xc.begin(),n,d,true,target.data(),n,k,false,rhs);
      for(int i=0;i<d*k;++i) rhs[i]=rho*(S[i]-U[i]+h*h*rhs[i]);
    } else {
      for(int i=0;i<n*k;++i) target[i]=M[i]-V[i];
      product(X.begin(),n,d,true,target.data(),n,k,false,rhs);
      for(int i=0;i<d*k;++i) rhs[i]=rho*(S[i]-U[i]+h*rhs[i]);
    }
    solver.solve(rhs,W);
    // Residual correction in the original matrix-free equation, entirely FP32.
    // Difference products avoid the rounded, ill-conditioned band matrix.
    for(int refinement_step=0;refinement_step<refinement;++refinement_step) {
      product(Xc.begin(),n,d,false,W.data(),d,k,false,linear_sample);
      product(Xc.begin(),n,d,true,linear_sample.data(),n,k,false,linear_cross);
      for(int i=0;i<d*k;++i)
        linear_residual[i]=rhs[i]-(1+rho)*W[i]-rho*h*h*linear_cross[i];
      if (H0.spatial) {
        for(int col=0;col<k;++col) for(int t=0;t<q;++t) for(int row=0;row<p;++row) {
          int i=col*d+row+p*t;
          if((row+1)%(p/H0.row_blocks)!=0) { float v=2*gamma*(W[i]-W[i+1]);
            linear_residual[i]-=v; linear_residual[i+1]+=v; }
          if(t+1<q) { float v=2*gamma*(W[i]-W[i+p]);
            linear_residual[i]-=v; linear_residual[i+p]+=v; }
        }
      } else for(int col=0;col<k;++col) for(int row=0;row<p;++row)
        for(int t=0;t<q-2;++t) {
          int i=col*d+row+p*t;
          float value=2*gamma*((W[i]-W[i+p])-(W[i+p]-W[i+2*p]));
          linear_residual[i]-=value; linear_residual[i+p]+=2*value;
          linear_residual[i+2*p]-=value;
        }
      solver.solve(linear_residual,correction);
      for(int i=0;i<d*k;++i) W[i]+=correction[i];
    }
    if(ovr) for(int col=0;col<k;++col) {
      float shift=0;
      for(int i=0;i<d;++i) shift+=(float)center[i]*W[i+d*col];
      b[col]=means[col]-(float)shift;
    }
    product(X.begin(),n,d,false,W.data(),d,k,false,raw);
    for(int col=0;col<k;++col) for(int i=0;i<n;++i) {
      int index=i+n*col;
      margin[index]=ovr ? h*Y[index]*(raw[index]+b[col]) : h*raw[index];
    }
    for(int i=0;i<d*k;++i) shrink_input[i]=alpha*W[i]+(1-alpha)*Sold[i]+U[i];
    spectral.shrink(shrink_input,S,k,tau/rho);
    for(int i=0;i<n*k;++i) a[i]=alpha*margin[i]+(1-alpha)*Mold[i]+V[i];
    if(ovr) for(int i=0;i<n*k;++i) {
      const float step=costs[i]/(rho*h);
      // In the interior hinge branch the exact proximal value is h. Computing
      // a+(h-a) loses FP32 margin digits when the multiplier is large.
      M[i]=a[i]>=h ? a[i] : (a[i]<=h-step ? a[i]+step : h);
    }
    else hinge(a,M,labels,costs,n,k,rho,h);
    for(int i=0;i<d*k;++i) U[i]+=alpha*W[i]+(1-alpha)*Sold[i]-S[i];
    for(int i=0;i<n*k;++i) {
      V[i]+=alpha*margin[i]+(1-alpha)*Mold[i]-M[i];
      dz[i]=ovr ? h*Y[i]*(M[i]-Mold[i]) : h*(M[i]-Mold[i]);
      yv[i]=ovr ? Y[i]*V[i] : V[i];
    }
    product(X.begin(),n,d,true,dz.data(),n,k,false,cross);
    product(X.begin(),n,d,true,yv.data(),n,k,false,adjoint);
    for(int i=0;i<d*k;++i) {
      cross[i]+=S[i]-Sold[i];
      adjoint[i]=U[i]+h*adjoint[i];
    }
    const int groups=ovr?k:1;
    for(int col=0;col<groups;++col) {
      int dbegin=ovr?d*col:0,nbegin=ovr?n*col:0,dlen=ovr?d:d*k,nlen=ovr?n:n*k;
      float residual=0;
      for(int i=dbegin;i<dbegin+dlen;++i) residual+=(float)(W[i]-S[i])*(W[i]-S[i]);
      for(int i=nbegin;i<nbegin+nlen;++i) residual+=(float)(margin[i]-M[i])*(margin[i]-M[i]);
      primal[col]=std::sqrt((float)residual);
      float intercept_dual=ovr ? total(dz,nbegin,nlen) : 0;
      float intercept_adj=ovr ? total(yv,nbegin,nlen) : 0;
      dual[col]=rho*std::sqrt(squares(cross,dbegin,dlen)+intercept_dual*intercept_dual);
      ptol[col]=std::sqrt((float)(dlen+nlen))*abs_tol+rel_tol*std::max(
        std::sqrt(squares(W,dbegin,dlen)+squares(margin,nbegin,nlen)),
        std::sqrt(squares(S,dbegin,dlen)+squares(M,nbegin,nlen)));
      dtol[col]=std::sqrt((float)(dlen+(ovr?1:0)))*abs_tol+rel_tol*rho*
        std::sqrt(squares(adjoint,dbegin,dlen)+h*h*intercept_adj*intercept_adj);
    }
    bool residual_ok=true;
    for(int col=0;col<groups;++col)
      residual_ok=residual_ok && primal[col]<=ptol[col] && dual[col]<=dtol[col];
    diagnostic_rho=rho;
    int period=interval==0 ? (ovr?25:10) : interval;
    if(iteration==1 || iteration%period==0 || (interval==0 && residual_ok) || iteration==limit) {
      ++checks;
      if(ovr) certificate.ovr(X,Y,costs,W,b,U,V,raw,rho,h,tau,gamma,p,q,k,spectral,H0);
      else {
        certificate.direct(X,costs,labels,W,U,V,raw,rho,h,tau,gamma,p,q,k,spectral,H0);
        hist_iteration.push_back(iteration); hist_objective.push_back(certificate.objective[0]);
        hist_dual.push_back(certificate.dual[0]); hist_gap.push_back(certificate.gap[0]);
        hist_primal.push_back(primal[0]); hist_residual.push_back(dual[0]); hist_rho.push_back(rho);
      }
      converged=residual_ok;
      for(int col=0;col<groups;++col) converged=converged && certificate.gap[col]<=gap_tol;
      if(residual_ok) {
        ++acceptance_checks;
        // FP32 is only a screening certificate. The unchanged FP64 certificate
        // and residuals must accept the candidate before the iteration loop stops.
        converged=as<bool>(accept_fp64(List::create(
          _["W"]=matrix_copy(W,d,k),_["S"]=matrix_copy(S,d,k),
          _["U"]=matrix_copy(U,d,k),_["M"]=matrix_copy(M,n,k),
          _["V"]=matrix_copy(V,n,k),_["b"]=wrap(b),
          _["S_old"]=matrix_copy(Sold,d,k),_["M_old"]=matrix_copy(Mold,n,k),
          _["rho"]=rho,_["diagnostic_rho"]=diagnostic_rho)));
        if(converged) break;
        ++rejected_stops;
      }
    }
    if(adaptive && iteration%50==0 && iteration<=1000 && changes<10) {
      float next=rho;
      if (ovr) {
        float pp=squares(primal,0,groups),dd=squares(dual,0,groups);
        if(pp>100*dd) next=std::min(rho*2,1e6f);
        if(dd>100*pp) next=std::max(rho/2,1e-6f);
      } else {
        if(primal[0]>10*dual[0]) next=std::min(rho*2,1e6f);
        if(dual[0]>10*primal[0]) next=std::max(rho/2,1e-6f);
      }
      if(next!=rho) {
        for(float& x:U) x=x*rho/next;
        for(float& x:V) x=x*rho/next;
        rho=next; ++changes; solver.reset(as<List>(factor(rho)));
      }
    }
  }
  return List::create(_["W"]=matrix_copy(W,d,k),_["S"]=matrix_copy(S,d,k),
    _["U"]=matrix_copy(U,d,k),_["M"]=matrix_copy(M,n,k),_["V"]=matrix_copy(V,n,k),
    _["b"]=wrap(b),_["S_old"]=matrix_copy(Sold,d,k),_["M_old"]=matrix_copy(Mold,n,k),
    _["rho"]=rho,_["diagnostic_rho"]=diagnostic_rho,_["iterations"]=std::min(iteration,limit),
    _["converged"]=converged,_["certificate_checks"]=checks,
    _["fp64_acceptance_checks"]=acceptance_checks,_["fp64_rejected_stops"]=rejected_stops,
    _["history"]=DataFrame::create(_["iteration"]=hist_iteration,
      _["objective"]=hist_objective,_["dual_objective"]=hist_dual,_["relative_gap"]=hist_gap,
      _["primal_residual"]=hist_primal,_["dual_residual"]=hist_residual,_["rho"]=hist_rho));
}


namespace {
// Exact path-graph eigenvectors, evaluated and applied entirely in FP32.
// Construct the orthonormal path-Laplacian eigenbasis for n adjacent positions.
Vec path_basis32(int n) {
  Vec basis(n*n);
  const float pi=std::acos(-1.0f);
  for(int col=0;col<n;++col) for(int row=0;row<n;++row)
    basis[row+n*col]=std::sqrt((col==0 ? 1.0f : 2.0f)/n)*
      std::cos(pi*col*(row+0.5f)/n);
  return basis;
}
// Build the inverse representation from p,q,gamma,diagonal; spatial selects grid adjacency and row_blocks prevents cross-channel edges.
Inverse make_inverse32(int p, int q, float gamma, float diagonal, bool spatial=false, int row_blocks=1) {
  if(p<1 || q<1 || gamma<0 || diagonal<=0) stop("Invalid inverse parameters");
  if(row_blocks<1 || p%row_blocks!=0) stop("Invalid colour-plane dimensions");
  Inverse result; result.row_blocks=row_blocks;
  result.p = p; result.q = q; result.scale = diagonal; result.spatial=spatial;
  result.isotropic = gamma == 0;
  if (result.isotropic) return result;
  if (spatial) {
    const int height=p/row_blocks;
    Vec block=path_basis32(height); result.row_basis.assign(p*p,0);
    for(int c=0;c<row_blocks;++c) for(int j=0;j<height;++j) for(int i=0;i<height;++i)
      result.row_basis[c*height+i+p*(c*height+j)]=block[i+height*j];
    result.column_basis=path_basis32(q);
    result.denominator.resize(p*q);
    const float pi=std::acos(-1.0f);
    for(int col=0;col<q;++col) for(int row=0;row<p;++row) {
      float a=std::sin(pi*(row%height)/(2*height)), b=std::sin(pi*col/(2*q));
      result.denominator[row+p*col]=diagonal+8*gamma*(a*a+b*b);
    }
    return result;
  }
  Vec main(q, diagonal), first(std::max(0, q-1), 0.0f), second(std::max(0, q-2), 0.0f);
  for (int j=0; j<q-2; ++j) {
    main[j] += 2*gamma; main[j+1] += 8*gamma; main[j+2] += 2*gamma;
    first[j] -= 4*gamma; first[j+1] -= 4*gamma; second[j] += 2*gamma;
  }
  result.diagonal.resize(q); result.lower1.resize(first.size()); result.lower2.resize(second.size());
  for (int j=0; j<q; ++j) {
    float a=0, b=0;
    if (j>1) b=result.lower2[j-2]=second[j-2]/result.diagonal[j-2];
    if (j>0) a=result.lower1[j-1]=(first[j-1]-(j>1 ? b*result.lower1[j-2] : 0.0f))/result.diagonal[j-1];
    float pivot=main[j]-a*a-b*b;
    if (!(pivot>0)) stop("FP32 band Cholesky has nonpositive pivot");
    result.diagonal[j]=std::sqrt(pivot);
  }
  return result;
}
// Convert inverse representation x to the R list expected by the backend.
List inverse_list32(const Inverse& x) {
  if (x.isotropic) return List::create(_["isotropic"]=true, _["diagonal"]=x.scale, _["spatial"]=x.spatial, _["row_blocks"]=x.row_blocks);
  if (x.spatial) return List::create(_["isotropic"]=false, _["spatial"]=true, _["row_blocks"]=x.row_blocks,
    _["row_basis"]=wrap(x.row_basis), _["column_basis"]=wrap(x.column_basis),
    _["denominator"]=wrap(x.denominator));
  return List::create(_["isotropic"]=false, _["diagonal"]=wrap(x.diagonal),
    _["lower1"]=wrap(x.lower1), _["lower2"]=wrap(x.lower2));
}
// Factor positive-definite column-major matrix x in place; dimension is its order.
void cholesky32(Vec& x, int dimension) {
  char upper='U'; __CLPK_integer n=dimension, info=0;
  spotrf_(&upper, &n, x.data(), &n, &info);
  if (info!=0) stop("FP32 dense Cholesky failed");
  for (int col=0; col<n; ++col) for (int row=col+1; row<n; ++row) x[row+n*col]=0.0f;
}
}

// Construct the quadratic inverse for p-by-q coefficients, penalty gamma and diagonal. spatial selects grid first differences; row_blocks partitions independent colour planes.
// [[Rcpp::export]]
List smm_inverse_fp32(int p, int q, float gamma, float diagonal, bool spatial=false, int row_blocks=1) {
  return inverse_list32(make_inverse32(p,q,gamma,diagonal,spatial,row_blocks));
}

// Factor the ADMM system for input design and optional Gram matrix gram_input. p,q specify coefficient shape; gamma,rho,h specify structure, ADMM penalty and score scale; spatial,row_blocks specify the adjacency.
// [[Rcpp::export]]
List smm_factor_fp32(SEXP input, SEXP gram_input, int p, int q, float gamma, float rho, float h, bool spatial=false, int row_blocks=1) {
  Matrix32 X(input);
  const int n=X.nrow(),d=X.ncol();
  if (n<1 || d!=p*q || rho<=0 || h<=0 || gamma<0) stop("Invalid FP32 factor inputs");
  const float diagonal=1+rho;
  if (d<=n) {
    Matrix32 gram(gram_input); Vec factor(gram.begin(),gram.end());
    if (factor.size()!=static_cast<size_t>(d*d)) stop("Invalid FP32 feature Gram size");
    for (float& value:factor) value*=rho*h*h;
    for (int i=0;i<d;++i) factor[i+d*i]+=diagonal;
    if (spatial) {
      auto edge=[&](int i,int j) {
        factor[i+d*i]+=2*gamma; factor[j+d*j]+=2*gamma;
        factor[i+d*j]-=2*gamma; factor[j+d*i]-=2*gamma;
      };
      for(int t=0;t<q;++t) for(int r=0;r<p;++r) {
        int i=r+p*t;
        if((r+1)%(p/row_blocks)!=0) edge(i,i+1);
        if(t+1<q) edge(i,i+p);
      }
    } else {
    const float c[3]={1.0f,-2.0f,1.0f};
    for(int t=0;t<q-2;++t) for(int a=0;a<3;++a) for(int b=0;b<3;++b)
      for(int r=0;r<p;++r) factor[r+p*(t+a)+d*(r+p*(t+b))]+=2*gamma*c[a]*c[b];
    }
    cholesky32(factor,d);
    return List::create(_["type"]=0, _["factor"]=packed_matrix(factor,d,d));
  }
  Vec factor(n*n);
  if(gamma==0) {
    Matrix32 gram(gram_input); factor.assign(gram.begin(),gram.end());
    if (factor.size()!=static_cast<size_t>(n*n)) stop("Invalid FP32 sample Gram size");
    for(float& value:factor) value/=diagonal;
    for(int i=0;i<n;++i) factor[i+n*i]+=1/(rho*h*h);
    cholesky32(factor,n);
    return List::create(_["type"]=1, _["diagonal"]=diagonal, _["factor"]=packed_matrix(factor,n,n));
  }
  Inverse inverse=make_inverse32(p,q,gamma,diagonal,spatial,row_blocks);
  Vec transpose(d*n),V(d*n);
  for(int j=0;j<n;++j) for(int i=0;i<d;++i) transpose[i+d*j]=X(j,i);
  inverse.apply(transpose,V,n);
  product(X.begin(),n,d,false,V.data(),d,n,false,factor);
  for(int i=0;i<n;++i) factor[i+n*i]+=1/(rho*h*h);
  cholesky32(factor,n);
  return List::create(_["type"]=2, _["factor"]=packed_matrix(factor,n,n),
    _["V"]=packed_matrix(V,d,n), _["inverse"]=inverse_list32(inverse));
}

// Return native precision/backend metadata; takes no arguments.
// [[Rcpp::export]]
List smm_fp32_info() {
  return List::create(_["scalar_bytes"]=sizeof(float), _["significand_bits"]=std::numeric_limits<float>::digits,
    _["matrix_product"]="cblas_sgemm", _["triangular_solve"]="cblas_strsm",
    _["svd"]="sgesdd", _["cholesky"]="spotrf", _["reductions"]="float", _["spatial_operator"]="FP32_four_neighbour_graph_DCT");
}
// Return the smaller primal or observation Gram matrix of input, a numeric or packed FP32 design.
// [[Rcpp::export]]
List smm_gram_fp32(SEXP input) {
  Matrix32 X(input);
  const int n=X.nrow(),d=X.ncol(),m=std::min(n,d);
  Vec result(m*m);
  product(X.begin(),n,d,d<=n,X.begin(),n,d,d>n,result);
  return packed_matrix(result,m,m);
}

// Apply the penalty inverse to input right-hand sides. p,q,gamma,diagonal,spatial,row_blocks define the temporal or disconnected spatial operator.
// [[Rcpp::export]]
NumericMatrix smm_apply_inverse_fp32(NumericMatrix input, int p, int q, float gamma, float diagonal, bool spatial=false, int row_blocks=1) {
  if(input.nrow()!=p*q) stop("Invalid inverse right-hand side");
  Vec x(input.begin(),input.end()), answer(x.size());
  make_inverse32(p,q,gamma,diagonal,spatial,row_blocks).apply(x,answer,input.ncol());
  return matrix_copy(answer,p*q,input.ncol());
}
// Solve for rhs using the precomputed factor and input design; p,q define the coefficient dimensions. Used by numerical validation.
// [[Rcpp::export]]
NumericMatrix smm_apply_factor_fp32(SEXP input, NumericMatrix rhs, List factor, int p, int q) {
  Matrix32 X(input);
  if(rhs.nrow()!=p*q) stop("Invalid solve right-hand side");
  LinearSolve solver(X.nrow(),p*q,rhs.ncol(),p,q,X.begin());
  solver.reset(factor);
  Vec b(rhs.begin(),rhs.end()), answer(b.size()); solver.solve(b,answer);
  return matrix_copy(answer,p*q,rhs.ncol());
}

// Pack numeric matrix X into a column-major FP32 representation for repeated calls.
// [[Rcpp::export]]
List smm_pack_fp32(NumericMatrix X) {
  Vec values(X.begin(), X.end());
  return packed_matrix(values, X.nrow(), X.ncol());
}
