#include "likelihood.h"

//Construtor
template<class V, class M>
Likelihood<V,M>::Likelihood(const char* prefix, const VectorSet<V, M>& domain, const double* data_mean, const double* t, const double* stdDevs, const double& poi, const double& s):BaseScalarFunction<V,M>(prefix,domain), m_data_mean(0), m_t(0), m_stdDevs(0){
    size_t const size = sizeof(data_mean)/sizeof(double);
    m_data_mean.assign(data_mean,data_mean+size); 
    size_t const sizet = s;
    m_t.assign(t,t+sizet);
    m_stdDevs.assign(stdDevs,stdDevs+size);
    //cout<<Likelihood::values_of_a.size()<<endl;
    this->poi = poi;
}

template<class V, class M>
double Likelihood<V, M>::lnValue(const V& domainVector) const{
    double a = domainVector[0];
    cout<<domainVector[0]<<" ";
    //Likelihood::values_of_a.push_back(a);
    double model = exp(-a*m_t[poi]);
    double error_A = abs((model - m_data_mean[0])/m_stdDevs[0]);
    cerr<<error_A<<" ";
    double misfitValue = error_A*error_A;
    return -1.5*misfitValue;
}


template<class V, class M>
Likelihood<V,M>::~Likelihood(){}

template<class V, class M>
double Likelihood<V,M>::actualValue(const V & domainVector, const V * domainDirection, V * gradVector, M * hessianMatrix, V * hessianEffect) const{
    return exp(this->lnValue(domainVector));
}

 template class Likelihood<GslVector,GslMatrix>;
