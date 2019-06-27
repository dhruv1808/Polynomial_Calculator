#include <iostream>
#include <cmath>
#include <algorithm>
#include "Polynomial.h"


void init_poly(poly_t &p, double const init_coeffs[],unsigned int const init_degree)
{
	if(p.a_coeffs!=nullptr)
	{
		delete[] p.a_coeffs;
		p.a_coeffs=nullptr;
	}
	p.a_coeffs=new double[init_degree+1];
	p.degree=init_degree;
	for(unsigned int i=0;i<=init_degree;i++)
	{
		p.a_coeffs[i]=init_coeffs[i];
		std::cout<<p.a_coeffs[i]<<" ";
	}
}

void destroy_poly(poly_t &p)
{
	delete[] p.a_coeffs;
	p.a_coeffs=nullptr;
}

unsigned int poly_degree(poly_t const &p)
{
	if(p.a_coeffs==nullptr)
	{
		throw 0;
	}
	return p.degree;
}

double poly_coeff(poly_t const &p, unsigned int n)
{
	if(p.a_coeffs==nullptr)
	{
		throw 0;
	}
	if(n>p.degree)
	{
		return 0;
	}
	else
	{
	return p.a_coeffs[n];
	}
}

double poly_val(poly_t const &p,double const x)
{
	if(p.a_coeffs==nullptr)
	{
		throw 0;
	}
	double sum;
	for(unsigned int i=0;i<=p.degree;i++)
	{
		sum=sum+(p.a_coeffs[i]*pow(x,i));
	}
	return sum;
}

void poly_add(poly_t &p, poly_t const &q)
{
    if(p.a_coeffs==nullptr || q.a_coeffs == nullptr)
    {
        throw 0;
    }
    int max = std::max(p.degree,q.degree);
    double *temp = new double [max+1];
    if(p.degree == q.degree)
    {
    	for(unsigned int i=0;i<=p.degree;i++)
    	{
    		temp[i] =p.a_coeffs[i] + q.a_coeffs[i];
    	}
    }
    else if(p.degree>q.degree)
    {
        for(unsigned int i=0;i<=p.degree;i++)
        {
            temp[i] = p.a_coeffs[i];
            std::cout<<temp[i]<<" ";
        }
        for(unsigned int i=0;i<=q.degree;i++)
        {
            temp[i] = temp[i]+q.a_coeffs[i];
            std::cout<<temp[i]<<" ";
        }
    }
    else
    {
        for(unsigned int i=0;i<=q.degree;i++)
        {
            temp[i] = q.a_coeffs[i];
            std::cout<<temp[i]<<" ";
        }
        for(unsigned int i=0;i<=p.degree;i++)
        {
            temp[i] =temp[i]+p.a_coeffs[i];
            std::cout<<temp[i]<<" ";
        }
    }
    std::cout<<max;
    delete[] p.a_coeffs;
    p.a_coeffs = temp;
    for(int i=0;i<=max;i++)
    {
        p.a_coeffs[i] = temp[i];
    }
    p.degree = max;
    if(p.a_coeffs[max]==0)
    {
        p.degree--;
    }
}

void poly_subtract(poly_t &p, poly_t const &q){
    if(p.a_coeffs==nullptr || q.a_coeffs == nullptr){
        throw 0;
    }
    int max = std::max(p.degree,q.degree);
    double *temp = new double [max+1];

    if(p.degree == q.degree){
        for(unsigned int i=0;i<=p.degree;i++){
            temp[i] =p.a_coeffs[i] - q.a_coeffs[i];
        }
    }
    else if(p.degree>q.degree){
        for(unsigned int i=0;i<=p.degree;i++){
            temp[i] = p.a_coeffs[i];
            std::cout<<temp[i]<<" ";
        }
        for(unsigned int i=0;i<=q.degree;i++){
            temp[i] -=q.a_coeffs[i];
            std::cout<<temp[i]<<" ";
        }
    }
    else{
        for(unsigned int i=0;i<=q.degree;i++){
            temp[i] = (-1)*q.a_coeffs[i];
            std::cout<<temp[i]<<" ";
        }
        for(unsigned int i=0;i<=p.degree;i++){
            temp[i] +=p.a_coeffs[i];
            std::cout<<temp[i]<<" ";
        }
    }
    std::cout<<max;
    delete[] p.a_coeffs;
    p.a_coeffs = temp;
    for(int i=0;i<=max;i++){
        p.a_coeffs[i] = temp[i];
    }
    p.degree = max;
    if(p.a_coeffs[max]==0){
        p.degree--;
    }

}


void poly_multiply(poly_t &p, poly_t const &q)
{
	 if(p.a_coeffs==nullptr||q.a_coeffs==nullptr)
	 {
	    throw 0;
	 }
	 double *temp = new double[p.degree + q.degree + 1]{};
	 if((q.a_coeffs[0] == 0 && q.degree==0) || (p.a_coeffs[0] == 0 && p.degree==0) )
	 {
	    delete[] p.a_coeffs;
	    p.a_coeffs = new double[1];
	    p.a_coeffs[0] = 0;
	    p.degree = 0;
	 }
	 else
	 {
	    for(unsigned int i= 0; i<=p.degree;i++)
	    {
	    	for(unsigned int j=0;j<=q.degree;j++)
	    	{
	    		temp[i+j] += p.a_coeffs[i]*q.a_coeffs[j];
	    	}
	    }
	 delete[] p.a_coeffs;
	 p.a_coeffs = nullptr;
	 p.a_coeffs = new double[p.degree+q.degree+1];
	 for(unsigned int j=0;j<=(p.degree+q.degree);j++)
	 	 {
	     	 p.a_coeffs[j] = temp[j];
	 	 }
	 p.degree = p.degree+q.degree;
	 }
	 std::cout<<p.degree;
	 delete[] temp;
	 temp = nullptr;
}


double poly_divide(poly_t &p, double r)
{
    if(p.a_coeffs==nullptr)
    	{
        	throw 0;
    	}
    double a = 0.0;
    double *temp = new double [p.degree+1];
    if(p.degree==0)
    	{
        	double temp1 = p.a_coeffs[0];
        	delete[] p.a_coeffs;
        	p.a_coeffs = new double[1]{0};
        	std::cout<<r;
        	a = temp1;
    	}
    else
    {
        std::cout<<"1 element"<<p.a_coeffs[p.degree]<<" ";
        temp[0] = p.a_coeffs[p.degree];
        for(int i =p.degree-1,j=1;i>=0;i--,j++)
        {
        	temp[j] = p.a_coeffs[i] + r*temp[j-1];
        }
        a = temp[p.degree];
        std::cout<<temp[0];
        delete[] p.a_coeffs;
        p.a_coeffs = nullptr;
        p.degree--;
        p.a_coeffs = new double[p.degree+1];
        for(std::size_t i = 0; i<=p.degree;i++)
        	{
            	p.a_coeffs[i] = temp[p.degree-i];
        	}
        std::cout<<"degree is "<<p.degree<<" res poly is ";
    }
    delete[] temp;
    temp = nullptr;
    return a;
}


void poly_diff(poly_t &p)
{
	if(p.a_coeffs==nullptr)
	{
		throw 0;
	}
	if(p.degree==0)
	{
		delete[] p.a_coeffs;
		p.a_coeffs=new double[1];
		p.a_coeffs[0]=0;
	}
	else
	{
		unsigned int a=p.degree;
		double *temp=new double[p.degree];
		for(unsigned int i=1,n=1;i<=a;i++,n++)
		{
			temp[i-1]=n*p.a_coeffs[i];
		}
		p.degree--;
		delete[] p.a_coeffs;
		p.a_coeffs=temp;
		for(unsigned int i=0;i<=p.degree;i++)
		{
			p.a_coeffs[i]=temp[i];
		}
	}
}

double poly_approx_int(poly_t const &p, double a, double b, unsigned int n)
{
	if(p.a_coeffs==nullptr)
	{
		throw 0;
	}
	double h=((b-a)/n);
	double sum=0;
	for(unsigned int i=0;i<=n;i++)
	{
		double x = a+(i*h);
		if( i==0 || i==n )
		{
			sum=sum+poly_val(p,x);
		}
		else
		{
			sum=sum+2*poly_val(p,x);
		}
	}
	return h*(sum/2);
}

int main()
{
	/*
	poly_t p{nullptr,0},q{nullptr,0};
	double init_coeffs[5] = {1,2,3,4,5};
	double init_coeffs1[5] = {1,2,-3,4,5};
	init_poly(p, init_coeffs, 4);
	std::cout<<std::endl;
	init_poly(q, init_coeffs1, 2);
	std::cout<<std::endl;
	std::cout<<poly_degree(p);
	std::cout<<std::endl;
	std::cout<<poly_coeff(p,2);
	std::cout<<std::endl;
	std::cout<<poly_val(p,1);
	std::cout<<std::endl;
	poly_add(p,q);
	std::cout<<std::endl;
	poly_subtract(p,q);
	std::cout<<std::endl;
	poly_multiply(p,q);
	std::cout<<std::endl;
	poly_diff(p);
	std::cout<<std::endl;
	std::cout<<poly_approx_int(p,1,2,3);
	*/
	return 0;
}
