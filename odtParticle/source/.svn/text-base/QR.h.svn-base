#ifndef QR_H
#define QR_H

#include<stdio.h>
#include<iostream>
#include<math.h>
#include<vector>

// Solves linear system, A*x=b using vector objects in the fuction QR
// QR(vector vector A, vector b, vector x)
// This algorithm crashes if A b and x do not have consistant dimensions, uses the vector .size() function.
//


/*multiply two matrixex*/
void multmat(std::vector<std::vector<double> > &a,std::vector <std::vector<double> > &b,std::vector<std::vector<double> > &c) {
    int N=a.size();
    int i,j,k;
    for(i=0;i<N;i++)
        for(j=0;j<N;j++)
        {
            c[i][j]=0;
            for(k=0;k<N;k++)
                c[i][j]+=a[i][k]*b[k][j];
        }
}


/*multiply two matrixex*/
void multmat2(std::vector<double>  &a,std::vector<std::vector<double> > &b,std::vector<double>  &c) { 
    int j,k;
    int N=a.size();

    for(j=0;j<N;j++)
    {
        c[j]=0;
        for(k=0;k<N;k++)
            c[j]+=a[k]*b[k][j];
    }
}

/*seek the sign of a variable*/
    int sign(double &x) { 
        if(x<=0)
            return(-1);
        else 
            return(1);
    }

/*compute the "||x-y||*||x-y||"*/
double twomo(std::vector<double> &x) {
    int i;
    int N=x.size();
    double sum=0;
    for(i=0;i<N;i++)
        sum+=pow(x[i],2);
    /*sum=sqrt(sum);*/
    return (sum);
}

/*acquire (x-y)(x-y)T*/
void vectormat(std::vector<double> &a,std::vector<std::vector<double> > &c) {

    int N=c.size();
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++)
            c[i][j]=a[i]*a[j];

}


void ETA_Aq::QR(std::vector<std::vector<double> > &a,std::vector<double> &b,std::vector<double> &x){

    int N=b.size();
    int i,j,k;

    std::vector< std::vector<double> >H(N, std::vector<double> (N)) ;
    std::vector< std::vector<double> >temp(N,std::vector<double> (N));
    std::vector< std::vector<double> >Q(N,std::vector<double> (N));

    std::vector<double> h(N);
    double sum;
    for(k=0;k<N-1;k++)
    {
        sum=0;
        for(i=k;i<N;i++)
            sum+=a[i][k]*a[i][k];
        sum=sqrt(sum);
        for(i=k;i<N;i++)
        {	  
            if(i==k)
                h[i]=sum*sign(a[k][k])+a[k][k];
            else
                h[i]=a[i][k];
        }
        if(k>0)
            for(i=0;i<k;i++)
                h[i]=0;
        vectormat(h,H);
        for(i=0;i<N;i++)
            for(j=0;j<N;j++)
                if(i==j)
                    H[i][j]=1-2*H[i][j]/twomo(h);
                else
                    H[i][j]=-2*H[i][j]/twomo(h);
        multmat(H,a,temp);
        for(i=0;i<N;i++)
            for(j=0;j<N;j++)
                a[i][j]=temp[i][j];
        if(k==0)
            for(i=0;i<N;i++)
                for(j=0;j<N;j++)
                    Q[i][j]=H[i][j];
        else
        {
            multmat(Q,H,temp);
            for(i=0;i<N;i++)
                for(j=0;j<N;j++)
                    Q[i][j]=temp[i][j];
        }

    }

    std::vector<double>  y(N,0) ; // Intermediate array
    multmat2(b,Q,y);


    // Solves Diagnolized matrix
    for (i=N-1; i>=0 ;i--){
        for   (j=N-1; j>i ; j--){
            y[i]-=a[i][j]*x[j];
        }
        x[i]=y[i]/a[i][i];
    }





    return;
}


////////////////////////////////////////////////////////////////////////////////////////////

void QRdecomp(std::vector<std::vector<double> > &a, std::vector<std::vector<double> > &Q){
    // makes Q and R out of a.  A becomes R.
    int N=a.size();
    int i,j,k;

    std::vector< std::vector<double> >H(N, std::vector<double> (N)) ;
    std::vector< std::vector<double> >temp(N,std::vector<double> (N));

    std::vector<double> h(N);
    double sum;
    for(k=0;k<N-1;k++)
    {
        sum=0;
        for(i=k;i<N;i++)
            sum+=a[i][k]*a[i][k];
        sum=sqrt(sum);
        for(i=k;i<N;i++)
        {	  
            if(i==k)
                h[i]=sum*sign(a[k][k])+a[k][k];
            else
                h[i]=a[i][k];
        }
        if(k>0)
            for(i=0;i<k;i++)
                h[i]=0;
        vectormat(h,H);
        for(i=0;i<N;i++)
            for(j=0;j<N;j++)
                if(i==j)
                    H[i][j]=1-2*H[i][j]/twomo(h);
                else
                    H[i][j]=-2*H[i][j]/twomo(h);
        multmat(H,a,temp);
        for(i=0;i<N;i++)
            for(j=0;j<N;j++)
                a[i][j]=temp[i][j];
        if(k==0)
            for(i=0;i<N;i++)
                for(j=0;j<N;j++)
                    Q[i][j]=H[i][j];
        else
        {
            multmat(Q,H,temp);
            for(i=0;i<N;i++)
                for(j=0;j<N;j++)
                    Q[i][j]=temp[i][j];
        }

    }



    return;
}

////////////////////////////////////////////////////////////////////////////////////////////

void solve_x_QR(std::vector<std::vector<double> > &Q, 
        std::vector<std::vector<double> > &R,
        std::vector<double> &b,
        std::vector<double> &x) {

    int N = b.size();

    std::vector<double>  y(N,0) ; // Intermediate array
    multmat2(b,Q,y);


    // Solves Diagnolized matrix
    for (int i=N-1; i>=0 ;i--){
        for   (int j=N-1; j>i ; j--){
            y[i]-=R[i][j]*x[j];
        }
        x[i]=y[i]/R[i][i];
    }
}

#endif


