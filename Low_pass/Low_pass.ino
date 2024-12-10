//LOW-PASS FILTER
//It will be used in previous project to filter signal
//It is more of education project and I know that there are dozens of libraries for filtering noisy signals 
#define PI 3.1415926535897932384626433832795
double xn1 = 0;
double yn1 = 0;
int k = 0;

// long nCr(int n, int r)
// {
//     return fact(n) / (fact(r) * fact(n - r));
// }

// // Returns factorial of n
// long fact(int n)
// {
//       if(n==0)
//       return 1;
//     long res = 1;
//     for (int i = 2; i <= n; i++)
//         res = res * i;
//     return res;
// }


template <int order>
class LowPassFilter
{
  private:
    double a[order+1];
    double b[order+1];
    double omega0;
    double dt;
    bool adapt;
    double tn1 = 0;
    double x[order+1]; // Raw values
    double y[order+1]; // Filtered values
  public:  
    LowPassFilter(double f0, double fs, bool adaptive){
      // f0: cutoff frequency (Hz)
      // fs: sample frequency (Hz)
      // adaptive: boolean flag, if set to 1, the code will automatically set the sample frequency based on the time history.
      
      omega0 = 2*PI*f0;
      dt = 1.0/fs;
      adapt = adaptive;
      tn1 = -dt;
      for(int k = 0; k < order+1; k++){
        x[k] = 0;
        y[k] = 0;        
      }
      setCoef();
    }

    void PolynomialMultiply(double A[], int nA, double B[], int nB){ // array1,length of array1, array2(always size of 3, except if order is odd)
        double result[nB+nA-1]; 

        for(int i = 0; i<nB+nA-1; i++)
          result[i] = 0.0;

        for(int i = 0; i<nA; i++)
          for(int j = 0; j<nB; j++)
            result[i+j]+=A[i]*B[j];
        
        for(int i = 0; i<nB+nA-1; i++)
          A[i] = result[i];
      return;
    }
    
    void setCoef(){
      if(adapt){
        double t = micros()/1.0e6;
        dt = t - tn1;
        tn1 = t;
      }
      
      double alpha = omega0*dt;
      double alphaSq = alpha*alpha;
      
      a[0] = 1;
      int length_a = 1;
      if(order%2 == 1){
        a[0]=alpha+2;
        a[1]=alpha-2;
        length_a = 2;
      }
      for(int k = 1; k<= (order - (order%2))/2; k++){
        double coef = 4.0*alpha*cos((2.0*k+order-1.0)/(2.0*order)*PI);
        double element_wo_coef = 4+alphaSq;
        double Bn[3] = {element_wo_coef - coef,-8.0 + 2.0*alphaSq,element_wo_coef+coef};
        PolynomialMultiply(a,length_a,Bn,3);
        length_a+=2;
      }
      double a0 = a[0];

      b[0] = alpha/a0;
      b[1] = alpha/a0;
      double temp[2] = {alpha,alpha};
      for(int i = 2; i<=order; i++){
        PolynomialMultiply(b,i,temp,2);
      }

      for(int i = 0; i<=order; i++){
        a[i]/=-a0;
      }

    }

     double filt(double xn){
      // Provide me with the current raw value: x
      // I will give you the current filtered value: y
      y[0] = 0.0;
      x[0] = xn;
      // Compute the filtered values
      for(int k = 0; k < order; k++){
        y[0] += a[k+1]*y[k+1] + b[k]*x[k];
      }
      y[0] += b[order]*x[order];
      for(int k = order; k > 0; k--){
        y[k] = y[k-1];
        x[k] = x[k-1];
      }
  
      // Return the filtered value    
      return y[0];
    }
};







LowPassFilter<6> lp(3,1000,false);

void setup() {
  Serial.begin(115200);
}

void loop() { 
 
  // Test signal
  double t = micros()/1.0e6;
  double xn = sin(2*PI*2*t) + 1.0*sin(2*PI*10*t);

  // Compute the filtered signal
  double yn = lp.filt(xn);


  if(k % 3 == 0){
    // This extra conditional statement is here to reduce
    // the number of times the data is sent through the Serial port
    // because sending data through the Serial port
    // messes with the sampling frequency
  
    // Output
    Serial.print(xn);
    Serial.print(" ");
    Serial.println(yn);
  }
  k = k+1;
}