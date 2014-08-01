/******************************************************************************
 * Version: 1.0
 * Last modified on: 30th July, 2014 
 * Developers: Michael G. Epitropakis, Xiaodong Li.
 *      email: mge_(AT)_cs_(DOT)_stir_(DOT)_ac_(DOT)_uk 
 *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
 * c# port by Keith Nelson : keith_(AT)_cmind_(DOT)_org
 * ***************************************************************************/
using System;
using System.Collections.Generic;

namespace cec2013
{
    public class Funcs {
	
	    public class CF1 : CFunc {

		    public CF1(IBoundsFn boundsFn) 
			    :base(boundsFn, 6) {
			
		        for (int i=0; i<nofunc_; ++i) {
		            sigma_[i] = 1;
		            bias_[i]  = 0;
		            weight_[i]= 0;
		        }
		        lambda_[0] = 1.0;
		        lambda_[1] = 1.0;
		        lambda_[2] = 8.0;
		        lambda_[3] = 8.0;
		        lambda_[4] = 1.0/5.0;
		        lambda_[5] = 1.0/5.0;
		        /* load optima */
		        if (this.getDimension() == 2 || this.getDimension() == 3 || this.getDimension() == 5
		                || this.getDimension() == 10 || this.getDimension() == 20 ) {
		            String fname = "../../data/CF1_M_D" + this.getDimension() + "_opt.dat";
		            //String root = System.getProperty("user.dir");
    //		        System.out.println(root+ " : " +fname);
		            loadOptima(fname);
		        } else {
		    	    Console.WriteLine("Error: NOT SUPPOSED TO BE HERE");
		            initOptimaRandomly();
		        }
		        /* M_ Identity matrices */
		        initRotmatIdentity();
		        /* Initialise functions of the composition */
                funcs_ = new List<Func>()
                {
		            new FGriewank(boundsFn),
		            new FGriewank(boundsFn),
		            new FWeierstrass(boundsFn),
		            new FWeierstrass(boundsFn),
		            new FSphere(boundsFn),
		            new FSphere(boundsFn)
                };
		        //TODO: calculate this correctly in the initialisation phase
		        CalculateFMaxi();
		    }

		    public override double doEvaluate(double[] x) {
			    return this.evaluateInner_(x);
		    }
		
	    }//CF1
	
	    public class CF2 : CFunc {

		    public CF2(IBoundsFn boundsFn)
			    :base( boundsFn, 8){
		        for (int i=0; i<nofunc_; ++i) {
		            sigma_[i] = 1.0;
		            bias_[i]  = 0.0;
		            weight_[i]= 0.0;
		        }
		        lambda_[0] = 1.0;
		        lambda_[1] = 1.0;
		        lambda_[2] = 10.0;
		        lambda_[3] = 10.0;
		        lambda_[4] = 1.0/10.0;
		        lambda_[5] = 1.0/10.0;
		        lambda_[6] = 1.0/7.0;
		        lambda_[7] = 1.0/7.0;
		        /* load optima */
		        if (this.getDimension() == 2 || this.getDimension() == 3 || this.getDimension() == 5
		                || this.getDimension() == 10 || this.getDimension() == 20 ) {
		            String fname = "../../data/CF2_M_D" + this.getDimension() + "_opt.dat";
		            loadOptima(fname);
		        } else {
		            initOptimaRandomly();
		        }
		    
		        /* M_ Identity matrices */
		        initRotmatIdentity();
		    
		        /* Initialise functions of the composition */
                funcs_ = new List<Func>()
                {
		            new FRastrigin(boundsFn),
		            new FRastrigin(boundsFn),
		            new FWeierstrass(boundsFn),
		            new FWeierstrass(boundsFn),
		            new FGriewank(boundsFn),
		            new FGriewank(boundsFn),
		            new FSphere(boundsFn),
		            new FSphere(boundsFn)
                };

		        CalculateFMaxi();
		    }

		    
		    public override double doEvaluate(double[] x) {
			    return this.evaluateInner_(x);
		    }
	    }//CF2
	
	    public class CF3 : CFunc {

		    public CF3(IBoundsFn boundsFn)
			    :base( boundsFn, 6){

			    for (int i=0; i<nofunc_; ++i) {
		            bias_[i]  = 0.0;
		            weight_[i]= 0.0;
		        }
		        sigma_[0] = 1.0;
		        sigma_[1] = 1.0;
		        sigma_[2] = 2.0;
		        sigma_[3] = 2.0;
		        sigma_[4] = 2.0;
		        sigma_[5] = 2.0;
		        lambda_[0] = 1.0/4.0;
		        lambda_[1] = 1.0/10.0;
		        lambda_[2] = 2.0;
		        lambda_[3] = 1.0;
		        lambda_[4] = 2.0;
		        lambda_[5] = 5.0;
		        /* load optima */
		        if (this.getDimension() == 2 || this.getDimension() == 3 || this.getDimension() == 5
		                || this.getDimension() == 10 || this.getDimension() == 20 ) {
		            String fname;
		            fname = "../../data/CF3_M_D" + this.getDimension() + "_opt.dat";
		            loadOptima(fname);
		            fname = "../../data/CF3_M_D" + this.getDimension() + ".dat";
		            loadRotationMatrix(fname);
		        } else {
		            initOptimaRandomly();
		            //TODO: Generate dimension independent rotation matrices
		            /* M_ Identity matrices */
		            initRotmatIdentity();
		        }
		        /* Initialise functions of the composition */
                funcs_ = new List<Func>()
                {
		            new FEF8F2(boundsFn),
		            new FEF8F2(boundsFn),
		            new FWeierstrass(boundsFn),
		            new FWeierstrass(boundsFn),
		            new FGriewank(boundsFn),
		            new FGriewank(boundsFn),
                };
		        CalculateFMaxi();
		    }

		    
		    public override double doEvaluate(double[] x) {
			    return this.evaluateInner_(x);
		    }
		
	    }//CF3
	
	    public class CF4 : CFunc {

		    public CF4(IBoundsFn boundsFn)
			    :base( boundsFn, 8){
			
		        for (int i=0; i<nofunc_; ++i) {
		            sigma_[i] = 1.0;
		            bias_[i]  = 0.0;
		            weight_[i]= 0.0;
		        }
		        sigma_[0] = 1.0;
		        sigma_[1] = 1.0;
		        sigma_[2] = 1.0;
		        sigma_[3] = 1.0;
		        sigma_[4] = 1.0;
		        sigma_[5] = 2.0;
		        sigma_[6] = 2.0;
		        sigma_[7] = 2.0;
		        lambda_[0] = 4.0;
		        lambda_[1] = 1.0;
		        lambda_[2] = 4.0;
		        lambda_[3] = 1.0;
		        lambda_[4] = 1.0/10.0;
		        lambda_[5] = 1.0/5.0;
		        lambda_[6] = 1.0/10.0;
		        lambda_[7] = 1.0/40.0;
		        /* load optima */
		        if (this.getDimension() == 2 || this.getDimension() == 3 || this.getDimension() == 5
		                || this.getDimension() == 10 || this.getDimension() == 20) {
		            String fname;
		            fname = "../../data/CF4_M_D" + this.getDimension() + "_opt.dat";
		            loadOptima(fname);
		            fname = "../../data/CF4_M_D" + this.getDimension() + ".dat";
		            loadRotationMatrix(fname);
		        } else {
		            initOptimaRandomly();
		            //TODO: make dimension independent rotation matrices
		            /* M_ Identity matrices */
		            initRotmatIdentity();
		        }
		        /* Initialise functions of the composition */
                funcs_ = new List<Func>()
                {
		            new FRastrigin(boundsFn),
		            new FRastrigin(boundsFn),
		            new FEF8F2(boundsFn),
		            new FEF8F2(boundsFn),
		            new FWeierstrass(boundsFn),
		            new FWeierstrass(boundsFn),
		            new FGriewank(boundsFn),
		            new FGriewank(boundsFn)
                };
		        CalculateFMaxi();
		    }

		    
		    public override double doEvaluate(double[] x) {
			    return this.evaluateInner_(x);
		    }
		
	    }//CF4
	
	    /**
	     *  F1: Five-Uneven-Peak Trap
	     *  Variable ranges: x in [0, 30]
	     *  No. of global peaks: 2
	     *  No. of local peaks:  3.
	     */
	    public class FiveUnevenPeakTrap : Func {

		    public FiveUnevenPeakTrap(IBoundsFn boundsFn)
			    :base( boundsFn){}

		    
		    public override double doEvaluate(double[] x) 
		    {
		        double result=-1.0;
		        if (x[0]>=0 && x[0]< 2.5) {
		            result = 80*(2.5-x[0]);
		        } else if (x[0] >= 2.5 && x[0] < 5.0) {
		            result = 64*(x[0]-2.5);
		        } else if (x[0] >= 5.0 && x[0] < 7.5) {
		            result = 64*(7.5-x[0]);
		        } else if (x[0] >= 7.5 && x[0] < 12.5) {
		            result = 28*(x[0]-7.5);
		        } else if (x[0] >= 12.5 && x[0] < 17.5) {
		            result = 28*(17.5-x[0]);
		        } else if (x[0] >= 17.5 && x[0] < 22.5) {
		            result = 32*(x[0]-17.5);
		        } else if (x[0] >= 22.5 && x[0] < 27.5) {
		            result = 32*(27.5-x[0]);
		        } else if (x[0] >= 27.5 && x[0] <= 30) {
		            result = 80*(x[0]-27.5);
		        }

		        return result;
		    }
		
	    }// FiveUnevenPeakTrap

	    /**
	     * F2: Equal Maxima
	     * Variable ranges: x in [0, 1]
	     * No. of global peaks: 5
	     * No. of local peaks:  0. 
	     */
	    public class EqualMaxima : Func {

		    public EqualMaxima(IBoundsFn boundsFn) 
			    :base( boundsFn){}
		    
		    public override double doEvaluate(double[] x) {
		        double s = Math.Sin(5.0 * Math.PI * x[0]);
		        return Math.Pow(s, 6);
		    }
	    }//EqualMaxima

	    /**
	     * F3: Uneven Decreasing Maxima
	     * Variable ranges: x in [0, 1]
	     * No. of global peaks: 1
	     * No. of local peaks:  4. 
	     */
	    public class UnevenDecreasingMaxima : Func {

		    public UnevenDecreasingMaxima(IBoundsFn boundsFn) 
			    :base( boundsFn){}
		    
		    public override double doEvaluate(double[] x) {
		        double tmp1 = -2*Math.Log(2.0)*((x[0]-0.08)/0.854)*((x[0]-0.08)/0.854);
		        double tmp2 = Math.Sin( 5*Math.PI*(Math.Pow(x[0],3.0/4.0)-0.05) );
		        return Math.Exp(tmp1) * Math.Pow(tmp2, 6);
		    }
	    }// UnevenDecreasingMaxima
	
	    /**
	     * F4: Himmelblau
	     * Variable ranges: x, y in [−6, 6
	     * No. of global peaks: 4
	     * No. of local peaks:  0.
	     */
	    public class Himmelblau : Func {

		    public Himmelblau(IBoundsFn boundsFn)
			    :base( boundsFn){}

		    
		    public override double doEvaluate(double[] x) {
		        return 200.0 - (x[0]*x[0] + x[1] - 11.0)*(x[0]*x[0] + x[1] - 11.0) -
			            (x[0] + x[1]*x[1] - 7.0)*(x[0] + x[1]*x[1] - 7.0);
		    }
	    }//Himmelblau
	
	    /**
	     * F5: Six-Hump Camel Back
	     * Variable ranges: x in [−1.9, 1.9]; y in [−1.1, 1.1]
	     * No. of global peaks: 2
	     * No. of local peaks:  2.
	     */
	    public class SixHumpCamelBack : Func {

		    public SixHumpCamelBack(IBoundsFn boundsFn)
                : base(boundsFn) {}
		    
		    public override double doEvaluate(double[] x) {
			    return -( (4.0 - 2.1*x[0]*x[0] + Math.Pow(x[0],4.0)/3.0)*x[0]*x[0] +
			            x[0]*x[1] + (4.0*x[1]*x[1] -4.0)*x[1]*x[1] );
		    }
		
	    }//SixHumpCamelBack
	
	    /**
	     * F6: Shubert
	     * Variable ranges: x_i in  [−10, 10]^n, i=1,2,...,n
	     * No. of global peaks: n*3^n
	     * No. of local peaks: many
	     */
	    public class Shubert : Func {

		    public Shubert(IBoundsFn boundsFn)
                : base(boundsFn) {
		    }

		    
		    public override double doEvaluate(double[] x) {
		        double result = 1, sum = 0;
		        for (int i=0; i<this.getDimension(); i++) {
		            sum=0;
		            for (int j=1; j<6; j++) {
		                sum = sum + j * Math.Cos((j+1) * x[i] + j);
		            }
		            result = result * sum;
		        }
		        return -result;
		    }
	    }//Shubert
	
	    /**
	     * F7: Vincent
	     * Variable range: x_i in [0.25, 10]^n, i=1,2,...,n
	     * No. of global optima: 6^n
	     * No. of local optima:  0.
	     */
	    public class Vincent : Func {

		    public Vincent(IBoundsFn boundsFn)
			    :base( boundsFn) {
		    }

		    
		    public override double doEvaluate(double[] x) {
		        double result = 0;
		        for (int i=0; i<this.getDimension(); i++){
		            if (x[i]<=0){
		        	    throw new ArgumentException();
		            }
		            result = result + Math.Sin(10 * Math.Log(x[i]) );
		        }
		        return result/(1.0*this.getDimension());
		    }
	    }//vincent
	
	    /**
	     * F8: Modified Rastrigin - All Global Optima
	     * Variable ranges: x_i in [0, 1]^n, i=1,2,...,n
	     * No. of global peaks: \prod_{i=1}^n k_i
	     * No. of local peaks:  0.
	     */
	    public class ModifiedRastriginAll : Func {

		    /* Modified Rastrigin -- All Global Optima */
		    static readonly double [] MPPF92 = {3, 4};
		    static readonly double [] MPPF98 = {1, 2, 1, 2, 1, 3, 1, 4};
            static readonly double[] MPPF916 = { 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 3, 1, 1, 1, 4 };

		    public ModifiedRastriginAll(IBoundsFn boundsFn)
			    :base( boundsFn) {
		    }

		    
		    public override double doEvaluate(double[] x) {
			     double result = 0;
			        for (int i=0; i<this.getDimension(); i++){
			            if (this.getDimension() == 2)  { 
			        	    result = result + 10+ 9*Math.Cos(2*Math.PI*MPPF92[i]*x[i]); 
			            }
			            if (this.getDimension() == 8)  { 
			        	    result = result + 10+ 9*Math.Cos(2*Math.PI*MPPF98[i]*x[i]); 
			            }
			            if (this.getDimension() == 16) { 
			        	    result = result + 10+ 9*Math.Cos(2*Math.PI*MPPF916[i]*x[i]); 
			            }
			        }
			        return -result;
		    }
		
	    }//ModifiedRastriginAll
	
	    /******************************************************************************
	     * Basic functions for composition 
	     *****************************************************************************/
	    /* Ackley's function */
	    public class FAckley : Func {

		    public FAckley(IBoundsFn boundsFn) 
			    :base( boundsFn) {
		    }

		    
		    public override double doEvaluate(double[] x) {
		        double sum1 = 0.0, sum2 =0.0, result;
		        for (int i=0; i<this.getDimension(); ++i) {
		            sum1 += x[i]*x[i];
		            sum2 += Math.Cos(2.0*Math.PI*x[i]);
		        }
		        sum1 = -0.2*Math.Sqrt(sum1/this.getDimension());
		        sum2 /= this.getDimension();
		        result = 20.0 + Math.E - 20.0*Math.Exp(sum1) - Math.Exp(sum2);
		        return result;
		    }
	    }//FAckley
	
	    /* Rastrigin's function */
	    public class FRastrigin : Func {

		    public FRastrigin(IBoundsFn boundsFn) 
			    :base( boundsFn) {
		    }

		    public override double doEvaluate(double[] x) {
			    double result = 0.0;
		        for (int i=0; i<this.getDimension(); ++i) {
		            result += (x[i]*x[i] - 10.0*Math.Cos(2.0*Math.PI*x[i]) + 10.0);
		        }
		        return result;
		    }
	    }//FRastrigin

	    /* Weierstrass's function */
	    public class FWeierstrass : Func {

		    public FWeierstrass(IBoundsFn boundsFn)
			    :base( boundsFn) {
		    }

		    public override double doEvaluate(double[] x) {
		        double result=0.0, sum=0.0, sum2=0.0, a=0.5, b=3.0;
		        int k_max=20;

		        for (int j=0; j<=k_max; ++j) {
		            sum2 += Math.Pow(a,j)*Math.Cos(2.0*Math.PI*Math.Pow(b,j)*(0.5));
		        }
		        for (int i=0; i<this.getDimension(); ++i) {
		            sum = 0.0;
		            for (int j=0; j<=k_max; ++j) {
		                sum += Math.Pow(a,j)*Math.Cos(2.0*Math.PI*Math.Pow(b,j)*(x[i]+0.5));
		            }
		            result += sum;
		        }
		        return result - sum2*this.getDimension();
		    }
	    }//FWeierstrass
	
	    /* Griewank's function */
	    public class FGriewank : Func {

		    public FGriewank(IBoundsFn boundsFn)
			    :base( boundsFn) {
		    }

		    public override double doEvaluate(double[] x) {
		        double sum=0.0, prod=1.0, result=0.0;

		        for (int i=0; i<this.getDimension(); ++i) {
		            sum  += x[i]*x[i]/4000.0;
		            prod *= Math.Cos( x[i]/Math.Sqrt((1.0+i)) );
		        }
		        result = 1.0 + sum - prod;
		        return result;
		    }
		
	    }//FGriewank

	    /* Sphere function */
	    public class FSphere : Func {

		    public FSphere(IBoundsFn boundsFn)
			    :base( boundsFn) {
			    // TODO Auto-generated constructor stub
		    }

		    public override double doEvaluate(double[] x) {
			    double result=0.0;
			    for (int i=0; i<this.getDimension(); ++i) {
				    result += x[i]*x[i];
			    }
			    return result;		
		    }
	    }//FSphere
	
	    /* Schwefel's function */
	    public class FSchwefel : Func {

		    public FSchwefel(IBoundsFn boundsFn)
			    :base( boundsFn) {
		    }

		    public override double doEvaluate(double[] x) {
		        double sum1=0.0, sum2=0.0;

		        for (int i=0; i<this.getDimension(); ++i) {
		            sum2 = 0.0;
		            for (int j=0; j<=i; ++j) {
		                sum2 += x[j];
		            }
		            sum1 += sum2*sum2;
		        }
		        return sum1;
		    }
	    }//FSchwefel

	    /* Rosenbrock's function */
	    public class FRosenbrock : Func {

		    public FRosenbrock(IBoundsFn boundsFn)
			    :base( boundsFn) {
		    }

		    public override double doEvaluate(double[] x) {
		        double result =0.0;

		        for (int i=0; i<this.getDimension()-1; ++i) {
		            result += 100.0*Math.Pow((x[i]*x[i]-x[i+1]),2.0) + 1.0*Math.Pow((x[i]-1.0),2.0);
		        }
		        return result;
		    }
	    }//FRosenbrock
	
	    /* FEF8F2 function */
	    public class FEF8F2 : Func {

		    public FEF8F2(IBoundsFn boundsFn) 
			    :base( boundsFn) {
		    }

		    
		    public override double doEvaluate(double[] xx) {
		        double result=0.0;
		        double x=0, y=0, f=0, f2=0;

		        for (int i=0; i<this.getDimension()-1; ++i) {
		            x = xx[i]   +1;
		            y = xx[i+1] +1;

		            f2 = 100.0*(x*x - y)*(x*x - y) + (1.0 - x)*(1.0 - x);
		            f  = 1.0 + f2*f2/4000.0 - Math.Cos(f2);

		            result += f;
		        }
		        /* do not forget the (dim-1,0) case! */
		        x = xx[this.getDimension()-1] +1;
		        y = xx[0]     +1;

		        f2 = 100.0*(x*x - y)*(x*x - y) + (1.0 - x)*(1.0 - x);
		        f  = 1.0 + f2*f2/4000.0 - Math.Cos(f2);

		        result += f;

		        return result;
		    }
	    }//FEF8F2
    }
}
