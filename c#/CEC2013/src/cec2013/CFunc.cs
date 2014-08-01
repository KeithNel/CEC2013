/******************************************************************************
 * Version: 1.0
 * Last modified on: 30th July, 2014 
 * Developers: Michael G. Epitropakis, Xiaodong Li.
 *      email: mge_(AT)_cs_(DOT)_stir_(DOT)_ac_(DOT)_uk 
 *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
 * c# port by Keith Nelson : keith_(AT)_cmind_(DOT)_org
 * ***************************************************************************/
using System;
using System.IO;
using System.Collections.Generic;

namespace cec2013
{
    public abstract class CFunc : Func 
    {
	    protected List< Func > funcs_;
	    //protected int dimension_;
	    protected readonly int nofunc_;
	    protected double C_;// = 2000.0;
	    protected double [] lambda_;
	    protected double [] sigma_;
	    protected double [] bias_;
	    protected double [][] O_;
	    protected double [][][] M_;
	    protected double [] weight_;
	    protected double [] fi_;
	    protected double [] z_;
	    protected double f_bias_;
	    protected double [] fmaxi_;
	    protected double [] tmpx_;

	    public CFunc(IBoundsFn boundsFn, int nofunc  ) 
                    :base(boundsFn)
        {
		    this.nofunc_ 		= nofunc;

		    if (base.getDimension() < 1 || nofunc_ < 1) 
			    throw new ArgumentException();

		    initializeCFunc();
	    }

	    private void initializeCFunc() {
		    C_ 		= 2000.0;
		    f_bias_ = 0;
		    lambda_ = new double[nofunc_];
		    sigma_  = new double[nofunc_];
		    bias_   = new double[nofunc_];
		    O_  	= new double[nofunc_][]; //[base.getDimension()];
		    M_		= new double[nofunc_][][]; //[base.getDimension()][base.getDimension()];
            int boundsCount = base.getDimension();
            for (int i=0; i<nofunc_;i++)
            {
                O_[i] = new double[boundsCount];
                M_[i] = new double[boundsCount][];
                for (int j=0; j<boundsCount;j++)
                    M_[i][j] = new double[boundsCount];
            }
		    weight_ = new double[nofunc_];
		    fi_ 	= new double[nofunc_];
		    z_  	= new double[boundsCount];
		    fmaxi_  = new double[nofunc_];
		    tmpx_   = new double[boundsCount];
	    }
	
	    protected void loadOptima(string filename) 
	    {
            string[] lines;
            try {
                if (!File.Exists(filename))
                {
                    Console.WriteLine("Error: loadOptima, Can not find file: " + filename);
                    return;
                }
                lines = File.ReadAllLines(filename);
            }
            catch(Exception e) {
                Console.WriteLine("Error: loadOptima, IOException: " + e.ToString());
                return;
            }

			for (int i=0; i<this.nofunc_; ++i) {
				double tmp = 0;
				string[] numbers = lines[i].Split('\t');
				for( int j = 0; j < base.getDimension(); ++j ) {
                    if (!double.TryParse(numbers[j].Trim(), out tmp))
                    {
                        Console.WriteLine("Error: loadOptima, NumberFormatException: line {0} dimension {1}", i, j);                         
                        return;
                    }
					O_[i][j] = tmp;
				}
			}
	    }
		
	    protected void loadRotationMatrix(string filename)
	    {	
            string[] lines;
            try {
                if (!File.Exists(filename))
                {
                    Console.WriteLine("Error: loadRotationMatrix, Can not find file: " + filename);
                    return;
                }
                lines = File.ReadAllLines(filename);
            }
            catch(Exception e) {
                Console.WriteLine("Error: loadRotationMatrix, IOException: " + e.ToString());
                return;
            }

			double tmp = -1;
            int lineNumber = 0;
			for (int i=0; i<nofunc_; ++i) {
                //System.out.println("Matrix: "+i);
				for (int j=0; j<base.getDimension(); ++j) {
					String[] numbers = lines[lineNumber++].Split('\t');
				    for (int k=0; k<base.getDimension(); ++k) {
                        if (!double.TryParse(numbers[k].Trim(), out tmp))
                        {
                            Console.WriteLine("Error: loadRotationMatrix, NumberFormatException at [{0}, {1}, {2}]", i, j, k);                         
                            return;
                        }
				        M_[i][j][k] = tmp;
                        //System.out.print(M_[i][j][k]+"\t");
				    }
				}
            }
	    }

	    void calculateWeights(double [] x)
	    {
		    double sum = 0, maxi = double.MinValue, maxindex = 0;

		    for (int i=0; i<nofunc_; ++i) {
			    sum = 0.0;
			    for (int j=0; j<base.getDimension(); ++j) {
				    sum += ( x[j] - O_[i][j] ) * ( x[j] - O_[i][j] );
			    }
			    weight_[i] = Math.Exp( -sum/(2.0 * base.getDimension() * sigma_[i] * sigma_[i]) );
			    if (i==0) { maxi = weight_[i]; }
			    if (weight_[i] > maxi) {
				    maxi = weight_[i];
				    maxindex = i;
			    }
			    //maxi = max(maxi, weight_[i]);
		    }
		    sum = 0.0;
		    for (int i=0; i<nofunc_; ++i) {
			    //if (weight_[i] != maxi) {
			    if (i != maxindex) {
				    weight_[i] *= (1.0 - Math.Pow(maxi, 10.0));
			    }
			    sum += weight_[i];
		    }
		    for (int i=0; i<nofunc_; ++i) {
			    if (sum == 0.0) {
				    weight_[i] = 1.0/(double)nofunc_;
			    } else {
				    weight_[i] /= sum;
			    }
		    }
    //		for (int i=0; i<nofunc_; ++i) {
    //			System.out.print(weight_[i] +"\t");
    //		}
    //		System.out.println("");
	    }

	    protected void initRotmatIdentity()
	    {
		    for (int i=0; i<nofunc_; ++i) {
    //			System.out.println("Matrix: "+ i);
			    for (int j=0; j<base.getDimension(); ++j) {
				    for (int k=0; k<base.getDimension(); ++k) {
					    M_[i][j][k] = (j==k ? 1 : 0 );
    //					System.out.print(M_[i][j][k] + "\t");
				    }
    //				System.out.println("");
			    }
    //			System.out.println("");
		    }
	    }
	    protected void initOptimaRandomly()
	    {
		    //List< ClosedInterval.Double > bounds = base.getBounds();
            List<ClosedInterval.Double> bounds = base.getBounds();
            Random random = new Random();
		    for (int i=0; i< nofunc_; ++i) {
			    for (int j=0; j< base.getDimension(); ++j) {
				    //O_[i][j] = lbound_[j] + (ubound_[j] - lbound_[j]) * Math.random();//nextUniform(0.0,1.0);
                    ClosedInterval.Double doubleCI = bounds[j];
				    O_[i][j] = doubleCI.getLower() + 
						    (doubleCI.getUpper() - doubleCI.getLower()) * random.NextDouble();
				    //	          System.out.print(O_[i][j] +"\t");
			    }
			    //	      System.out.println("");
		    }
	    }

        void transformToZ(double[] x, int index)
	    {
		    /* Calculate z_i = (x - o_i)/\lambda_i */
		    for (int i=0; i<base.getDimension(); ++i) {
			    tmpx_[i] = (x[i] - O_[index][i])/lambda_[index];
		    }
		    /* Multiply z_i * M_i */
		    for (int i=0; i<base.getDimension(); ++i) {
			    z_[i] = 0;
			    for (int j=0; j<base.getDimension(); ++j) {
				    /* in MATLAB: M.M1*tmpx' */
				    //z_[i] += M_[index][i][j] * tmpx_[j];

				    /* in MATLAB: tmpx*M.M1 */
				    z_[i] += M_[index][j][i] * tmpx_[j];
			    }
			    //System.out.println("i: "+i+" "+tmpx_[i]+" "+x[i]+" "+O_[index][i]+" "+lambda_[index]+" "+z_[i]);
		    }
	    }

	    void transformToZNoshift(double[] x, int index)
	    {
		    /* Calculate z_i = (x - o_i)/\lambda_i */
		    for (int i=0; i<base.getDimension(); ++i) {
			    //tmpx_[i] = (x[i] - O_[index][i])/lambda_[index];
			    tmpx_[i] = (x[i])/lambda_[index];
		    }
		    /* Multiply z_i * M_i */
		    for (int i=0; i<base.getDimension(); ++i) {
			    z_[i] = 0;
			    for (int j=0; j<base.getDimension(); ++j) {
				    /* in MATLAB: M.M1*tmpx' */
				    //z_[i] += M_[index][i][j] * tmpx_[j];

				    /* in MATLAB: tmpx*M.M1 */
				    z_[i] += M_[index][j][i] * tmpx_[j];
			    }
		    }
	    }
	
	    protected void CalculateFMaxi()
	    {
		    /* functions */
		    double [] x5 = new double[base.getDimension()];
		    for (int i=0; i<base.getDimension(); ++i) { x5[i] = 5 ; }
		    int index=0;
		    foreach( Func f in funcs_ ) {
			    transformToZNoshift(x5, index);
    //			for (int i=0; i<z_.length; i++)
    //				System.out.println("z_: "+this.z_[i]);
			    fmaxi_[index] = f.cfuncEvaluate(this.z_);
    //			System.out.println("FMAXI: "+index+" : "+fmaxi_[index]);
			    index++;
		    }
	    }
	
	    protected double evaluateInner_(double[] x)
	    {
	        double result = 0;
	        calculateWeights(x);
	        int index=0;
	        foreach( Func f in funcs_ ) {
	            transformToZ(x, index);
	            fi_[index] = f.cfuncEvaluate(this.z_);
    //	        System.out.println("Func: "+i+" : "+fi_[i]);
	            index++;
	        }
	        for (int i=0; i<nofunc_; ++i) {
	            result += weight_[i]*( C_ * fi_[i] / fmaxi_[i] + bias_[i] );
	        }
	        //Assuming maximisation
	        return -1.0*result + f_bias_;
	    }

	    //public abstract double doEvaluate(double [] x);

    }
}
