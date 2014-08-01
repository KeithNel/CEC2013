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
    public class ExampleUsage
    {
        static void Main(string[] args)
        {
		    /** Always Initialise a CEC2013 object **/
		    CEC2013 comp = new CEC2013();
		    /** Get the function list **/
		    List< Func > funcs = comp.getFunctions();
		
		    /*********************************************************************
		     * Example usage on how to iterate through the functions 
		     *********************************************************************/
		    int index = 1;
		    foreach ( Func f in funcs ) {

                double[] x = new double[f.getDimension()];
                for (int i = 0; i < f.getDimension(); i++)
                    x[i] = 1.0;
			
			    Console.WriteLine( "f"+index+"(1..1)= " + f.evaluate(x) );
			    index++;
		    }
		    /*********************************************************************
		     * Example usage on how to get iterate through the functions and 
		     * get the optimisation box (bounds) of each function 
		     *********************************************************************/
    //		Console.Writeline("Optimisation bounds for each function: ");
    //		List< ClosedInterval.Double > bounds;
    //		index = 1;
    //		for( Func f : funcs ) {
    //			bounds = f.getBounds();
    //			Console.Writeline("\nF"+index + ": ");
    //			for (int i=0; i<f.getDimension(); ++i) {
    //				Console.Writeline("i: "+i+" ["
    //						+ bounds.get(i).getLower() + "," 
    //						+ bounds.get(i).getUpper() + "]");
    //			}
    //			index++;
    //		}
		
		    /*********************************************************************
		     * Demonstration on using howManyGlobalOptimaInPopulation method
		     *********************************************************************/
		    /* Create and initialise a population inside bounds */
		    int NP = 30; // population size
		    List<List<double>> population;
		    List< ClosedInterval.Double > bounds;
		    double accuracy = 0.001;

            Random random = new Random();

		    /* Be careful findex has to start from value equal to 1 */
		    int findex = 1;
		    foreach ( Func f in funcs ) {
			    bounds = f.getBounds();
			    population = new List<List<double>>();

			    /* Initialise the population randomly within bounds */
			    for (int i=0; i<NP; ++i) {
				    List<double> indi = new List<double>();
				    for(int j=0; j<f.getDimension(); ++j) {
                        double tmp = bounds[j].getLower() + random.NextDouble() * 
							    (bounds[j].getUpper() - bounds[j].getLower());
					    indi.Add(tmp);
				    }
				    population.Add(indi);
			    }
			    /* Count how many global optima exist in the population */
			    List<List<double>> seeds = new List<List<double>>();
			    int numberOfOptimaFound = CEC2013.howManyGlobalOptimaInPopulation(
							    population, seeds, findex, accuracy, CEC2013.getRho(findex));

                Console.WriteLine("F{0}: {1} of {2} global optimizers found in the current population.", findex, numberOfOptimaFound, CEC2013.getNoGoptima(findex));
			
			    if (numberOfOptimaFound != 0) {
				    /* Print seeds */
				    foreach (List<double> indi in seeds){
					    for (int i=0; i<indi.Count; ++i)
                            Console.Write(indi[i] + "\t");
					    Console.WriteLine(" Fitness: "+ f.evaluate(indi));
				    }
			    }
			    findex++;
		    }
        }
    }
}
